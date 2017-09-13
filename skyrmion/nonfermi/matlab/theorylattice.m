clear;
data = textread('./data/lattice_4.out');
L = 60;
num = data(0*(L*L)+1:0*(L*L)+((L*L)),1);
theta = data(0*(L*L)+1:0*(L*L)+((L*L)),2);
fai = data(0*(L*L)+1:0*(L*L)+((L*L)),3);
N = L*L;
% ep = zeros(L,L,3,3);
eigenvalue = zeros(N,3);
eigenvector = zeros(N,3,3);
A = zeros(3,3);
out_matrix = zeros(N,5);
for ii = 1:L
    for jj = 1:L
        for  m= 1:N
            for n = 1:N
                i_para_m = rem( m-1,L );
                j_para_m = fix( (m-1)/L );
                i_para_n = rem( n-1,L );
                j_para_n = fix( (n-1)/L );
                S_m(1) = sin( theta(m) )*cos( fai(m) );
                S_m(2) = sin( theta(m) )*sin( fai(m) );
                S_m(3) = cos( theta(m) );
                S_n(1) = sin( theta(n) )*cos( fai(n) );
                S_n(2) = sin( theta(n) )*sin( fai(n) );
                S_n(3) = cos( theta(n) );
                for mm = 1:3
                    for nn = 1:3
                        % ep(ii,jj,mm,nn) = ep(ii,jj,mm,nn)+1/N^2*exp( sqrt(-1)*( ii*(i_para_m-i_para_n)+jj*(j_para_m-j_para_n) ) )*S_m(mm)*S_n(nn);
                        A(mm,nn) =  1/N^2*exp( sqrt(-1)*2*pi*( ii*(i_para_m-i_para_n)+jj*(j_para_m-j_para_n) ) )*S_m(mm)*S_n(nn);
                    end
                end
                [V, D] = eig( A );
                % disp('run')
                eigenvalue( (ii-1)*L+jj,1 ) = D(1,1);
                eigenvalue( (ii-1)*L+jj,2 ) = D(2,2);
                eigenvalue( (ii-1)*L+jj,3 ) = D(3,3);
                for i_i = 1:3
                    for j_j = 1:3
                        eigenvector( (ii-1)*L+jj, i_i, j_j ) = V(i_i,j_j);
                    end
                end
            end
        end
    end
end
disp('>>Inint success!')
min_value = min(min(eigenvalue));
[x1,y1] = find(eigenvalue == max(max(eigenvalue)));
temp = eigenvalue;
temp(x1,y1) = min_value;
[x2,y2] = find(eigenvalue == max(max(temp)));
temp1 = temp;
temp1(x2,y2) = min_value;
[x3,y3] = find(eigenvalue == max(max(temp1)));
for ii = 1:3
    S1(ii) = eigenvector(x1, ii, y1);
    S2(ii) = eigenvector(x2, ii, y2);
    S3(ii) = eigenvector(x3, ii, y3);
end
disp('>>caculate S(x) success!')
lam = 55.0/89;
for ii = 0:L-1
    for jj = 0:L-1
        x( ii*L+jj+1 ) = ii+jj/2;
        y( ii*L+jj+1 ) = sqrt(3.0)/2*jj;
        aa = ( S1(1)*exp( sqrt(-1)*lam*(ii+jj/2) )+S2(1)*exp( sqrt(-1)*lam*(ii/2+jj) )+S3(1)*exp( sqrt(-1)*lam*(-ii/2+jj) ) );
        Sx( ii*L+jj+1 ) = 1/N*( aa + conj(aa) );
        aa = ( S1(2)*exp( sqrt(-1)*lam*(ii+jj/2) )+S2(2)*exp( sqrt(-1)*lam*(ii/2+jj) )+S3(2)*exp( sqrt(-1)*lam*(-ii/2+jj) ) );
        Sy( ii*L+jj+1 ) = 1/N*( aa + conj(aa) );
        aa = ( S1(3)*exp( sqrt(-1)*lam*(ii+jj/2) )+S2(3)*exp( sqrt(-1)*lam*(ii/2+jj) )+S3(3)*exp( sqrt(-1)*lam*(-ii/2+jj) ) );
        Sz( ii*L+jj+1 ) = 1/N*( aa + conj(aa) );
        
%         x_new(ii*L+jj+1) = x(ii*L+jj+1)-Sx( ii*L+jj+1 )/2;
%         y_new(ii*L+jj+1) = y(ii*L+jj+1)-Sy( ii*L+jj+1 )/2;
%         z_new(ii*L+jj+1) = 0-Sz( ii*L+jj+1 )/2;

    end
end
% for ii = 1:L
%     for jj = 1:L
%         map_x(ii,jj) = x((ii-1)*L+jj);
%         map_y(ii,jj) = y((ii-1)*L+jj);
%         map(ii,jj) = Sz((ii-1)*L+jj);

%     end
% end
for ii = 1:N
    out_matrix(ii,1) = x(ii);
    out_matrix(ii,2) = y(ii);
    out_matrix(ii,3) = Sx(ii);
    out_matrix(ii,4) = Sy(ii);
    out_matrix(ii,5) = Sz(ii);
end
dlmwrite('out.dat', out_matrix,'\t');
