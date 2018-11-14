% Calculate homogenous joint transformation matrices for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% 
% Output:
% T_mdh [4x4x3]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_mdh = S3RRP1_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:05
% EndTime: 2018-11-14 10:15:05
% DurationCPUTime: 0.02s
% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (8->4), ass. (0->5)
t15 = cos(qJ(1));
t14 = cos(qJ(2));
t13 = sin(qJ(1));
t12 = sin(qJ(2));
t1 = [t15, -t13, 0, 0; t13, t15, 0, 0; 0, 0, 1, pkin(3); 0, 0, 0, 1; t14, -t12, 0, pkin(1); t12, t14, 0, 0; 0, 0, 1, pkin(4); 0, 0, 0, 1; 1, 0, 0, pkin(2); 0, 0, -1, -qJ(3); 0, 1, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,3);             % numerisch
else,                         T_mdh = sym('xx', [4,4,3]); end % symbolisch

for i = 1:3
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
