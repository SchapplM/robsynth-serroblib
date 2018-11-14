% Calculate homogenous joint transformation matrices for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% 
% Output:
% T_mdh [4x4x3]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_mdh = S3RRR1_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:42
% EndTime: 2018-11-14 10:15:42
% DurationCPUTime: 0.02s
% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (12->6), ass. (0->7)
t20 = cos(qJ(1));
t19 = cos(qJ(2));
t18 = cos(qJ(3));
t17 = sin(qJ(1));
t16 = sin(qJ(2));
t15 = sin(qJ(3));
t1 = [t20, -t17, 0, 0; t17, t20, 0, 0; 0, 0, 1, pkin(3); 0, 0, 0, 1; t19, -t16, 0, pkin(1); t16, t19, 0, 0; 0, 0, 1, pkin(4); 0, 0, 0, 1; t18, -t15, 0, pkin(2); t15, t18, 0, 0; 0, 0, 1, pkin(5); 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,3);             % numerisch
else,                         T_mdh = sym('xx', [4,4,3]); end % symbolisch

for i = 1:3
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
