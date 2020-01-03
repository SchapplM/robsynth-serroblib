% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPPRR12_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:51
% EndTime: 2019-12-31 18:06:51
% DurationCPUTime: 0.09s
% Computational Cost: add. (83->38), mult. (66->36), div. (0->0), fcn. (103->8), ass. (0->28)
t14 = sin(qJ(1));
t8 = pkin(8) + qJ(4);
t3 = cos(t8);
t32 = t14 * t3;
t16 = cos(qJ(1));
t31 = t16 * t3;
t10 = sin(pkin(8));
t30 = t14 * t10;
t13 = sin(qJ(5));
t29 = t14 * t13;
t15 = cos(qJ(5));
t28 = t14 * t15;
t27 = t16 * t13;
t26 = t16 * t15;
t9 = pkin(5) + 0;
t25 = t14 * pkin(1) + 0;
t24 = pkin(2) + t9;
t23 = t16 * pkin(1) + t14 * qJ(2) + 0;
t22 = -pkin(3) * t10 - qJ(2);
t11 = cos(pkin(8));
t21 = t11 * pkin(3) + t24;
t2 = sin(t8);
t20 = pkin(4) * t2 - pkin(7) * t3;
t12 = -pkin(6) - qJ(3);
t19 = -t14 * t12 + t25;
t18 = -t16 * qJ(2) + t25;
t17 = pkin(3) * t30 - t16 * t12 + t23;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t14, 0, 0; t14, t16, 0, 0; 0, 0, 1, t9; 0, 0, 0, 1; 0, -t16, t14, t23; 0, -t14, -t16, t18; 1, 0, 0, t9; 0, 0, 0, 1; t30, t14 * t11, t16, t16 * qJ(3) + t23; -t16 * t10, -t16 * t11, t14, t14 * qJ(3) + t18; t11, -t10, 0, t24; 0, 0, 0, 1; t14 * t2, t32, t16, t17; -t16 * t2, -t31, t14, t22 * t16 + t19; t3, -t2, 0, t21; 0, 0, 0, 1; t2 * t28 + t27, -t2 * t29 + t26, -t32, t20 * t14 + t17; -t2 * t26 + t29, t2 * t27 + t28, t31, (-t20 + t22) * t16 + t19; t3 * t15, -t3 * t13, t2, t3 * pkin(4) + t2 * pkin(7) + t21; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
