% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-29 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPPRR10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 16:22:49
% EndTime: 2019-12-29 16:22:49
% DurationCPUTime: 0.30s
% Computational Cost: add. (89->44), mult. (119->46), div. (0->0), fcn. (167->8), ass. (0->32)
t19 = sin(qJ(1));
t16 = sin(pkin(8));
t35 = qJ(3) * t16;
t17 = cos(pkin(8));
t5 = t19 * t17;
t39 = pkin(2) * t5 + t19 * t35;
t18 = sin(qJ(4));
t38 = t16 * t18;
t37 = t19 * t16;
t21 = cos(qJ(1));
t36 = t21 * t16;
t6 = t21 * t17;
t14 = pkin(5) + 0;
t34 = t19 * pkin(1) + 0;
t33 = t16 * pkin(2) + t14;
t32 = t21 * pkin(1) + t19 * qJ(2) + 0;
t31 = t34 + t39;
t15 = qJ(4) + qJ(5);
t8 = sin(t15);
t9 = cos(t15);
t30 = t16 * t9 - t17 * t8;
t29 = t16 * t8 + t17 * t9;
t20 = cos(qJ(4));
t28 = t16 * t20 - t17 * t18;
t27 = t17 * t20 + t38;
t26 = pkin(2) * t6 + t21 * t35 + t32;
t25 = -t21 * qJ(2) + t34;
t7 = t20 * pkin(4) + pkin(3);
t24 = pkin(4) * t38 + t17 * t7;
t23 = -t17 * qJ(3) + t33;
t22 = -pkin(7) - pkin(6);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t19, 0, 0; t19, t21, 0, 0; 0, 0, 1, t14; 0, 0, 0, 1; t6, -t36, t19, t32; t5, -t37, -t21, t25; t16, t17, 0, t14; 0, 0, 0, 1; t6, t19, t36, t26; t5, -t21, t37, t25 + t39; t16, 0, -t17, t23; 0, 0, 0, 1; t27 * t21, t28 * t21, -t19, pkin(3) * t6 - t19 * pkin(6) + t26; t27 * t19, t28 * t19, t21, pkin(3) * t5 + (pkin(6) - qJ(2)) * t21 + t31; t28, -t27, 0, t16 * pkin(3) + t23; 0, 0, 0, 1; t29 * t21, t30 * t21, -t19, t19 * t22 + t24 * t21 + t26; t29 * t19, t30 * t19, t21, (-qJ(2) - t22) * t21 + t24 * t19 + t31; t30, -t29, 0, t16 * t7 + (-pkin(4) * t18 - qJ(3)) * t17 + t33; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
