% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-29 19:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPRR11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 19:14:58
% EndTime: 2019-12-29 19:14:58
% DurationCPUTime: 0.28s
% Computational Cost: add. (89->44), mult. (119->46), div. (0->0), fcn. (167->8), ass. (0->32)
t18 = sin(qJ(1));
t17 = sin(qJ(2));
t35 = qJ(3) * t17;
t20 = cos(qJ(2));
t5 = t18 * t20;
t39 = pkin(2) * t5 + t18 * t35;
t16 = sin(qJ(4));
t38 = t17 * t16;
t37 = t18 * t17;
t21 = cos(qJ(1));
t36 = t21 * t17;
t6 = t21 * t20;
t14 = pkin(5) + 0;
t34 = t18 * pkin(1) + 0;
t33 = t17 * pkin(2) + t14;
t32 = t21 * pkin(1) + t18 * pkin(6) + 0;
t31 = t34 + t39;
t15 = qJ(4) + qJ(5);
t8 = sin(t15);
t9 = cos(t15);
t30 = t17 * t9 - t20 * t8;
t29 = t17 * t8 + t20 * t9;
t19 = cos(qJ(4));
t28 = -t20 * t16 + t17 * t19;
t27 = t20 * t19 + t38;
t26 = -t21 * pkin(6) + t34;
t25 = pkin(2) * t6 + t21 * t35 + t32;
t7 = t19 * pkin(4) + pkin(3);
t24 = pkin(4) * t38 + t20 * t7;
t23 = -t20 * qJ(3) + t33;
t22 = -pkin(8) - pkin(7);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t18, 0, 0; t18, t21, 0, 0; 0, 0, 1, t14; 0, 0, 0, 1; t6, -t36, t18, t32; t5, -t37, -t21, t26; t17, t20, 0, t14; 0, 0, 0, 1; t6, t18, t36, t25; t5, -t21, t37, t26 + t39; t17, 0, -t20, t23; 0, 0, 0, 1; t27 * t21, t28 * t21, -t18, pkin(3) * t6 - t18 * pkin(7) + t25; t27 * t18, t28 * t18, t21, pkin(3) * t5 + (-pkin(6) + pkin(7)) * t21 + t31; t28, -t27, 0, t17 * pkin(3) + t23; 0, 0, 0, 1; t29 * t21, t30 * t21, -t18, t18 * t22 + t21 * t24 + t25; t29 * t18, t30 * t18, t21, (-pkin(6) - t22) * t21 + t24 * t18 + t31; t30, -t29, 0, t17 * t7 + (-pkin(4) * t16 - qJ(3)) * t20 + t33; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
