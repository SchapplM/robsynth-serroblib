% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:30
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRRP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:29:54
% EndTime: 2018-11-23 18:29:54
% DurationCPUTime: 0.15s
% Computational Cost: add. (207->65), mult. (179->71), div. (0->0), fcn. (247->10), ass. (0->39)
t35 = -pkin(9) - pkin(8);
t29 = sin(qJ(3));
t48 = t29 * pkin(3);
t32 = cos(qJ(3));
t18 = t32 * pkin(3) + pkin(2);
t28 = qJ(3) + qJ(4);
t21 = qJ(5) + t28;
t16 = sin(t21);
t30 = sin(qJ(2));
t47 = t30 * t16;
t31 = sin(qJ(1));
t46 = t31 * t29;
t14 = t31 * t30;
t33 = cos(qJ(2));
t45 = t31 * t33;
t34 = cos(qJ(1));
t15 = t34 * t30;
t44 = t34 * t33;
t26 = pkin(6) + 0;
t43 = t31 * pkin(1) + 0;
t42 = t34 * pkin(1) + t31 * pkin(7) + 0;
t27 = -pkin(10) + t35;
t20 = cos(t28);
t9 = pkin(4) * t20 + t18;
t41 = t33 * t27 + t30 * t9 + t26;
t40 = pkin(2) * t33 + pkin(8) * t30;
t39 = t18 * t33 - t30 * t35;
t38 = -t34 * pkin(7) + t43;
t19 = sin(t28);
t10 = pkin(4) * t19 + t48;
t37 = t31 * t10 - t27 * t15 + t9 * t44 + t42;
t36 = -t27 * t14 + t9 * t45 + (-pkin(7) - t10) * t34 + t43;
t17 = cos(t21);
t11 = t30 * t17;
t4 = t31 * t16 + t17 * t44;
t3 = t16 * t44 - t31 * t17;
t2 = -t34 * t16 + t17 * t45;
t1 = t16 * t45 + t34 * t17;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t31, 0, 0; t31, t34, 0, 0; 0, 0, 1, t26; 0, 0, 0, 1; t44, -t15, t31, t42; t45, -t14, -t34, t38; t30, t33, 0, t26; 0, 0, 0, 1; t32 * t44 + t46, -t29 * t44 + t31 * t32, t15, t40 * t34 + t42; -t34 * t29 + t32 * t45, -t29 * t45 - t34 * t32, t14, t40 * t31 + t38; t30 * t32, -t30 * t29, -t33, t30 * pkin(2) - t33 * pkin(8) + t26; 0, 0, 0, 1; t31 * t19 + t20 * t44, -t19 * t44 + t31 * t20, t15, pkin(3) * t46 + t39 * t34 + t42; -t34 * t19 + t20 * t45, -t19 * t45 - t34 * t20, t14 (-pkin(7) - t48) * t34 + t39 * t31 + t43; t30 * t20, -t30 * t19, -t33, t30 * t18 + t33 * t35 + t26; 0, 0, 0, 1; t4, -t3, t15, t37; t2, -t1, t14, t36; t11, -t47, -t33, t41; 0, 0, 0, 1; t4, t15, t3, t4 * pkin(5) + t3 * qJ(6) + t37; t2, t14, t1, t2 * pkin(5) + t1 * qJ(6) + t36; t11, -t33, t47 (pkin(5) * t17 + qJ(6) * t16) * t30 + t41; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
