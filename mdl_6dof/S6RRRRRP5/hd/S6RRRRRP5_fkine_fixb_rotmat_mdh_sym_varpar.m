% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRP5
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
% Datum: 2018-11-23 18:29
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRRP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:29:04
% EndTime: 2018-11-23 18:29:04
% DurationCPUTime: 0.15s
% Computational Cost: add. (196->68), mult. (164->74), div. (0->0), fcn. (228->10), ass. (0->41)
t35 = -pkin(9) - pkin(8);
t29 = sin(qJ(3));
t21 = t29 * pkin(3);
t32 = cos(qJ(3));
t16 = t32 * pkin(3) + pkin(2);
t28 = qJ(3) + qJ(4);
t20 = qJ(5) + t28;
t14 = sin(t20);
t30 = sin(qJ(2));
t46 = t30 * t14;
t31 = sin(qJ(1));
t45 = t31 * t29;
t33 = cos(qJ(2));
t44 = t31 * t33;
t34 = cos(qJ(1));
t43 = t34 * t33;
t17 = sin(t28);
t8 = pkin(4) * t17 + t21;
t27 = -pkin(10) + t35;
t26 = pkin(6) + 0;
t42 = t31 * pkin(1) + 0;
t18 = cos(t28);
t7 = pkin(4) * t18 + t16;
t41 = t34 * pkin(1) + t31 * pkin(7) + 0;
t40 = pkin(2) * t33 + pkin(8) * t30;
t19 = -qJ(6) + t27;
t15 = cos(t20);
t5 = pkin(5) * t15 + t7;
t39 = -t19 * t30 + t33 * t5;
t38 = -t27 * t30 + t33 * t7;
t37 = t16 * t33 - t30 * t35;
t36 = -t34 * pkin(7) + t42;
t13 = t34 * t30;
t12 = t31 * t30;
t9 = t30 * t15;
t6 = pkin(5) * t14 + t8;
t4 = t31 * t14 + t15 * t43;
t3 = -t14 * t43 + t31 * t15;
t2 = -t34 * t14 + t15 * t44;
t1 = -t14 * t44 - t34 * t15;
t10 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t31, 0, 0; t31, t34, 0, 0; 0, 0, 1, t26; 0, 0, 0, 1; t43, -t13, t31, t41; t44, -t12, -t34, t36; t30, t33, 0, t26; 0, 0, 0, 1; t32 * t43 + t45, -t29 * t43 + t31 * t32, t13, t40 * t34 + t41; -t34 * t29 + t32 * t44, -t29 * t44 - t34 * t32, t12, t40 * t31 + t36; t30 * t32, -t30 * t29, -t33, t30 * pkin(2) - t33 * pkin(8) + t26; 0, 0, 0, 1; t31 * t17 + t18 * t43, -t17 * t43 + t31 * t18, t13, pkin(3) * t45 + t37 * t34 + t41; -t34 * t17 + t18 * t44, -t17 * t44 - t34 * t18, t12 (-pkin(7) - t21) * t34 + t37 * t31 + t42; t30 * t18, -t30 * t17, -t33, t30 * t16 + t33 * t35 + t26; 0, 0, 0, 1; t4, t3, t13, t31 * t8 + t38 * t34 + t41; t2, t1, t12 (-pkin(7) - t8) * t34 + t38 * t31 + t42; t9, -t46, -t33, t33 * t27 + t30 * t7 + t26; 0, 0, 0, 1; t4, t3, t13, t31 * t6 + t39 * t34 + t41; t2, t1, t12 (-pkin(7) - t6) * t34 + t39 * t31 + t42; t9, -t46, -t33, t33 * t19 + t30 * t5 + t26; 0, 0, 0, 1;];
T_ges = t10;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
