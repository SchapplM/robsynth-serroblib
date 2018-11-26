% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 18:39
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:39:27
% EndTime: 2018-11-23 18:39:27
% DurationCPUTime: 0.16s
% Computational Cost: add. (206->73), mult. (164->84), div. (0->0), fcn. (228->12), ass. (0->38)
t33 = -pkin(9) - pkin(8);
t27 = sin(qJ(3));
t19 = t27 * pkin(3);
t30 = cos(qJ(3));
t13 = t30 * pkin(3) + pkin(2);
t26 = qJ(3) + qJ(4);
t15 = sin(t26);
t4 = pkin(4) * t15 + t19;
t29 = sin(qJ(1));
t43 = t29 * t27;
t31 = cos(qJ(2));
t42 = t29 * t31;
t32 = cos(qJ(1));
t41 = t32 * t31;
t25 = -pkin(10) + t33;
t24 = pkin(6) + 0;
t16 = cos(t26);
t3 = pkin(4) * t16 + t13;
t40 = t29 * pkin(1) + 0;
t17 = qJ(5) + t26;
t39 = t32 * pkin(1) + t29 * pkin(7) + 0;
t28 = sin(qJ(2));
t38 = pkin(2) * t31 + pkin(8) * t28;
t12 = cos(t17);
t1 = pkin(5) * t12 + t3;
t18 = -pkin(11) + t25;
t37 = t1 * t31 - t18 * t28;
t36 = -t25 * t28 + t3 * t31;
t35 = t13 * t31 - t28 * t33;
t34 = -t32 * pkin(7) + t40;
t14 = qJ(6) + t17;
t11 = sin(t17);
t10 = t32 * t28;
t9 = t29 * t28;
t6 = cos(t14);
t5 = sin(t14);
t2 = pkin(5) * t11 + t4;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t32, -t29, 0, 0; t29, t32, 0, 0; 0, 0, 1, t24; 0, 0, 0, 1; t41, -t10, t29, t39; t42, -t9, -t32, t34; t28, t31, 0, t24; 0, 0, 0, 1; t30 * t41 + t43, -t27 * t41 + t29 * t30, t10, t38 * t32 + t39; -t32 * t27 + t30 * t42, -t27 * t42 - t32 * t30, t9, t38 * t29 + t34; t28 * t30, -t28 * t27, -t31, t28 * pkin(2) - t31 * pkin(8) + t24; 0, 0, 0, 1; t29 * t15 + t16 * t41, -t15 * t41 + t29 * t16, t10, pkin(3) * t43 + t35 * t32 + t39; -t32 * t15 + t16 * t42, -t15 * t42 - t32 * t16, t9 (-pkin(7) - t19) * t32 + t35 * t29 + t40; t28 * t16, -t28 * t15, -t31, t28 * t13 + t31 * t33 + t24; 0, 0, 0, 1; t29 * t11 + t12 * t41, -t11 * t41 + t29 * t12, t10, t29 * t4 + t36 * t32 + t39; -t32 * t11 + t12 * t42, -t11 * t42 - t32 * t12, t9 (-pkin(7) - t4) * t32 + t36 * t29 + t40; t28 * t12, -t28 * t11, -t31, t31 * t25 + t28 * t3 + t24; 0, 0, 0, 1; t29 * t5 + t6 * t41, t29 * t6 - t5 * t41, t10, t29 * t2 + t37 * t32 + t39; -t32 * t5 + t6 * t42, -t32 * t6 - t5 * t42, t9 (-pkin(7) - t2) * t32 + t37 * t29 + t40; t28 * t6, -t28 * t5, -t31, t28 * t1 + t31 * t18 + t24; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
