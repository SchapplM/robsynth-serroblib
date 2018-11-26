% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2018-11-23 17:16
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRRP8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:15:31
% EndTime: 2018-11-23 17:15:31
% DurationCPUTime: 0.16s
% Computational Cost: add. (207->65), mult. (179->71), div. (0->0), fcn. (247->10), ass. (0->39)
t29 = sin(pkin(10));
t48 = t29 * pkin(3);
t30 = cos(pkin(10));
t16 = t30 * pkin(3) + pkin(2);
t27 = pkin(10) + qJ(4);
t21 = qJ(5) + t27;
t14 = sin(t21);
t32 = sin(qJ(2));
t47 = t32 * t14;
t33 = sin(qJ(1));
t46 = t33 * t29;
t17 = t33 * t32;
t34 = cos(qJ(2));
t45 = t33 * t34;
t35 = cos(qJ(1));
t18 = t35 * t32;
t44 = t35 * t34;
t31 = -pkin(8) - qJ(3);
t28 = pkin(6) + 0;
t43 = t33 * pkin(1) + 0;
t42 = t35 * pkin(1) + t33 * pkin(7) + 0;
t26 = -pkin(9) + t31;
t20 = cos(t27);
t9 = pkin(4) * t20 + t16;
t41 = t34 * t26 + t32 * t9 + t28;
t40 = pkin(2) * t34 + qJ(3) * t32;
t39 = t16 * t34 - t31 * t32;
t38 = -t35 * pkin(7) + t43;
t19 = sin(t27);
t10 = pkin(4) * t19 + t48;
t37 = t33 * t10 - t26 * t18 + t9 * t44 + t42;
t36 = -t26 * t17 + t9 * t45 + (-pkin(7) - t10) * t35 + t43;
t15 = cos(t21);
t11 = t32 * t15;
t4 = t33 * t14 + t15 * t44;
t3 = t14 * t44 - t33 * t15;
t2 = -t35 * t14 + t15 * t45;
t1 = t14 * t45 + t35 * t15;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t35, -t33, 0, 0; t33, t35, 0, 0; 0, 0, 1, t28; 0, 0, 0, 1; t44, -t18, t33, t42; t45, -t17, -t35, t38; t32, t34, 0, t28; 0, 0, 0, 1; t30 * t44 + t46, -t29 * t44 + t33 * t30, t18, t40 * t35 + t42; -t35 * t29 + t30 * t45, -t29 * t45 - t35 * t30, t17, t40 * t33 + t38; t32 * t30, -t32 * t29, -t34, t32 * pkin(2) - t34 * qJ(3) + t28; 0, 0, 0, 1; t33 * t19 + t20 * t44, -t19 * t44 + t33 * t20, t18, pkin(3) * t46 + t39 * t35 + t42; -t35 * t19 + t20 * t45, -t19 * t45 - t35 * t20, t17 (-pkin(7) - t48) * t35 + t39 * t33 + t43; t32 * t20, -t32 * t19, -t34, t32 * t16 + t34 * t31 + t28; 0, 0, 0, 1; t4, -t3, t18, t37; t2, -t1, t17, t36; t11, -t47, -t34, t41; 0, 0, 0, 1; t4, t18, t3, t4 * pkin(5) + t3 * qJ(6) + t37; t2, t17, t1, t2 * pkin(5) + t1 * qJ(6) + t36; t11, -t34, t47 (pkin(5) * t15 + qJ(6) * t14) * t32 + t41; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
