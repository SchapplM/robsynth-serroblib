% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2018-11-23 15:17
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRRPRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:17:03
% EndTime: 2018-11-23 15:17:03
% DurationCPUTime: 0.24s
% Computational Cost: add. (609->86), mult. (635->104), div. (0->0), fcn. (717->18), ass. (0->56)
t37 = sin(pkin(12));
t61 = t37 * pkin(4) + pkin(8);
t36 = pkin(12) + qJ(5);
t28 = sin(t36);
t69 = pkin(5) * t28 + t61;
t40 = cos(pkin(12));
t27 = t40 * pkin(4) + pkin(3);
t68 = cos(qJ(3));
t38 = sin(pkin(11));
t39 = sin(pkin(6));
t67 = t38 * t39;
t41 = cos(pkin(11));
t66 = t41 * t39;
t42 = -pkin(9) - qJ(4);
t65 = cos(pkin(6));
t64 = pkin(6) - qJ(2);
t63 = pkin(6) + qJ(2);
t62 = qJ(1) + 0;
t60 = t39 * t68;
t59 = t41 * pkin(1) + pkin(7) * t67 + 0;
t58 = t65 * pkin(7) + t62;
t57 = cos(t63);
t56 = sin(t64);
t52 = sin(t63) / 0.2e1;
t17 = t52 - t56 / 0.2e1;
t45 = cos(qJ(2));
t10 = -t38 * t17 + t41 * t45;
t55 = t10 * pkin(2) + t59;
t53 = cos(t64) / 0.2e1;
t18 = t53 - t57 / 0.2e1;
t54 = t18 * pkin(2) + t58;
t51 = t38 * pkin(1) - pkin(7) * t66 + 0;
t44 = sin(qJ(2));
t46 = t53 + t57 / 0.2e1;
t9 = t38 * t46 + t41 * t44;
t50 = t9 * pkin(8) + t55;
t8 = t41 * t17 + t38 * t45;
t49 = t8 * pkin(2) + t51;
t16 = t52 + t56 / 0.2e1;
t48 = -t16 * pkin(8) + t54;
t7 = t38 * t44 - t41 * t46;
t47 = t7 * pkin(8) + t49;
t43 = sin(qJ(3));
t35 = -pkin(10) + t42;
t30 = qJ(6) + t36;
t29 = cos(t36);
t26 = cos(t30);
t25 = sin(t30);
t15 = pkin(5) * t29 + t27;
t12 = t18 * t68 + t65 * t43;
t11 = t18 * t43 - t65 * t68;
t4 = t10 * t68 + t43 * t67;
t3 = t10 * t43 - t38 * t60;
t2 = -t43 * t66 + t8 * t68;
t1 = t41 * t60 + t8 * t43;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t41, -t38, 0, 0; t38, t41, 0, 0; 0, 0, 1, t62; 0, 0, 0, 1; t10, -t9, t67, t59; t8, -t7, -t66, t51; t18, t16, t65, t58; 0, 0, 0, 1; t4, -t3, t9, t50; t2, -t1, t7, t47; t12, -t11, -t16, t48; 0, 0, 0, 1; t9 * t37 + t4 * t40, -t4 * t37 + t9 * t40, t3, t4 * pkin(3) + t3 * qJ(4) + t50; t2 * t40 + t7 * t37, -t2 * t37 + t7 * t40, t1, t2 * pkin(3) + t1 * qJ(4) + t47; t12 * t40 - t16 * t37, -t12 * t37 - t16 * t40, t11, t12 * pkin(3) + t11 * qJ(4) + t48; 0, 0, 0, 1; t9 * t28 + t4 * t29, -t4 * t28 + t9 * t29, t3, t4 * t27 - t3 * t42 + t61 * t9 + t55; t2 * t29 + t7 * t28, -t2 * t28 + t7 * t29, t1, -t1 * t42 + t2 * t27 + t61 * t7 + t49; t12 * t29 - t16 * t28, -t12 * t28 - t16 * t29, t11, -t11 * t42 + t12 * t27 - t61 * t16 + t54; 0, 0, 0, 1; t9 * t25 + t4 * t26, -t4 * t25 + t9 * t26, t3, t4 * t15 - t3 * t35 + t69 * t9 + t55; t2 * t26 + t7 * t25, -t2 * t25 + t7 * t26, t1, -t1 * t35 + t2 * t15 + t69 * t7 + t49; t12 * t26 - t16 * t25, -t12 * t25 - t16 * t26, t11, -t11 * t35 + t12 * t15 - t69 * t16 + t54; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
