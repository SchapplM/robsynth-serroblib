% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2018-11-23 15:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRRRPR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:24:10
% EndTime: 2018-11-23 15:24:10
% DurationCPUTime: 0.23s
% Computational Cost: add. (609->86), mult. (635->104), div. (0->0), fcn. (717->18), ass. (0->56)
t41 = sin(qJ(4));
t61 = t41 * pkin(4) + pkin(8);
t36 = qJ(4) + pkin(12);
t28 = sin(t36);
t69 = pkin(5) * t28 + t61;
t44 = cos(qJ(4));
t27 = t44 * pkin(4) + pkin(3);
t68 = cos(qJ(3));
t37 = sin(pkin(11));
t38 = sin(pkin(6));
t67 = t37 * t38;
t39 = cos(pkin(11));
t66 = t39 * t38;
t40 = -qJ(5) - pkin(9);
t65 = cos(pkin(6));
t64 = pkin(6) - qJ(2);
t63 = pkin(6) + qJ(2);
t62 = qJ(1) + 0;
t60 = t38 * t68;
t59 = t39 * pkin(1) + pkin(7) * t67 + 0;
t58 = t65 * pkin(7) + t62;
t57 = cos(t63);
t56 = sin(t64);
t52 = sin(t63) / 0.2e1;
t17 = t52 - t56 / 0.2e1;
t45 = cos(qJ(2));
t10 = -t37 * t17 + t39 * t45;
t55 = t10 * pkin(2) + t59;
t53 = cos(t64) / 0.2e1;
t18 = t53 - t57 / 0.2e1;
t54 = t18 * pkin(2) + t58;
t51 = t37 * pkin(1) - pkin(7) * t66 + 0;
t43 = sin(qJ(2));
t46 = t53 + t57 / 0.2e1;
t9 = t37 * t46 + t39 * t43;
t50 = t9 * pkin(8) + t55;
t8 = t39 * t17 + t37 * t45;
t49 = t8 * pkin(2) + t51;
t16 = t52 + t56 / 0.2e1;
t48 = -t16 * pkin(8) + t54;
t7 = t37 * t43 - t39 * t46;
t47 = t7 * pkin(8) + t49;
t42 = sin(qJ(3));
t35 = -pkin(10) + t40;
t30 = qJ(6) + t36;
t29 = cos(t36);
t26 = cos(t30);
t25 = sin(t30);
t15 = pkin(5) * t29 + t27;
t12 = t18 * t68 + t65 * t42;
t11 = t18 * t42 - t65 * t68;
t4 = t10 * t68 + t42 * t67;
t3 = t10 * t42 - t37 * t60;
t2 = -t42 * t66 + t8 * t68;
t1 = t39 * t60 + t8 * t42;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t39, -t37, 0, 0; t37, t39, 0, 0; 0, 0, 1, t62; 0, 0, 0, 1; t10, -t9, t67, t59; t8, -t7, -t66, t51; t18, t16, t65, t58; 0, 0, 0, 1; t4, -t3, t9, t50; t2, -t1, t7, t47; t12, -t11, -t16, t48; 0, 0, 0, 1; t4 * t44 + t9 * t41, -t4 * t41 + t9 * t44, t3, t4 * pkin(3) + t3 * pkin(9) + t50; t2 * t44 + t7 * t41, -t2 * t41 + t7 * t44, t1, t2 * pkin(3) + t1 * pkin(9) + t47; t12 * t44 - t16 * t41, -t12 * t41 - t16 * t44, t11, t12 * pkin(3) + t11 * pkin(9) + t48; 0, 0, 0, 1; t9 * t28 + t4 * t29, -t4 * t28 + t9 * t29, t3, t4 * t27 - t3 * t40 + t61 * t9 + t55; t2 * t29 + t7 * t28, -t2 * t28 + t7 * t29, t1, -t1 * t40 + t2 * t27 + t61 * t7 + t49; t12 * t29 - t16 * t28, -t12 * t28 - t16 * t29, t11, -t11 * t40 + t12 * t27 - t61 * t16 + t54; 0, 0, 0, 1; t9 * t25 + t4 * t26, -t4 * t25 + t9 * t26, t3, t4 * t15 - t3 * t35 + t69 * t9 + t55; t2 * t26 + t7 * t25, -t2 * t25 + t7 * t26, t1, -t1 * t35 + t2 * t15 + t69 * t7 + t49; t12 * t26 - t16 * t25, -t12 * t25 - t16 * t26, t11, -t11 * t35 + t12 * t15 - t69 * t16 + t54; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
