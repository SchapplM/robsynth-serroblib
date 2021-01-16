% Calculate minimal parameter regressor of gravitation load for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:06:04
% EndTime: 2021-01-16 02:06:06
% DurationCPUTime: 0.52s
% Computational Cost: add. (372->96), mult. (635->160), div. (0->0), fcn. (767->14), ass. (0->58)
t38 = sin(qJ(3));
t40 = cos(qJ(3));
t61 = cos(pkin(6));
t35 = sin(pkin(6));
t39 = sin(qJ(2));
t65 = t35 * t39;
t77 = -t38 * t65 + t61 * t40;
t34 = sin(pkin(10));
t54 = t34 * t61;
t60 = cos(pkin(10));
t72 = cos(qJ(2));
t13 = t39 * t54 - t60 * t72;
t64 = t35 * t40;
t76 = t13 * t38 + t34 * t64;
t32 = qJ(3) + pkin(11);
t28 = sin(t32);
t30 = cos(t32);
t10 = t28 * t65 - t30 * t61;
t47 = t61 * t60;
t15 = t34 * t72 + t39 * t47;
t53 = t35 * t60;
t4 = t15 * t28 + t30 * t53;
t66 = t34 * t35;
t6 = t13 * t28 + t30 * t66;
t46 = g(1) * t6 - g(2) * t4 - g(3) * t10;
t74 = g(3) * t35;
t16 = t39 * t60 + t54 * t72;
t26 = pkin(3) * t40 + pkin(2);
t37 = qJ(4) + pkin(8);
t73 = -t13 * t37 - t16 * t26;
t31 = pkin(12) + qJ(6);
t27 = sin(t31);
t70 = t27 * t30;
t29 = cos(t31);
t69 = t29 * t30;
t33 = sin(pkin(12));
t68 = t30 * t33;
t36 = cos(pkin(12));
t67 = t30 * t36;
t56 = t35 * t72;
t63 = t26 * t56 + t37 * t65;
t62 = qJ(5) * t28;
t57 = t30 * t72;
t14 = t34 * t39 - t47 * t72;
t55 = -t14 * t26 + t15 * t37;
t3 = t13 * t30 - t28 * t66;
t51 = t76 * pkin(3);
t49 = -pkin(4) * t30 - t62;
t48 = t77 * pkin(3);
t11 = t28 * t61 + t30 * t65;
t5 = t15 * t30 - t28 * t53;
t45 = g(1) * t3 - g(2) * t5 - g(3) * t11;
t44 = -t15 * t38 - t40 * t53;
t43 = -g(1) * t13 + g(2) * t15 + g(3) * t65;
t42 = t44 * pkin(3);
t41 = g(1) * t16 + g(2) * t14 - g(3) * t56;
t1 = t41 * t28;
t2 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t41, t43, 0, 0, 0, 0, 0, t41 * t40, -t41 * t38, t41 * t30, -t1, -t43, -g(1) * t73 - g(2) * t55 - g(3) * t63, -g(1) * (-t13 * t33 - t16 * t67) - g(2) * (-t14 * t67 + t15 * t33) - (t33 * t39 + t36 * t57) * t74, -g(1) * (-t13 * t36 + t16 * t68) - g(2) * (t14 * t68 + t15 * t36) - (-t33 * t57 + t36 * t39) * t74, t1, -g(1) * (t16 * t49 + t73) - g(2) * (t14 * t49 + t55) - g(3) * ((pkin(4) * t57 + t62 * t72) * t35 + t63), 0, 0, 0, 0, 0, -g(1) * (-t13 * t27 - t16 * t69) - g(2) * (-t14 * t69 + t15 * t27) - (t27 * t39 + t29 * t57) * t74, -g(1) * (-t13 * t29 + t16 * t70) - g(2) * (t14 * t70 + t15 * t29) - (-t27 * t57 + t29 * t39) * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t76 - g(2) * t44 - g(3) * t77, -g(1) * (t13 * t40 - t38 * t66) - g(2) * (-t15 * t40 + t38 * t53) - g(3) * (-t38 * t61 - t39 * t64), -t46, -t45, 0, -g(1) * t51 - g(2) * t42 - g(3) * t48, -t46 * t36, t46 * t33, t45, -g(1) * (pkin(4) * t6 - qJ(5) * t3 + t51) - g(2) * (-t4 * pkin(4) + t5 * qJ(5) + t42) - g(3) * (-pkin(4) * t10 + qJ(5) * t11 + t48), 0, 0, 0, 0, 0, -t46 * t29, t46 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, 0, 0, -t41, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t16 * t29 + t27 * t3) - g(2) * (t14 * t29 - t27 * t5) - g(3) * (-t11 * t27 - t29 * t56), -g(1) * (-t16 * t27 + t29 * t3) - g(2) * (-t14 * t27 - t29 * t5) - g(3) * (-t11 * t29 + t27 * t56);];
taug_reg = t2;
