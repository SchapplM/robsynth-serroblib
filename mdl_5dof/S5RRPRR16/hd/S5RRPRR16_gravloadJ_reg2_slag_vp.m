% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR16_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:47:26
% EndTime: 2019-12-31 20:47:27
% DurationCPUTime: 0.42s
% Computational Cost: add. (295->99), mult. (757->155), div. (0->0), fcn. (910->10), ass. (0->62)
t41 = sin(qJ(2));
t42 = sin(qJ(1));
t45 = cos(qJ(2));
t46 = cos(qJ(1));
t67 = cos(pkin(5));
t63 = t46 * t67;
t23 = t41 * t63 + t42 * t45;
t64 = t42 * t67;
t25 = -t41 * t64 + t46 * t45;
t87 = -g(1) * t25 - g(2) * t23;
t39 = sin(qJ(5));
t43 = cos(qJ(5));
t22 = t42 * t41 - t45 * t63;
t40 = sin(qJ(4));
t44 = cos(qJ(4));
t38 = sin(pkin(5));
t74 = t38 * t46;
t54 = -t22 * t40 + t44 * t74;
t86 = t23 * t43 + t39 * t54;
t85 = -t23 * t39 + t43 * t54;
t82 = g(3) * t38;
t81 = t22 * pkin(8);
t24 = t46 * t41 + t45 * t64;
t80 = t24 * pkin(8);
t77 = t38 * t41;
t76 = t38 * t42;
t75 = t38 * t45;
t73 = t39 * t40;
t72 = t39 * t41;
t71 = t40 * t43;
t70 = t41 * t43;
t69 = pkin(2) * t75 + qJ(3) * t77;
t68 = t46 * pkin(1) + pkin(7) * t76;
t66 = pkin(8) * t75 + t69;
t65 = -t42 * pkin(1) + pkin(7) * t74;
t16 = t22 * pkin(2);
t62 = t23 * qJ(3) - t16;
t18 = t24 * pkin(2);
t61 = t25 * qJ(3) - t18;
t53 = t22 * t44 + t40 * t74;
t7 = -t24 * t44 + t40 * t76;
t60 = g(1) * t53 + g(2) * t7;
t59 = pkin(4) * t40 - pkin(9) * t44;
t58 = g(1) * t22 - g(2) * t24;
t6 = g(1) * t23 - g(2) * t25;
t57 = g(1) * t46 + g(2) * t42;
t56 = t25 * pkin(2) + t24 * qJ(3) + t68;
t20 = -t67 * t40 - t44 * t75;
t52 = g(1) * t7 - g(2) * t53 - g(3) * t20;
t21 = -t40 * t75 + t67 * t44;
t8 = t24 * t40 + t44 * t76;
t51 = g(1) * t8 - g(2) * t54 + g(3) * t21;
t50 = -t23 * pkin(2) - t22 * qJ(3) + t65;
t4 = -g(1) * t24 - g(2) * t22 + g(3) * t75;
t49 = g(3) * t77 - t87;
t48 = pkin(3) * t76 + t25 * pkin(8) + t56;
t47 = pkin(3) * t74 - t23 * pkin(8) + t50;
t26 = t57 * t38;
t3 = t49 * t44;
t2 = t25 * t39 + t8 * t43;
t1 = t25 * t43 - t8 * t39;
t5 = [0, 0, 0, 0, 0, 0, g(1) * t42 - g(2) * t46, t57, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t58, -t26, -g(1) * t65 - g(2) * t68, 0, 0, 0, 0, 0, 0, -t26, -t6, t58, -g(1) * t50 - g(2) * t56, 0, 0, 0, 0, 0, 0, -g(1) * t54 - g(2) * t8, t60, t6, -g(1) * t47 - g(2) * t48, 0, 0, 0, 0, 0, 0, -g(1) * t85 - g(2) * t2, g(1) * t86 - g(2) * t1, -t60, -g(1) * (pkin(4) * t54 + pkin(9) * t53 + t47) - g(2) * (t8 * pkin(4) + t7 * pkin(9) + t48); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t49, -g(1) * t61 - g(2) * t62 - g(3) * t69, 0, 0, 0, 0, 0, 0, -t49 * t40, -t3, -t4, -g(1) * (t61 - t80) - g(2) * (t62 - t81) - g(3) * t66, 0, 0, 0, 0, 0, 0, -g(1) * (-t24 * t39 + t25 * t71) - g(2) * (-t22 * t39 + t23 * t71) - (t39 * t45 + t40 * t70) * t82, -g(1) * (-t24 * t43 - t25 * t73) - g(2) * (-t22 * t43 - t23 * t73) - (-t40 * t72 + t43 * t45) * t82, t3, -g(1) * (-t18 - t80) - g(2) * (-t16 - t81) - g(3) * (t59 * t77 + t66) + t87 * (qJ(3) + t59); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t51, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t43, -t52 * t39, -t51, -g(1) * (-t7 * pkin(4) + t8 * pkin(9)) - g(2) * (pkin(4) * t53 - pkin(9) * t54) - g(3) * (t20 * pkin(4) + t21 * pkin(9)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t86 - g(3) * (-t21 * t39 + t38 * t70), g(1) * t2 - g(2) * t85 - g(3) * (-t21 * t43 - t38 * t72), 0, 0;];
taug_reg = t5;
