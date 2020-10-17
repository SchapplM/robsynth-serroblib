% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:16:28
% EndTime: 2019-05-06 18:16:30
% DurationCPUTime: 0.50s
% Computational Cost: add. (357->107), mult. (912->128), div. (0->0), fcn. (1040->8), ass. (0->67)
t52 = sin(qJ(4));
t56 = cos(qJ(2));
t53 = sin(qJ(2));
t85 = cos(qJ(4));
t77 = t53 * t85;
t33 = -t56 * t52 + t77;
t54 = sin(qJ(1));
t24 = t33 * t54;
t51 = sin(qJ(5));
t55 = cos(qJ(5));
t66 = pkin(5) * t55 + qJ(6) * t51;
t104 = t24 * t66;
t57 = cos(qJ(1));
t83 = t56 * t57;
t26 = t52 * t83 - t57 * t77;
t103 = t26 * t66;
t32 = t53 * t52 + t56 * t85;
t102 = t32 * t66;
t34 = g(1) * t57 + g(2) * t54;
t101 = t34 * t53;
t62 = g(1) * t26 - g(2) * t24 + g(3) * t32;
t100 = t62 * t51;
t99 = t62 * t55;
t44 = t53 * qJ(3);
t82 = t56 * pkin(2) + t44;
t73 = -t32 * pkin(4) + t33 * pkin(9);
t27 = t32 * t57;
t74 = -t26 * pkin(4) + t27 * pkin(9);
t25 = t32 * t54;
t75 = t24 * pkin(4) + t25 * pkin(9);
t97 = pkin(2) * t53;
t95 = g(1) * t54;
t92 = g(3) * t33;
t46 = t56 * pkin(3);
t81 = t57 * pkin(1) + t54 * pkin(7);
t80 = qJ(3) * t56;
t78 = t46 + t82;
t48 = t57 * pkin(7);
t76 = -t57 * pkin(8) + t48;
t72 = pkin(2) * t83 + t57 * t44 + t81;
t37 = t54 * t80;
t71 = t37 - t75;
t39 = t57 * t80;
t70 = t39 - t74;
t69 = pkin(3) * t83 + t72;
t10 = t25 * t51 - t57 * t55;
t14 = t27 * t51 + t54 * t55;
t68 = -g(1) * t10 + g(2) * t14;
t9 = -g(1) * t24 - g(2) * t26;
t67 = -g(2) * t57 + t95;
t11 = t25 * t55 + t57 * t51;
t65 = -pkin(1) - t82;
t64 = t78 - t73;
t63 = -t25 * pkin(4) + t24 * pkin(9) + t76;
t8 = g(1) * t27 + g(2) * t25 + t92;
t1 = g(1) * t14 + g(2) * t10 + t51 * t92;
t15 = t27 * t55 - t54 * t51;
t61 = g(1) * t15 + g(2) * t11 + t55 * t92;
t60 = t27 * pkin(4) + t26 * pkin(9) + t69;
t59 = (pkin(2) + pkin(3)) * t101;
t58 = (-g(1) * (t65 - t46) + g(2) * pkin(8)) * t54;
t29 = t67 * t56;
t28 = t67 * t53;
t23 = g(3) * t53 + t34 * t56;
t22 = -g(3) * t56 + t101;
t6 = g(1) * t11 - g(2) * t15;
t2 = [0, 0, 0, 0, 0, 0, t67, t34, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, -t34, -g(1) * (-t54 * pkin(1) + t48) - g(2) * t81, 0, 0, 0, 0, 0, 0, t29, -t34, t28, -g(1) * t48 - g(2) * t72 - t65 * t95, 0, 0, 0, 0, 0, 0, g(1) * t25 - g(2) * t27, -t9, t34, -g(1) * t76 - g(2) * t69 + t58, 0, 0, 0, 0, 0, 0, t6, t68, t9, -g(1) * t63 - g(2) * t60 + t58, 0, 0, 0, 0, 0, 0, t6, t9, -t68, -g(1) * (-pkin(5) * t11 - qJ(6) * t10 + t63) - g(2) * (t15 * pkin(5) + t14 * qJ(6) + t60) + t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, -t23, -g(1) * (-t57 * t97 + t39) - g(2) * (-t54 * t97 + t37) - g(3) * t82, 0, 0, 0, 0, 0, 0, -t62, -t8, 0, -g(1) * t39 - g(2) * t37 - g(3) * t78 + t59, 0, 0, 0, 0, 0, 0, -t99, t100, t8, -g(1) * t70 - g(2) * t71 - g(3) * t64 + t59, 0, 0, 0, 0, 0, 0, -t99, t8, -t100, -g(1) * (t70 + t103) - g(2) * (t71 - t104) - g(3) * (t64 + t102) + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t8, 0, 0, 0, 0, 0, 0, 0, 0, t99, -t100, -t8, -g(1) * t74 - g(2) * t75 - g(3) * t73, 0, 0, 0, 0, 0, 0, t99, -t8, t100, -g(1) * (t74 - t103) - g(2) * (t75 + t104) - g(3) * (t73 - t102); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t61, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t61, -g(1) * (-t14 * pkin(5) + t15 * qJ(6)) - g(2) * (-t10 * pkin(5) + t11 * qJ(6)) - (-pkin(5) * t51 + qJ(6) * t55) * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
