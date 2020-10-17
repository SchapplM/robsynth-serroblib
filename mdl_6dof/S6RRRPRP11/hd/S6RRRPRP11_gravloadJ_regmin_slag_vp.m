% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:18:40
% EndTime: 2019-05-07 09:18:42
% DurationCPUTime: 0.48s
% Computational Cost: add. (384->113), mult. (979->165), div. (0->0), fcn. (1195->10), ass. (0->61)
t51 = sin(qJ(3));
t55 = cos(qJ(3));
t94 = pkin(3) * t55 + qJ(4) * t51 + pkin(2);
t52 = sin(qJ(2));
t53 = sin(qJ(1));
t56 = cos(qJ(2));
t74 = cos(pkin(6));
t86 = cos(qJ(1));
t64 = t74 * t86;
t31 = t52 * t64 + t53 * t56;
t48 = sin(pkin(6));
t71 = t48 * t86;
t13 = t31 * t51 + t55 * t71;
t30 = t53 * t52 - t56 * t64;
t50 = sin(qJ(5));
t54 = cos(qJ(5));
t93 = t13 * t50 + t30 * t54;
t92 = t13 * t54 - t30 * t50;
t83 = t48 * t52;
t28 = t51 * t83 - t74 * t55;
t69 = t53 * t74;
t33 = -t52 * t69 + t86 * t56;
t81 = t48 * t55;
t17 = t33 * t51 - t53 * t81;
t32 = t86 * t52 + t56 * t69;
t5 = t17 * t54 - t32 * t50;
t78 = t50 * t56;
t91 = -g(2) * t92 - g(3) * (t28 * t54 + t48 * t78) - g(1) * t5;
t88 = g(3) * t48;
t46 = t54 * pkin(5) + pkin(4);
t87 = pkin(9) + t46;
t82 = t48 * t53;
t80 = t48 * t56;
t79 = t50 * t51;
t77 = t51 * t54;
t76 = t54 * t56;
t73 = t94 * t30;
t72 = t94 * t32;
t70 = pkin(5) * t50 + qJ(4);
t14 = t31 * t55 - t51 * t71;
t18 = t33 * t55 + t51 * t82;
t68 = t86 * pkin(1) + t33 * pkin(2) + t18 * pkin(3) + pkin(8) * t82;
t67 = -g(1) * t13 + g(2) * t17;
t66 = -g(1) * t14 + g(2) * t18;
t65 = g(1) * t30 - g(2) * t32;
t63 = g(3) * (pkin(9) * t83 + t94 * t80);
t62 = -t53 * pkin(1) - t31 * pkin(2) - pkin(3) * t14 + pkin(8) * t71;
t49 = -qJ(6) - pkin(10);
t61 = pkin(5) * t79 - t49 * t55;
t2 = g(1) * t17 + g(2) * t13 + g(3) * t28;
t29 = t74 * t51 + t52 * t81;
t60 = g(1) * t18 + g(2) * t14 + g(3) * t29;
t59 = -g(1) * t32 - g(2) * t30 + g(3) * t80;
t58 = g(1) * t33 + g(2) * t31 + g(3) * t83;
t23 = t28 * pkin(3);
t11 = t17 * pkin(3);
t9 = t13 * pkin(3);
t8 = t59 * t55;
t7 = t59 * t51;
t6 = t17 * t50 + t32 * t54;
t1 = [0, g(1) * t53 - g(2) * t86, g(1) * t86 + g(2) * t53, 0, 0, 0, 0, 0, g(1) * t31 - g(2) * t33, -t65, 0, 0, 0, 0, 0, -t66, t67, t65, t66, -t67, -g(1) * (-t30 * pkin(9) - qJ(4) * t13 + t62) - g(2) * (t32 * pkin(9) + t17 * qJ(4) + t68) 0, 0, 0, 0, 0, g(1) * t93 - g(2) * t6, g(1) * t92 - g(2) * t5, -t66, -g(1) * (-t13 * t70 + t14 * t49 - t87 * t30 + t62) - g(2) * (t70 * t17 - t18 * t49 + t87 * t32 + t68); 0, 0, 0, 0, 0, 0, 0, 0, -t59, t58, 0, 0, 0, 0, 0, -t8, t7, -t58, t8, -t7, -g(1) * (t33 * pkin(9) - t72) - g(2) * (t31 * pkin(9) - t73) - t63, 0, 0, 0, 0, 0, -g(1) * (-t32 * t79 + t33 * t54) - g(2) * (-t30 * t79 + t31 * t54) - (t51 * t78 + t52 * t54) * t88, -g(1) * (-t32 * t77 - t33 * t50) - g(2) * (-t30 * t77 - t31 * t50) - (-t50 * t52 + t51 * t76) * t88, -t8, -g(1) * (-t61 * t32 + t87 * t33 - t72) - g(2) * (-t61 * t30 + t87 * t31 - t73) - t63 - (t46 * t52 + t61 * t56) * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t60, 0, -t2, -t60, -g(1) * (t18 * qJ(4) - t11) - g(2) * (t14 * qJ(4) - t9) - g(3) * (t29 * qJ(4) - t23) 0, 0, 0, 0, 0, -t60 * t50, -t60 * t54, t2, -g(1) * (t17 * t49 + t70 * t18 - t11) - g(2) * (t13 * t49 + t70 * t14 - t9) - g(3) * (t28 * t49 + t70 * t29 - t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, g(1) * t6 + g(2) * t93 - g(3) * (-t28 * t50 + t48 * t76) 0, t91 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60;];
taug_reg  = t1;
