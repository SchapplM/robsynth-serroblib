% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:02:56
% EndTime: 2019-05-06 12:02:58
% DurationCPUTime: 0.60s
% Computational Cost: add. (502->130), mult. (988->190), div. (0->0), fcn. (1180->12), ass. (0->71)
t51 = sin(pkin(11));
t100 = pkin(4) * t51;
t56 = sin(qJ(2));
t57 = sin(qJ(1));
t59 = cos(qJ(2));
t60 = cos(qJ(1));
t81 = cos(pkin(6));
t77 = t60 * t81;
t31 = t57 * t56 - t59 * t77;
t32 = t56 * t77 + t57 * t59;
t54 = -pkin(9) - qJ(4);
t105 = t32 * t100 + t31 * t54;
t78 = t57 * t81;
t33 = t60 * t56 + t59 * t78;
t34 = -t56 * t78 + t60 * t59;
t104 = t34 * t100 + t33 * t54;
t103 = -g(1) * t34 - g(2) * t32;
t55 = sin(qJ(6));
t58 = cos(qJ(6));
t50 = pkin(11) + qJ(5);
t47 = sin(t50);
t48 = cos(t50);
t52 = sin(pkin(6));
t87 = t52 * t60;
t68 = -t31 * t47 + t48 * t87;
t102 = t32 * t58 + t55 * t68;
t101 = -t32 * t55 + t58 * t68;
t97 = g(3) * t52;
t96 = t31 * t51;
t93 = t33 * t51;
t92 = t47 * t55;
t91 = t47 * t58;
t90 = t52 * t56;
t89 = t52 * t57;
t88 = t52 * t59;
t86 = t54 * t59;
t85 = t55 * t56;
t84 = t56 * t58;
t83 = pkin(2) * t88 + qJ(3) * t90;
t82 = t60 * pkin(1) + pkin(8) * t89;
t80 = t90 * t100 + t83;
t79 = -t57 * pkin(1) + pkin(8) * t87;
t27 = t31 * pkin(2);
t76 = t32 * qJ(3) - t27;
t29 = t33 * pkin(2);
t75 = t34 * qJ(3) - t29;
t67 = t31 * t48 + t47 * t87;
t8 = -t33 * t48 + t47 * t89;
t74 = g(1) * t67 + g(2) * t8;
t73 = pkin(5) * t47 - pkin(10) * t48;
t72 = g(1) * t31 - g(2) * t33;
t7 = g(1) * t32 - g(2) * t34;
t71 = g(1) * t60 + g(2) * t57;
t70 = t34 * pkin(2) + t33 * qJ(3) + t82;
t17 = -t81 * t47 - t48 * t88;
t66 = g(1) * t8 - g(2) * t67 - g(3) * t17;
t18 = -t47 * t88 + t81 * t48;
t9 = t33 * t47 + t48 * t89;
t65 = g(1) * t9 - g(2) * t68 + g(3) * t18;
t64 = -t32 * pkin(2) - t31 * qJ(3) + t79;
t4 = -g(1) * t33 - g(2) * t31 + g(3) * t88;
t63 = g(3) * t90 - t103;
t53 = cos(pkin(11));
t46 = t53 * pkin(4) + pkin(3);
t62 = pkin(4) * t93 - t34 * t54 + t46 * t89 + t70;
t61 = -pkin(4) * t96 + t32 * t54 + t46 * t87 + t64;
t35 = t71 * t52;
t3 = t34 * t55 + t9 * t58;
t2 = t34 * t58 - t9 * t55;
t1 = t63 * t48;
t5 = [0, 0, 0, 0, 0, 0, g(1) * t57 - g(2) * t60, t71, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t72, -t35, -g(1) * t79 - g(2) * t82, 0, 0, 0, 0, 0, 0, -t35, -t7, t72, -g(1) * t64 - g(2) * t70, 0, 0, 0, 0, 0, 0, -g(1) * (t53 * t87 - t96) - g(2) * (t53 * t89 + t93) -g(1) * (-t31 * t53 - t51 * t87) - g(2) * (t33 * t53 - t51 * t89) t7, -g(1) * (pkin(3) * t87 - t32 * qJ(4) + t64) - g(2) * (pkin(3) * t89 + t34 * qJ(4) + t70) 0, 0, 0, 0, 0, 0, -g(1) * t68 - g(2) * t9, t74, t7, -g(1) * t61 - g(2) * t62, 0, 0, 0, 0, 0, 0, -g(1) * t101 - g(2) * t3, g(1) * t102 - g(2) * t2, -t74, -g(1) * (pkin(5) * t68 + pkin(10) * t67 + t61) - g(2) * (t9 * pkin(5) + t8 * pkin(10) + t62); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t63, -g(1) * t75 - g(2) * t76 - g(3) * t83, 0, 0, 0, 0, 0, 0, -t63 * t51, -t63 * t53, -t4, -g(1) * (-t33 * qJ(4) + t75) - g(2) * (-t31 * qJ(4) + t76) - g(3) * (qJ(4) * t88 + t83) 0, 0, 0, 0, 0, 0, -t63 * t47, -t1, -t4, -g(1) * (t75 + t104) - g(2) * (t76 + t105) - g(3) * (-t52 * t86 + t80) 0, 0, 0, 0, 0, 0, -g(1) * (-t33 * t55 + t34 * t91) - g(2) * (-t31 * t55 + t32 * t91) - (t47 * t84 + t55 * t59) * t97, -g(1) * (-t33 * t58 - t34 * t92) - g(2) * (-t31 * t58 - t32 * t92) - (-t47 * t85 + t58 * t59) * t97, t1, -g(1) * (-t29 + t104) - g(2) * (-t27 + t105) - g(3) * t80 - (t73 * t56 - t86) * t97 + t103 * (qJ(3) + t73); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t65, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t58, -t66 * t55, -t65, -g(1) * (-t8 * pkin(5) + t9 * pkin(10)) - g(2) * (pkin(5) * t67 - pkin(10) * t68) - g(3) * (t17 * pkin(5) + t18 * pkin(10)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t102 - g(3) * (-t18 * t55 + t52 * t84) g(1) * t3 - g(2) * t101 - g(3) * (-t18 * t58 - t52 * t85) 0, 0;];
taug_reg  = t5;
