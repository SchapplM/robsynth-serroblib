% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x34]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:18:00
% EndTime: 2019-03-09 07:18:04
% DurationCPUTime: 1.33s
% Computational Cost: add. (2744->171), mult. (5617->281), div. (0->0), fcn. (5531->8), ass. (0->117)
t76 = sin(qJ(4));
t77 = sin(qJ(3));
t80 = cos(qJ(4));
t81 = cos(qJ(3));
t52 = t76 * t81 + t80 * t77;
t75 = sin(qJ(5));
t79 = cos(qJ(5));
t95 = t76 * t77 - t80 * t81;
t38 = t52 * t75 + t79 * t95;
t37 = t38 ^ 2;
t122 = qJD(5) * t79;
t129 = t76 * t79;
t147 = ((t75 * t80 + t129) * qJD(4) + t76 * t122) * pkin(3);
t140 = -qJD(3) - qJD(4);
t42 = t140 * t52;
t99 = t79 * t52 - t75 * t95;
t139 = -t99 * qJD(5) + t79 * t42;
t90 = t95 * qJD(3);
t83 = t95 * qJD(4) + t90;
t14 = t75 * t83 + t139;
t133 = t14 * t38;
t78 = cos(qJ(6));
t70 = qJD(6) * t78;
t146 = t38 * t70;
t74 = sin(qJ(6));
t121 = qJD(6) * t74;
t145 = t38 * t121;
t82 = -pkin(1) - pkin(7);
t135 = pkin(8) - t82;
t54 = t135 * t77;
t55 = t135 * t81;
t97 = -t76 * t54 + t80 * t55;
t28 = t95 * pkin(9) - t97;
t144 = t70 * t99;
t137 = -t38 * qJD(5) + t42 * t75;
t16 = -t79 * t83 + t137;
t143 = t99 * t16;
t142 = t121 * t99;
t120 = t77 * qJD(3);
t49 = t135 * t120;
t50 = qJD(3) * t55;
t100 = t76 * t49 - t80 * t50;
t141 = t28 * qJD(4) + t100;
t73 = t78 ^ 2;
t125 = t74 ^ 2 - t73;
t106 = t125 * qJD(6);
t136 = 2 * qJD(2);
t98 = t54 * t80 + t55 * t76;
t29 = -pkin(9) * t52 - t98;
t22 = -t28 * t79 + t29 * t75;
t119 = t81 * qJD(3);
t23 = t75 * t28 + t79 * t29;
t26 = t98 * qJD(4) + t49 * t80 + t50 * t76;
t85 = pkin(9) * t42 - t26;
t5 = t23 * qJD(5) + t79 * t85 + ((-t80 * t119 + t76 * t120) * pkin(9) + t141) * t75;
t134 = t22 * t70 + t5 * t74;
t132 = t38 * t74;
t131 = t38 * t78;
t128 = t78 * t14;
t123 = qJD(5) * t75;
t68 = pkin(3) * t80 + pkin(4);
t33 = t68 * t123 + t147;
t117 = pkin(3) * t75 * t76;
t46 = -t68 * t79 - pkin(5) + t117;
t127 = t33 * t74 + t46 * t70;
t112 = pkin(4) * t123;
t67 = -pkin(4) * t79 - pkin(5);
t126 = t74 * t112 + t67 * t70;
t124 = pkin(3) * qJD(4);
t63 = t77 * pkin(3) + qJ(2);
t62 = pkin(3) * t119 + qJD(2);
t118 = qJ(2) * qJD(3);
t116 = pkin(5) * t121;
t115 = pkin(5) * t70;
t114 = t76 * t124;
t113 = t80 * t124;
t111 = pkin(4) * t122;
t109 = t74 * t70;
t108 = 0.4e1 * t74 * t131;
t43 = t46 * t121;
t107 = -t33 * t78 + t43;
t45 = pkin(4) * t52 + t63;
t24 = pkin(5) * t99 + pkin(10) * t38 + t45;
t105 = t23 * t78 + t24 * t74;
t104 = t23 * t74 - t24 * t78;
t103 = -t99 ^ 2 - t37;
t47 = pkin(3) * t129 + t68 * t75 + pkin(10);
t102 = t38 * t46 + t47 * t99;
t66 = pkin(4) * t75 + pkin(10);
t101 = t38 * t67 + t66 * t99;
t56 = t67 * t121;
t94 = -t78 * t112 + t56;
t93 = t74 * t14 - t146;
t92 = -t128 - t145;
t91 = -t16 * t78 + t142;
t32 = -t68 * t122 - t79 * t113 + (qJD(4) + qJD(5)) * t117;
t84 = t140 * t95;
t13 = -t79 * t84 - t137;
t15 = t75 * t84 - t139;
t89 = t13 * t99 - t15 * t38 + t133 - t143;
t88 = t14 * t46 - t16 * t47 + t32 * t99 - t33 * t38;
t87 = t14 * t67 - t16 * t66 + (-t38 * t75 - t79 * t99) * qJD(5) * pkin(4);
t30 = -t83 * pkin(4) + t62;
t61 = 0.2e1 * t109;
t51 = -0.2e1 * t106;
t25 = t97 * qJD(4) - t100;
t20 = t22 * t121;
t12 = -t15 * t78 + t145;
t11 = t15 * t74 + t146;
t10 = t16 * t74 + t144;
t8 = t106 * t38 + t74 * t128;
t7 = t16 * pkin(5) - t14 * pkin(10) + t30;
t6 = qJD(6) * t108 - t125 * t14;
t4 = t29 * t123 + t75 * t85 - t79 * (pkin(9) * t90 + t141) - t28 * t122;
t2 = -t105 * qJD(6) + t4 * t74 + t7 * t78;
t1 = t104 * qJD(6) + t4 * t78 - t7 * t74;
t3 = [0, 0, 0, 0, t136, qJ(2) * t136, -0.2e1 * t77 * t119, 0.2e1 * (t77 ^ 2 - t81 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t77 + 0.2e1 * t81 * t118, 0.2e1 * qJD(2) * t81 - 0.2e1 * t77 * t118, -0.2e1 * t95 * t42, -0.2e1 * t42 * t52 - 0.2e1 * t83 * t95, 0, 0, 0, 0.2e1 * t62 * t52 + 0.2e1 * t63 * t84, 0.2e1 * t42 * t63 - 0.2e1 * t62 * t95, -0.2e1 * t133, -0.2e1 * t14 * t99 + 0.2e1 * t16 * t38, 0, 0, 0, 0.2e1 * t16 * t45 + 0.2e1 * t30 * t99, 0.2e1 * t14 * t45 - 0.2e1 * t30 * t38, -0.2e1 * t37 * t109 - 0.2e1 * t73 * t133, 0.2e1 * t37 * t106 + t14 * t108, -0.2e1 * t16 * t131 - 0.2e1 * t92 * t99, 0.2e1 * t16 * t132 - 0.2e1 * t93 * t99, 0.2e1 * t143, -0.2e1 * t104 * t16 - 0.2e1 * t5 * t132 + 0.2e1 * t2 * t99 + 0.2e1 * t93 * t22, 0.2e1 * t1 * t99 - 0.2e1 * t105 * t16 - 0.2e1 * t5 * t131 - 0.2e1 * t92 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103 * t70 + t89 * t74, -t103 * t121 + t89 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t119, 0, -t82 * t120, -t82 * t119, 0, 0, t42, t83, 0, t26, t25, 0, 0, t14, -t16, 0, -t5, t4, t8, t6, t10, -t91, 0, t20 + (-t102 * qJD(6) - t5) * t78 + t88 * t74, t102 * t121 + t88 * t78 + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t119, 0, 0, 0, 0, 0, t42, t83, 0, 0, 0, 0, 0, -t15, t13, 0, 0, 0, 0, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t114, -0.2e1 * t113, 0, 0, 0, 0, 0, -0.2e1 * t33, 0.2e1 * t32, t61, t51, 0, 0, 0, 0.2e1 * t107, 0.2e1 * t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t83, 0, t26, t25, 0, 0, t14, -t16, 0, -t5, t4, t8, t6, t10, -t91, 0, t20 + (-t101 * qJD(6) - t5) * t78 + t87 * t74, t101 * t121 + t87 * t78 + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t83, 0, 0, 0, 0, 0, -t15, t13, 0, 0, 0, 0, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, -t113, 0, 0, 0, 0, 0 (-pkin(4) - t68) * t123 - t147, t32 - t111, t61, t51, 0, 0, 0, t43 + t56 + (-t33 - t112) * t78, t126 + t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t112, -0.2e1 * t111, t61, t51, 0, 0, 0, 0.2e1 * t94, 0.2e1 * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t16, 0, -t5, t4, t8, t6, t10, -t91, 0, t20 + (-pkin(5) * t14 - pkin(10) * t16) * t74 + (-t5 + (pkin(5) * t38 - pkin(10) * t99) * qJD(6)) * t78, t92 * pkin(5) + t91 * pkin(10) + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t13, 0, 0, 0, 0, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32, t61, t51, 0, 0, 0, t107 - t116, -t115 + t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t111, t61, t51, 0, 0, 0, t94 - t116, -t115 + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t51, 0, 0, 0, -0.2e1 * t116, -0.2e1 * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t93, t16, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t74 - t144, t13 * t78 + t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t121, 0, t32 * t74 - t47 * t70, t47 * t121 + t32 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t121, 0, -t74 * t111 - t66 * t70, -t78 * t111 + t66 * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t121, 0, -pkin(10) * t70, pkin(10) * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
