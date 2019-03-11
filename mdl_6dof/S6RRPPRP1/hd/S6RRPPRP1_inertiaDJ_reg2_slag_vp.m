% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRP1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:27:47
% EndTime: 2019-03-09 08:27:55
% DurationCPUTime: 2.34s
% Computational Cost: add. (4472->230), mult. (10092->409), div. (0->0), fcn. (9863->8), ass. (0->130)
t98 = sin(pkin(10));
t96 = t98 ^ 2;
t100 = cos(pkin(10));
t97 = t100 ^ 2;
t154 = t97 + t96;
t179 = 0.2e1 * t154 * qJD(4);
t102 = sin(qJ(5));
t169 = cos(qJ(5));
t130 = qJD(5) * t169;
t144 = qJD(5) * t102;
t136 = t98 * t144;
t137 = t169 * t98;
t101 = cos(pkin(9));
t104 = cos(qJ(2));
t145 = qJD(2) * t104;
t103 = sin(qJ(2));
t146 = qJD(2) * t103;
t99 = sin(pkin(9));
t77 = t101 * t145 - t99 * t146;
t83 = t101 * t103 + t99 * t104;
t30 = t77 * t137 - t83 * t136 + (t102 * t77 + t83 * t130) * t100;
t148 = t102 * t100;
t84 = t137 + t148;
t49 = t84 * t83;
t78 = -t100 * t130 + t136;
t157 = -t84 * t30 + t78 * t49;
t174 = t169 * t100 - t102 * t98;
t79 = t84 * qJD(5);
t29 = -t174 * t77 + t83 * t79;
t50 = t174 * t83;
t175 = t174 * t29 + t50 * t79;
t178 = t175 - t157;
t177 = t175 + t157;
t76 = t83 * qJD(2);
t81 = -t101 * t104 + t99 * t103;
t120 = t83 * t76 + t77 * t81;
t176 = -0.2e1 * t120;
t138 = pkin(2) * t146;
t111 = t76 * pkin(3) - t83 * qJD(4) + t138;
t108 = -t77 * qJ(4) + t111;
t159 = -qJ(3) - pkin(7);
t129 = qJD(2) * t159;
t112 = -t103 * qJD(3) + t104 * t129;
t75 = t104 * qJD(3) + t103 * t129;
t47 = t101 * t75 + t99 * t112;
t162 = t98 * t47;
t19 = t100 * t108 - t162;
t153 = t100 * t47;
t20 = t98 * t108 + t153;
t173 = t20 * t100 - t19 * t98;
t155 = -t174 * t78 - t84 * t79;
t171 = 2 * qJD(6);
t170 = t76 * pkin(5);
t46 = -t101 * t112 + t99 * t75;
t87 = t159 * t103;
t88 = t159 * t104;
t64 = -t101 * t87 - t99 * t88;
t167 = t64 * t46;
t165 = t174 * t79;
t164 = t83 * t98;
t161 = t98 * t76;
t160 = t98 * t77;
t95 = -t104 * pkin(2) - pkin(1);
t55 = t81 * pkin(3) - t83 * qJ(4) + t95;
t65 = -t101 * t88 + t99 * t87;
t31 = t100 * t55 - t98 * t65;
t22 = -t100 * t83 * pkin(8) + t81 * pkin(4) + t31;
t32 = t100 * t65 + t98 * t55;
t28 = -pkin(8) * t164 + t32;
t9 = t102 * t22 + t169 * t28;
t152 = t100 * t77;
t149 = t76 * qJ(6);
t147 = t81 * qJD(6);
t143 = 0.2e1 * t49 * t30;
t56 = 0.2e1 * t81 * t76;
t142 = -0.2e1 * t165;
t141 = 0.2e1 * t83 * t77;
t140 = -0.2e1 * pkin(1) * qJD(2);
t139 = t98 * t152;
t94 = -t101 * pkin(2) - pkin(3);
t134 = t103 * t145;
t133 = t99 * pkin(2) + qJ(4);
t128 = pkin(8) + t133;
t117 = t128 * t98;
t113 = t169 * t117;
t131 = qJD(4) * t169;
t80 = t128 * t100;
t40 = qJD(5) * t113 - t100 * t131 + (qJD(4) * t98 + qJD(5) * t80) * t102;
t116 = t102 * t117;
t41 = qJD(4) * t148 - qJD(5) * t116 + t80 * t130 + t98 * t131;
t53 = t102 * t80 + t113;
t54 = t169 * t80 - t116;
t132 = -t54 * t40 + t53 * t41;
t35 = pkin(4) * t160 + t46;
t45 = pkin(4) * t164 + t64;
t127 = t29 * t49 - t50 * t30;
t125 = t81 * t30 + t76 * t49;
t124 = -t174 * t30 + t49 * t79;
t123 = t40 * t81 - t54 * t76;
t122 = -t41 * t81 - t53 * t76;
t121 = t46 * t83 + t64 * t77;
t36 = t84 * t76 - t78 * t81;
t119 = -t174 * t76 + t81 * t79;
t61 = t84 * t78;
t118 = -0.2e1 * t61 - 0.2e1 * t165;
t86 = -t100 * pkin(4) + t94;
t8 = -t102 * t28 + t169 * t22;
t107 = (-pkin(8) - qJ(4)) * t77 + t111;
t105 = t76 * pkin(4) + t107 * t100 - t162;
t106 = t107 * t98 + t153;
t3 = -t102 * t105 - t169 * t106 - t22 * t130 + t28 * t144;
t115 = -t53 * t29 - t54 * t30 + t40 * t49 + t41 * t50;
t114 = -t174 * t41 - t40 * t84 + t53 * t79 - t54 * t78;
t44 = t79 * pkin(5) + t78 * qJ(6) - t84 * qJD(6);
t110 = -0.2e1 * t174 * t40 + 0.2e1 * t41 * t84 - 0.2e1 * t53 * t78 - 0.2e1 * t54 * t79;
t109 = -qJD(4) * t81 - t133 * t76 + t94 * t77;
t4 = -t102 * t106 + t169 * t105 - t28 * t130 - t22 * t144;
t71 = t100 * t76;
t57 = -0.2e1 * t61;
t48 = -pkin(5) * t174 - t84 * qJ(6) + t86;
t16 = t49 * pkin(5) - t50 * qJ(6) + t45;
t15 = -0.2e1 * t50 * t29;
t11 = -t29 * t84 - t50 * t78;
t10 = -0.2e1 * t29 * t81 + 0.2e1 * t50 * t76;
t7 = -t81 * pkin(5) - t8;
t6 = t81 * qJ(6) + t9;
t5 = t30 * pkin(5) + t29 * qJ(6) - t50 * qJD(6) + t35;
t2 = -t4 - t170;
t1 = t147 - t3 + t149;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t134, 0.2e1 * (-t103 ^ 2 + t104 ^ 2) * qJD(2), 0, -0.2e1 * t134, 0, 0, t103 * t140, t104 * t140, 0, 0, t141, t176, 0, t56, 0, 0, 0.2e1 * t81 * t138 + 0.2e1 * t95 * t76, 0.2e1 * t83 * t138 + 0.2e1 * t95 * t77, -0.2e1 * t47 * t81 - 0.2e1 * t65 * t76 + 0.2e1 * t121, 0.2e1 * t95 * t138 + 0.2e1 * t65 * t47 + 0.2e1 * t167, t97 * t141, -0.4e1 * t83 * t139, 0.2e1 * t120 * t100, t96 * t141, t98 * t176, t56, 0.2e1 * t121 * t98 + 0.2e1 * t19 * t81 + 0.2e1 * t31 * t76, 0.2e1 * t121 * t100 - 0.2e1 * t20 * t81 - 0.2e1 * t32 * t76, 0.2e1 * (-t20 * t83 - t32 * t77) * t98 + 0.2e1 * (-t19 * t83 - t31 * t77) * t100, 0.2e1 * t31 * t19 + 0.2e1 * t32 * t20 + 0.2e1 * t167, t15, 0.2e1 * t127, t10, t143, -0.2e1 * t125, t56, 0.2e1 * t45 * t30 + 0.2e1 * t35 * t49 + 0.2e1 * t4 * t81 + 0.2e1 * t8 * t76, -0.2e1 * t45 * t29 + 0.2e1 * t3 * t81 + 0.2e1 * t35 * t50 - 0.2e1 * t9 * t76, 0.2e1 * t8 * t29 + 0.2e1 * t3 * t49 - 0.2e1 * t9 * t30 - 0.2e1 * t4 * t50, -0.2e1 * t9 * t3 + 0.2e1 * t45 * t35 + 0.2e1 * t8 * t4, t15, t10, -0.2e1 * t127, t56, 0.2e1 * t125, t143, 0.2e1 * t16 * t30 - 0.2e1 * t2 * t81 + 0.2e1 * t5 * t49 - 0.2e1 * t7 * t76, -0.2e1 * t1 * t49 + 0.2e1 * t2 * t50 - 0.2e1 * t7 * t29 - 0.2e1 * t6 * t30, 0.2e1 * t1 * t81 + 0.2e1 * t16 * t29 - 0.2e1 * t5 * t50 + 0.2e1 * t6 * t76, 0.2e1 * t6 * t1 + 0.2e1 * t16 * t5 + 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, 0, -t146, 0, -pkin(7) * t145, pkin(7) * t146, 0, 0, 0, 0, t77, 0, -t76, 0, -t46, -t47 (-t101 * t77 - t76 * t99) * pkin(2) (-t101 * t46 + t47 * t99) * pkin(2), t139 (-t96 + t97) * t77, t161, -t139, t71, 0, -t46 * t100 + t109 * t98, t109 * t100 + t46 * t98, t173, t46 * t94 + (t32 * t100 - t31 * t98) * qJD(4) + t173 * t133, t11, -t178, t36, t124, -t119, 0, -t174 * t35 + t86 * t30 + t45 * t79 + t122, -t86 * t29 + t35 * t84 - t45 * t78 + t123, -t174 * t3 - t4 * t84 + t8 * t78 - t9 * t79 + t115, -t3 * t54 + t35 * t86 - t4 * t53 - t9 * t40 - t8 * t41, t11, t36, t178, 0, t119, t124, t16 * t79 - t174 * t5 + t48 * t30 + t44 * t49 + t122, t1 * t174 + t2 * t84 - t6 * t79 - t7 * t78 + t115, t16 * t78 + t48 * t29 - t44 * t50 - t5 * t84 - t123, t1 * t54 + t16 * t44 + t2 * t53 - t6 * t40 + t7 * t41 + t5 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, t133 * t179, t57, 0.2e1 * t155, 0, t142, 0, 0, 0.2e1 * t86 * t79, -0.2e1 * t86 * t78, t110, 0.2e1 * t132, t57, 0, -0.2e1 * t155, 0, 0, t142, -0.2e1 * t174 * t44 + 0.2e1 * t48 * t79, t110, -0.2e1 * t44 * t84 + 0.2e1 * t48 * t78, 0.2e1 * t48 * t44 + 0.2e1 * t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t77, 0, t138, 0, 0, 0, 0, 0, 0, t71, -t161, -t154 * t77, t19 * t100 + t20 * t98, 0, 0, 0, 0, 0, 0, -t119, -t36, t177, t174 * t4 - t3 * t84 - t9 * t78 - t8 * t79, 0, 0, 0, 0, 0, 0, -t119, t177, t36, t1 * t84 - t174 * t2 - t6 * t78 + t7 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, t152, 0, t46, 0, 0, 0, 0, 0, 0, t30, -t29, 0, t35, 0, 0, 0, 0, 0, 0, t30, 0, t29, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t78, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, t78, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, -t30, t76, t4, t3, 0, 0, 0, -t29, 0, t76, t30, 0, t4 + 0.2e1 * t170, pkin(5) * t29 - t30 * qJ(6) - t49 * qJD(6), 0.2e1 * t147 - t3 + 0.2e1 * t149, -t2 * pkin(5) + t1 * qJ(6) + t6 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, 0, -t79, 0, -t41, t40, 0, 0, 0, -t78, 0, 0, t79, 0, -t41, pkin(5) * t78 - t79 * qJ(6) + qJD(6) * t174, -t40, -t41 * pkin(5) - t40 * qJ(6) + t54 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t78, 0, 0, 0, 0, 0, 0, 0, 0, -t79, 0, -t78, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, qJ(6) * t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t29, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, 0, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t12;
