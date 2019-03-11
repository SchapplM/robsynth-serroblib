% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:51:33
% EndTime: 2019-03-09 01:51:36
% DurationCPUTime: 1.09s
% Computational Cost: add. (749->207), mult. (1482->277), div. (0->0), fcn. (744->4), ass. (0->118)
t67 = cos(qJ(4));
t112 = t67 * qJD(1);
t47 = qJD(6) + t112;
t139 = qJD(6) - t47;
t64 = sin(qJ(6));
t114 = t64 * qJD(4);
t65 = sin(qJ(4));
t119 = qJD(1) * t65;
t66 = cos(qJ(6));
t138 = t119 * t66 - t114;
t54 = qJ(2) * qJD(1) + qJD(3);
t42 = -pkin(7) * qJD(1) + t54;
t22 = (pkin(5) * qJD(1) - t42) * t67;
t108 = qJD(5) + t22;
t79 = (t47 + t112) * t65;
t63 = pkin(1) + qJ(3);
t137 = qJD(1) * t63;
t121 = t67 * t42;
t106 = qJD(1) * qJD(2);
t48 = t65 * t106;
t14 = -t48 + (-qJD(5) - t121) * qJD(4);
t118 = qJD(4) * t65;
t27 = t42 * t118;
t17 = -t106 * t67 + t27;
t90 = qJD(4) * pkin(4) - qJD(5);
t24 = -t90 - t121;
t107 = qJD(4) * qJ(5);
t32 = t65 * t42;
t26 = -t32 - t107;
t71 = -t14 * t65 - t17 * t67 + (t24 * t65 - t26 * t67) * qJD(4);
t113 = t66 * qJD(4);
t30 = t119 * t64 + t113;
t105 = qJD(1) * qJD(4);
t94 = t67 * t105;
t12 = qJD(6) * t30 - t66 * t94;
t58 = qJD(3) * qJD(1);
t136 = 0.2e1 * t58;
t68 = -pkin(4) - pkin(8);
t4 = t48 + (qJD(5) - t22) * qJD(4);
t135 = t4 * t64;
t134 = t4 * t66;
t62 = -pkin(7) + qJ(2);
t133 = pkin(5) - t62;
t11 = qJD(6) * t138 + t64 * t94;
t132 = t11 * t66;
t131 = t26 * t65;
t130 = t138 * t47;
t129 = t138 * t65;
t128 = t30 * t47;
t127 = t47 * t67;
t126 = t47 * t68;
t60 = t65 ^ 2;
t70 = qJD(1) ^ 2;
t125 = t60 * t70;
t69 = qJD(4) ^ 2;
t124 = t62 * t69;
t123 = t65 * t11;
t122 = t66 * t67;
t31 = pkin(4) * t112 + qJ(5) * t119;
t120 = t69 + t70;
t117 = qJD(4) * t67;
t116 = qJD(6) * t66;
t21 = -pkin(5) * t119 + t32;
t15 = t21 + t107;
t115 = t15 * qJD(6);
t111 = t67 * qJD(2);
t110 = t67 * qJD(5);
t109 = pkin(4) * t119 - qJD(2);
t104 = t64 * t127;
t103 = t47 * t122;
t102 = t67 * t70 * t65;
t95 = t65 * t105;
t101 = pkin(4) * t94 + qJ(5) * t95 + t58;
t100 = qJD(6) * t64 * t47;
t99 = t65 * t116;
t56 = 0.2e1 * t106;
t97 = pkin(4) * t117 + t107 * t65 + qJD(3);
t77 = (qJD(4) * pkin(8) - qJD(5)) * t67;
t3 = qJD(1) * t77 + t101;
t6 = t27 + (-pkin(5) * t118 - t111) * qJD(1);
t96 = -t64 * t3 + t6 * t66;
t93 = qJD(4) * t133;
t43 = -qJD(2) + t137;
t89 = (qJD(2) - t43) * qJD(1);
t88 = -0.2e1 * t94;
t87 = -qJ(5) * t67 + t63;
t38 = t64 * t95;
t86 = -t116 * t47 + t38;
t72 = pkin(8) * t65 + t87;
t7 = qJD(1) * t72 + t109;
t8 = qJD(4) * t68 + t108;
t2 = t64 * t8 + t66 * t7;
t85 = t64 * t7 - t66 * t8;
t57 = t65 * pkin(4);
t25 = t57 + t72;
t35 = t133 * t67;
t82 = t25 * t66 + t35 * t64;
t16 = qJD(1) * t87 + t109;
t33 = t57 + t87;
t81 = qJD(1) * t33 + qJD(2) + t16;
t80 = qJD(2) + t43 + t137;
t78 = -qJD(1) * t60 + t127;
t76 = t47 * (qJD(6) * t67 + qJD(1));
t74 = t136 - t124;
t20 = t97 - t110;
t9 = -qJD(1) * t110 + t101;
t73 = -qJD(1) * t20 + t124 - t9;
t61 = t67 ^ 2;
t55 = t61 * t70;
t37 = t120 * t67;
t36 = t120 * t65;
t34 = t133 * t65;
t23 = pkin(8) * t112 + t31;
t19 = t65 * qJD(2) - t67 * t93;
t18 = -t65 * t93 - t111;
t13 = t77 + t97;
t10 = t16 * t112;
t1 = [0, 0, 0, 0, t56, qJ(2) * t56, t56, t136, t54 * qJD(2) + t43 * qJD(3) + (qJ(2) * qJD(2) + qJD(3) * t63) * qJD(1), t65 * t88, 0.2e1 * (t60 - t61) * t105, -t69 * t65, -t69 * t67, 0, t117 * t80 + t65 * t74, -t118 * t80 + t67 * t74 (-t60 - t61) * t106 - t71, -t117 * t81 + t65 * t73, t118 * t81 + t67 * t73, t16 * t20 + t9 * t33 + (-t24 * t67 - t131) * qJD(2) + t71 * t62, t64 * t123 + (t114 * t67 + t99) * t30 (t138 * t64 + t30 * t66) * t117 + (t132 - t12 * t64 + (t138 * t66 - t30 * t64) * qJD(6)) * t65, t47 * t99 + t11 * t67 + (-t30 * t65 + t64 * t78) * qJD(4), -t65 * t100 - t12 * t67 + (t66 * t78 - t129) * qJD(4), -qJD(4) * t79 (-t64 * t13 + t66 * t18) * t47 - t19 * t138 - t34 * t12 + (-t113 * t15 + t96) * t67 + (-t2 * t67 - t47 * t82) * qJD(6) + (t64 * t115 - t134 + (-(-t25 * t64 + t35 * t66) * qJD(1) + t85) * qJD(4)) * t65, -t34 * t11 + t19 * t30 + (-(qJD(6) * t35 + t13) * t47 - (qJD(6) * t8 + t3) * t67) * t66 + (-(-qJD(6) * t25 + t18) * t47 + (qJD(4) * t15 + qJD(6) * t7 - t6) * t67) * t64 + (t66 * t115 + t135 + (qJD(1) * t82 + t2) * qJD(4)) * t65; 0, 0, 0, 0, -t70, -t70 * qJ(2), -t70, 0, -qJD(1) * t54 - t58, 0, 0, 0, 0, 0, t88, 0.2e1 * t95, t55 + t125, 0.2e1 * t94, -0.2e1 * t95 (t131 + (qJD(5) + t24) * t67) * qJD(1) - t101, 0, 0, 0, 0, 0 (t103 + t129) * qJD(1) - t86, -t100 + (-t104 + (-t30 - t113) * t65) * qJD(1); 0, 0, 0, 0, 0, 0, 0, -t70, t89, 0, 0, 0, 0, 0, -t36, -t37, 0, t36, t37, -t16 * qJD(1) + t71, 0, 0, 0, 0, 0, t65 * t12 + t64 * t76 + (-t138 * t67 + t66 * t79) * qJD(4), t123 + t66 * t76 + (t30 * t67 - t64 * t79) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t55 - t125, 0, 0, 0, t67 * t89, t119 * t43 - t48 ((-t26 - t107) * t67 + (t24 + t90) * t65) * qJD(1), t10 + (t31 * t65 - t111) * qJD(1), 0.2e1 * qJD(4) * qJD(5) + t48 + (-t16 * t65 + t31 * t67) * qJD(1), -t24 * t32 - t17 * pkin(4) - t14 * qJ(5) - t16 * t31 + (-qJD(5) + t121) * t26, -t128 * t64 + t132 (-t12 - t128) * t66 + (-t11 - t130) * t64, -t100 + (-t104 + (t30 - t113) * t65) * qJD(1) (-t103 + t129) * qJD(1) + t86, t47 * t119, qJ(5) * t12 + t135 - (t21 * t66 - t23 * t64) * t47 - t108 * t138 + (-t126 * t64 + t15 * t66) * qJD(6) + (t15 * t122 + (-t113 * t68 - t85) * t65) * qJD(1), qJ(5) * t11 + t134 + (t21 * t64 + t23 * t66) * t47 + t108 * t30 + (-t126 * t66 - t15 * t64) * qJD(6) + (-t2 * t65 + (t118 * t68 - t15 * t67) * t64) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t55 - t69, qJD(4) * t26 + t10 + t17, 0, 0, 0, 0, 0, -t100 + qJD(4) * t138 + (-t113 * t65 - t104) * qJD(1), -t47 ^ 2 * t66 - qJD(4) * t30 + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30 * t138, -t138 ^ 2 + t30 ^ 2, t11 - t130, -t12 + t128, -t95, -t139 * t2 - t15 * t30 + t96, -t138 * t15 + t139 * t85 - t66 * t3 - t64 * t6;];
tauc_reg  = t1;
