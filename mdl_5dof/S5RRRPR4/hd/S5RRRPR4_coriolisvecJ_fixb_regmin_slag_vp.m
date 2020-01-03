% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:39
% EndTime: 2019-12-31 21:11:44
% DurationCPUTime: 1.41s
% Computational Cost: add. (1073->195), mult. (1846->264), div. (0->0), fcn. (1039->6), ass. (0->138)
t92 = sin(qJ(5));
t93 = sin(qJ(3));
t95 = cos(qJ(5));
t96 = cos(qJ(3));
t54 = -t96 * t92 + t93 * t95;
t139 = qJD(3) * t96;
t140 = qJD(3) * t93;
t19 = t54 * qJD(5) + t92 * t139 - t95 * t140;
t53 = t93 * t92 + t96 * t95;
t87 = qJD(1) + qJD(2);
t37 = t53 * t87;
t105 = t53 * qJD(5);
t127 = t87 * t139;
t128 = t87 * t140;
t7 = -t87 * t105 + t95 * t127 + t92 * t128;
t86 = qJD(3) - qJD(5);
t174 = -t37 * t86 + t7;
t142 = pkin(1) * qJD(2);
t122 = qJD(1) * t142;
t97 = cos(qJ(2));
t116 = t97 * t122;
t68 = t96 * t116;
t88 = qJD(3) * qJD(4);
t146 = t68 + t88;
t143 = pkin(1) * qJD(1);
t94 = sin(qJ(2));
t130 = t94 * t143;
t59 = t87 * pkin(7) + t130;
t21 = -t59 * t140 + t146;
t67 = t93 * t116;
t26 = t59 * t139 + t67;
t172 = t21 * t96 + t26 * t93;
t39 = t54 * t87;
t8 = t19 * t87;
t171 = t39 * t86 + t8;
t82 = t93 * qJ(4);
t84 = t96 * pkin(3);
t170 = t82 + t84;
t90 = t93 ^ 2;
t91 = t96 ^ 2;
t169 = (t90 + t91) * t87;
t47 = t93 * t59;
t167 = qJD(4) + t47;
t166 = qJD(5) + t86;
t165 = pkin(3) + pkin(4);
t164 = pkin(7) - pkin(8);
t120 = pkin(2) + t82;
t129 = t97 * t143;
t12 = t129 + (t165 * t96 + t120) * t87;
t141 = qJ(4) * t96;
t107 = -t165 * t93 + t141;
t117 = t94 * t122;
t80 = t93 * qJD(4);
t9 = -t117 + (t107 * qJD(3) + t80) * t87;
t163 = t12 * t19 + t9 * t53;
t20 = t53 * qJD(3) - t105;
t162 = t12 * t20 + t9 * t54;
t99 = qJD(3) ^ 2;
t161 = pkin(7) * t99;
t76 = t94 * pkin(1) + pkin(7);
t160 = -pkin(8) + t76;
t158 = t39 * t37;
t156 = t86 * t97;
t155 = t87 * t93;
t154 = t87 * t96;
t48 = t96 * t59;
t151 = t99 * t93;
t81 = t99 * t96;
t60 = -t87 * pkin(2) - t129;
t150 = t93 * t117 + t60 * t139;
t138 = qJD(3) * t97;
t125 = t93 * t138;
t149 = t125 * t143 + t130 * t154;
t145 = t90 - t91;
t89 = qJD(3) * qJ(4);
t41 = t48 + t89;
t137 = t41 * qJD(3);
t136 = -qJD(1) - t87;
t135 = -pkin(8) * t155 + t167;
t134 = pkin(2) + t170;
t85 = t87 ^ 2;
t133 = t93 * t85 * t96;
t42 = pkin(3) * t140 - t96 * t89 - t80;
t132 = t97 * t142;
t131 = t94 * t142;
t72 = t164 * t96;
t123 = t37 ^ 2 - t39 ^ 2;
t77 = -t97 * pkin(1) - pkin(2);
t52 = t160 * t96;
t113 = pkin(3) * t93 - t141;
t13 = t117 + (t113 * qJD(3) - t80) * t87;
t121 = -t13 - t161;
t119 = -qJD(3) * pkin(3) + qJD(4);
t118 = t86 ^ 2;
t32 = -pkin(8) * t154 + t48;
t115 = (-qJD(2) + t87) * t143;
t114 = t136 * t142;
t18 = -t165 * qJD(3) + t135;
t25 = t32 + t89;
t112 = t95 * t18 - t92 * t25;
t111 = -t92 * t18 - t95 * t25;
t35 = t119 + t47;
t110 = t35 * t93 + t41 * t96;
t49 = t77 - t170;
t34 = t131 + t42;
t109 = -t34 * t87 - t76 * t99 - t13;
t108 = t49 * t87 - t132;
t33 = -pkin(4) * t140 - t42;
t106 = -t93 * t137 + t35 * t139 + t172;
t104 = t96 * t138 - t94 * t155;
t10 = (pkin(8) * t87 - t59) * t140 + t146;
t11 = -pkin(8) * t127 + t26;
t102 = -t95 * t10 - t92 * t11 + t12 * t37;
t101 = -t92 * t10 + t95 * t11 - t12 * t39;
t100 = (t35 * t96 - t41 * t93) * qJD(3) + t172;
t83 = t96 * pkin(4);
t78 = pkin(8) * t140;
t71 = t164 * t93;
t58 = 0.2e1 * t93 * t127;
t56 = qJD(3) * t72;
t55 = -pkin(7) * t140 + t78;
t51 = t160 * t93;
t50 = t83 + t134;
t45 = t60 * t140;
t43 = t113 * t87;
t40 = -t49 + t83;
t36 = -0.2e1 * t145 * t87 * qJD(3);
t30 = t107 * t87;
t28 = qJD(3) * t52 + t93 * t132;
t27 = t96 * t132 - t76 * t140 + t78;
t24 = -t129 + (-t120 - t84) * t87;
t23 = t33 - t131;
t17 = t24 * t140;
t15 = t20 * t86;
t14 = t19 * t86;
t2 = t39 * t20 + t7 * t54;
t1 = -t39 * t19 - t20 * t37 - t7 * t53 - t54 * t8;
t3 = [0, 0, 0, 0, t94 * t114, t97 * t114, t58, t36, t81, -t151, 0, t77 * t128 - t76 * t81 + t45 + (t136 * t96 * t94 - t125) * t142, -t104 * t142 + t77 * t127 + t76 * t151 + t150, t108 * t140 + t109 * t96 + t17, t132 * t169 + t106, t109 * t93 + (-t108 - t24) * t139, t100 * t76 + t110 * t132 + t13 * t49 + t24 * t34, t2, t1, -t15, t14, 0, t23 * t37 + t40 * t8 - (-t92 * t27 + t95 * t28 + (-t51 * t92 - t52 * t95) * qJD(5)) * t86 + t163, t23 * t39 + t40 * t7 + (t95 * t27 + t92 * t28 + (t51 * t95 - t52 * t92) * qJD(5)) * t86 + t162; 0, 0, 0, 0, t94 * t115, t97 * t115, t58, t36, t81, -t151, 0, -pkin(2) * t128 + t45 + (-t117 - t161) * t96 + t149, -pkin(2) * t127 + pkin(7) * t151 + t104 * t143 + t150, -t134 * t128 + t17 + (-t42 * t87 + t121) * t96 + t149, -t129 * t169 + t106, (t134 * t87 - t129 - t24) * t139 + ((-t42 + t130) * t87 + t121) * t93, -t13 * t134 + t24 * t42 + (-t110 * t97 - t24 * t94) * t143 + t100 * pkin(7), t2, t1, -t15, t14, 0, t33 * t37 + t50 * t8 - (-t92 * t55 + t95 * t56 + (-t71 * t92 - t72 * t95) * qJD(5)) * t86 + (t156 * t54 + t94 * t37) * t143 + t163, t33 * t39 + t50 * t7 + (t95 * t55 + t92 * t56 + (t71 * t95 - t72 * t92) * qJD(5)) * t86 + (-t156 * t53 + t94 * t39) * t143 + t162; 0, 0, 0, 0, 0, 0, -t133, t145 * t85, 0, 0, 0, -t60 * t155 - t67, -t60 * t154 - t68, -t67 + (-t24 * t93 + t43 * t96) * t87, ((t41 - t89) * t93 + (t119 - t35) * t96) * t87, t68 + 0.2e1 * t88 + (t24 * t96 + t43 * t93) * t87, -t26 * pkin(3) + t21 * qJ(4) + t167 * t41 - t24 * t43 - t35 * t48, -t158, t123, -t174, t171, 0, -t30 * t37 + (t135 * t92 + t95 * t32) * t86 + (-(-qJ(4) * t95 + t165 * t92) * t86 - t111) * qJD(5) - t101, -t30 * t39 + (t135 * t95 - t92 * t32) * t86 + ((-qJ(4) * t92 - t165 * t95) * t86 + t112) * qJD(5) - t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, 0, -t90 * t85 - t99, t155 * t24 - t137 + t26, 0, 0, 0, 0, 0, -t118 * t92 - t155 * t37, -t118 * t95 - t155 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, -t123, t174, -t171, 0, t166 * t111 + t101, -t166 * t112 + t102;];
tauc_reg = t3;
