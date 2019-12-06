% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:18
% EndTime: 2019-12-05 16:52:24
% DurationCPUTime: 1.40s
% Computational Cost: add. (1768->216), mult. (4469->280), div. (0->0), fcn. (3037->6), ass. (0->131)
t86 = sin(qJ(4));
t89 = cos(qJ(3));
t155 = t86 * t89;
t164 = cos(qJ(4));
t87 = sin(qJ(3));
t61 = t164 * t87 + t155;
t83 = qJD(3) + qJD(4);
t171 = t83 * t61;
t29 = t171 * qJD(2);
t55 = t61 * qJD(2);
t15 = -t55 * t83 + t29;
t156 = t86 * t87;
t103 = t164 * t89 - t156;
t90 = cos(qJ(2));
t100 = t90 * t103;
t167 = -pkin(7) - pkin(6);
t66 = t167 * t87;
t67 = t167 * t89;
t104 = t164 * t66 + t86 * t67;
t126 = qJD(3) * t167;
t63 = t87 * t126;
t150 = -qJD(1) * t100 + t104 * qJD(4) + t126 * t155 + t164 * t63;
t133 = qJD(2) * qJD(3);
t170 = -0.2e1 * t133;
t169 = t150 * t83;
t88 = sin(qJ(2));
t136 = t88 * qJD(1);
t139 = qJD(3) * t87;
t102 = pkin(3) * t139 - t136;
t168 = t55 ^ 2;
t135 = t90 * qJD(1);
t142 = qJD(2) * pkin(6);
t70 = t136 + t142;
t38 = -t70 * t139 + (-pkin(7) * t139 + t89 * t135) * qJD(2);
t138 = qJD(3) * t89;
t39 = -t70 * t138 + (-pkin(7) * t138 - t87 * t135) * qJD(2);
t123 = -t164 * t39 + t86 * t38;
t119 = pkin(7) * qJD(2) + t70;
t50 = t119 * t89;
t127 = t164 * t50;
t49 = t119 * t87;
t45 = qJD(3) * pkin(3) - t49;
t24 = t86 * t45 + t127;
t3 = t24 * qJD(4) + t123;
t166 = t3 * t104;
t51 = t61 * t88;
t165 = t3 * t51;
t122 = qJD(2) * t164;
t112 = t89 * t122;
t141 = qJD(2) * t87;
t128 = t86 * t141;
t53 = -t112 + t128;
t81 = -t89 * pkin(3) - pkin(2);
t58 = t81 * qJD(2) - t135;
t20 = t53 * pkin(4) - t55 * qJ(5) + t58;
t163 = t20 * t55;
t157 = t86 * t50;
t23 = t164 * t45 - t157;
t162 = t23 * t83;
t161 = t55 * t53;
t159 = t58 * t55;
t158 = t83 * t171;
t154 = t87 * t88;
t91 = qJD(3) ^ 2;
t153 = t91 * t87;
t152 = t91 * t89;
t92 = qJD(2) ^ 2;
t151 = t92 * t90;
t121 = t164 * qJD(3);
t113 = t89 * t121;
t42 = -t164 * t67 + t86 * t66;
t149 = t42 * qJD(4) - t167 * t113 - t61 * t135 + t86 * t63;
t120 = t164 * qJD(4);
t26 = -t164 * t49 - t157;
t148 = pkin(3) * t120 + qJD(5) - t26;
t147 = t83 * t112;
t124 = t87 * t133;
t59 = pkin(3) * t124 + qJD(2) * t136;
t84 = t87 ^ 2;
t85 = t89 ^ 2;
t146 = t84 - t85;
t145 = t84 + t85;
t144 = t91 + t92;
t143 = qJD(2) * pkin(2);
t140 = qJD(2) * t88;
t137 = qJD(4) * t86;
t134 = qJD(5) - t23;
t131 = t87 * t92 * t89;
t129 = pkin(3) * t141;
t21 = t53 ^ 2 - t168;
t125 = t149 * t83;
t110 = t83 * t156;
t36 = -t89 * t120 + t110 - t113;
t118 = -pkin(4) * t171 - t36 * qJ(5) + t61 * qJD(5) - t102;
t117 = -t45 * t120 + t50 * t137 - t164 * t38 - t86 * t39;
t71 = -t135 - t143;
t116 = -t71 - t135;
t115 = t90 * t170;
t114 = t89 * t124;
t25 = -t86 * t49 + t127;
t111 = pkin(3) * t137 - t25;
t30 = t55 * pkin(4) + t53 * qJ(5);
t109 = -t103 * t29 + t171 * t53;
t108 = qJD(2) * t116;
t11 = t87 * t90 * t122 + (qJD(2) * t86 * t90 + (t121 + t120) * t88) * t89 - t83 * t86 * t154;
t107 = -t11 * t83 + t53 * t140 - t90 * t29;
t106 = -t20 * t53 - t117;
t105 = t58 * t53 + t117;
t10 = qJD(2) * t100 - t171 * t88;
t28 = qJD(2) * t110 - t147;
t52 = t103 * t88;
t101 = -t10 * t53 + t11 * t55 - t51 * t28 - t52 * t29;
t99 = qJD(3) * (-t116 - t143);
t98 = t10 * t83 - t55 * t140 - t90 * t28;
t97 = -t103 * t28 - t171 * t55 - t61 * t29 + t36 * t53;
t96 = t104 * t28 + t149 * t55 - t150 * t53 - t42 * t29 + t3 * t61;
t95 = -t123 + (-qJD(4) + t83) * t24;
t93 = t25 * t83 + (-t127 + (-pkin(3) * t83 - t45) * t86) * qJD(4) - t123;
t82 = t83 * qJD(5);
t80 = -t164 * pkin(3) - pkin(4);
t77 = t86 * pkin(3) + qJ(5);
t32 = t36 * t83;
t31 = -pkin(4) * t103 - t61 * qJ(5) + t81;
t27 = t30 + t129;
t19 = t83 * qJ(5) + t24;
t18 = -t83 * pkin(4) + t134;
t14 = t147 + (-t128 + t53) * t83;
t5 = -t28 * t61 - t55 * t36;
t4 = t29 * pkin(4) + t28 * qJ(5) - t55 * qJD(5) + t59;
t1 = t82 - t117;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92 * t88, -t151, 0, 0, 0, 0, 0, 0, 0, 0, -t144 * t89 * t88 + t87 * t115, t89 * t115 + t144 * t154, t145 * t151, (t71 * t88 + (-t136 + (t70 + t136) * t145) * t90) * qJD(2), 0, 0, 0, 0, 0, 0, t107, -t98, t101, t24 * t10 - t23 * t11 - t117 * t52 + t140 * t58 - t59 * t90 + t165, 0, 0, 0, 0, 0, 0, t107, t101, t98, t1 * t52 + t19 * t10 + t18 * t11 + t140 * t20 - t4 * t90 + t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t114, t146 * t170, t152, -0.2e1 * t114, -t153, 0, -pkin(6) * t152 + t87 * t99, pkin(6) * t153 + t89 * t99, 0, ((-t71 - t143) * t88 + (t142 - t70) * t90 * t145) * qJD(1), t5, t97, -t32, t109, -t158, 0, t102 * t53 - t103 * t59 + t171 * t58 + t81 * t29 - t125, t102 * t55 - t81 * t28 - t58 * t36 + t59 * t61 - t169, -t103 * t117 - t171 * t24 + t23 * t36 + t96, t102 * t58 - t117 * t42 - t149 * t23 + t150 * t24 + t59 * t81 - t166, t5, -t32, -t97, 0, t158, t109, -t103 * t4 - t118 * t53 + t171 * t20 + t31 * t29 - t125, t1 * t103 - t171 * t19 - t18 * t36 + t96, t118 * t55 + t20 * t36 + t31 * t28 - t4 * t61 + t169, t1 * t42 - t118 * t20 + t149 * t18 + t150 * t19 + t4 * t31 - t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, t146 * t92, 0, t131, 0, 0, t87 * t108, t89 * t108, 0, 0, t161, -t21, t14, -t161, -t15, 0, -t129 * t53 - t159 + t93, t26 * t83 + (-t120 * t83 - t141 * t55) * pkin(3) + t105, (t24 - t25) * t55 + (-t23 + t26) * t53 + (t164 * t28 - t29 * t86 + (-t164 * t53 + t55 * t86) * qJD(4)) * pkin(3), t23 * t25 - t24 * t26 + (-t58 * t141 - t164 * t3 - t117 * t86 + (t164 * t24 - t23 * t86) * qJD(4)) * pkin(3), t161, t14, t21, 0, t15, -t161, -t27 * t53 - t163 + t93, -t80 * t28 - t77 * t29 + (t111 + t19) * t55 + (t18 - t148) * t53, t148 * t83 + t27 * t55 + t106 + t82, t1 * t77 + t111 * t18 + t148 * t19 - t20 * t27 + t3 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, -t21, t14, -t161, -t15, 0, t95 - t159, t105 + t162, 0, 0, t161, t14, t21, 0, t15, -t161, -t30 * t53 - t163 + t95, pkin(4) * t28 - t29 * qJ(5) + (t19 - t24) * t55 + (t18 - t134) * t53, t30 * t55 + t106 - t162 + 0.2e1 * t82, -t3 * pkin(4) + t1 * qJ(5) + t134 * t19 - t18 * t24 - t20 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, t14, -t83 ^ 2 - t168, -t19 * t83 + t163 + t3;];
tauc_reg = t2;
