% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:14
% EndTime: 2019-12-05 18:12:23
% DurationCPUTime: 1.87s
% Computational Cost: add. (2675->206), mult. (7420->283), div. (0->0), fcn. (6028->8), ass. (0->132)
t105 = sin(qJ(5));
t106 = sin(qJ(4));
t109 = cos(qJ(4));
t104 = cos(pkin(9));
t110 = cos(qJ(3));
t151 = t104 * t110;
t103 = sin(pkin(9));
t107 = sin(qJ(3));
t152 = t103 * t107;
t123 = -t151 + t152;
t83 = t123 * qJD(1);
t89 = t103 * t110 + t104 * t107;
t84 = t89 * qJD(1);
t126 = -t106 * t83 + t109 * t84;
t108 = cos(qJ(5));
t64 = t106 * t84 + t109 * t83;
t59 = t108 * t64;
t27 = t105 * t126 + t59;
t143 = -qJD(4) - qJD(5);
t99 = qJD(3) - t143;
t159 = t27 * t99;
t145 = qJD(5) * t105;
t146 = qJD(4) * t109;
t147 = qJD(4) * t106;
t140 = qJD(1) * t152;
t148 = qJD(3) * t110;
t95 = t104 * qJD(1) * t148;
t79 = -qJD(3) * t140 + t95;
t86 = t89 * qJD(3);
t80 = qJD(1) * t86;
t23 = -t106 * t80 + t109 * t79 - t83 * t146 - t84 * t147;
t24 = qJD(4) * t126 + t106 * t79 + t109 * t80;
t6 = -qJD(5) * t59 - t105 * t24 + t108 * t23 - t126 * t145;
t175 = t6 + t159;
t129 = t105 * t64 - t108 * t126;
t179 = t129 * t27;
t114 = qJD(5) * t129 - t105 * t23 - t108 * t24;
t160 = t129 * t99;
t173 = t114 - t160;
t176 = t129 ^ 2 - t27 ^ 2;
t158 = pkin(6) + qJ(2);
t93 = t158 * t103;
t90 = qJD(1) * t93;
t94 = t158 * t104;
t91 = qJD(1) * t94;
t125 = t107 * t90 - t110 * t91;
t55 = -pkin(7) * t83 - t125;
t52 = t109 * t55;
t168 = -t107 * t91 - t110 * t90;
t54 = -pkin(7) * t84 + t168;
t53 = qJD(3) * pkin(3) + t54;
t128 = -t106 * t53 - t52;
t182 = pkin(8) * t64;
t13 = -t128 - t182;
t97 = -pkin(2) * t104 - pkin(1);
t92 = t97 * qJD(1) + qJD(2);
t70 = pkin(3) * t83 + t92;
t40 = pkin(4) * t64 + t70;
t185 = t13 * t145 + t40 * t27;
t119 = t89 * qJD(2);
t118 = qJD(1) * t119;
t39 = -pkin(7) * t79 + qJD(3) * t125 - t118;
t136 = t106 * t39 - t55 * t147;
t38 = -pkin(7) * t80 - qJD(2) * t83 + t168 * qJD(3);
t164 = (qJD(4) * t53 + t38) * t109 + t136;
t2 = -pkin(8) * t24 + t164;
t137 = -t106 * t38 + t109 * t39;
t116 = qJD(4) * t128 + t137;
t3 = -pkin(8) * t23 + t116;
t172 = -t105 * t2 + t108 * t3 + t40 * t129;
t102 = qJD(3) + qJD(4);
t155 = t102 * t64;
t183 = t23 + t155;
t181 = t70 * t64;
t180 = t126 * t64;
t156 = t102 * t126;
t178 = -t24 + t156;
t177 = t126 ^ 2 - t64 ^ 2;
t174 = (-t13 * t99 - t3) * t105 + t185;
t161 = pkin(4) * t126;
t170 = pkin(8) * t126;
t169 = t70 * t126;
t165 = qJD(5) - t99;
t163 = pkin(3) * t84;
t162 = pkin(3) * t99;
t50 = t106 * t55;
t157 = t109 * t54 - t50;
t153 = t108 * t13;
t150 = t106 * t108;
t149 = t103 ^ 2 + t104 ^ 2;
t144 = qJD(1) * qJD(2);
t135 = t109 * t53 - t50;
t12 = t135 - t170;
t10 = pkin(4) * t102 + t12;
t142 = -pkin(4) * t99 - t10;
t139 = -pkin(3) * t102 - t53;
t134 = -t106 * t54 - t52;
t131 = t149 * qJD(1) ^ 2;
t130 = -t105 * t10 - t153;
t68 = t106 * t89 + t109 * t123;
t69 = -t106 * t123 + t109 * t89;
t36 = t105 * t69 + t108 * t68;
t37 = -t105 * t68 + t108 * t69;
t61 = -pkin(7) * t89 - t107 * t94 - t110 * t93;
t124 = t107 * t93 - t110 * t94;
t62 = -pkin(7) * t123 - t124;
t127 = -t106 * t61 - t109 * t62;
t74 = pkin(3) * t123 + t97;
t122 = 0.2e1 * t149 * t144;
t117 = -t93 * t148 + qJD(2) * t151 + (-qJD(2) * t103 - qJD(3) * t94) * t107;
t44 = -pkin(7) * t86 + t117;
t112 = qJD(3) * t124 - t119;
t85 = t123 * qJD(3);
t45 = pkin(7) * t85 + t112;
t120 = t106 * t45 + t109 * t44 + t61 * t146 - t62 * t147;
t115 = qJD(4) * t127 - t106 * t44 + t109 * t45;
t98 = pkin(3) * t109 + pkin(4);
t48 = pkin(4) * t68 + t74;
t46 = t163 + t161;
t35 = qJD(4) * t69 - t106 * t85 + t109 * t86;
t34 = -qJD(4) * t68 - t106 * t86 - t109 * t85;
t19 = pkin(3) * t86 + pkin(4) * t35;
t18 = pkin(3) * t80 + pkin(4) * t24;
t17 = -pkin(8) * t68 - t127;
t16 = -pkin(8) * t69 - t106 * t62 + t109 * t61;
t15 = t157 - t170;
t14 = t134 + t182;
t9 = qJD(5) * t37 + t105 * t34 + t108 * t35;
t8 = -qJD(5) * t36 - t105 * t35 + t108 * t34;
t5 = -pkin(8) * t34 + t115;
t4 = -pkin(8) * t35 + t120;
t1 = [0, 0, 0, 0, 0, t122, qJ(2) * t122, t79 * t89 - t84 * t85, -t123 * t79 - t80 * t89 + t83 * t85 - t84 * t86, -t85 * qJD(3), -t86 * qJD(3), 0, qJD(3) * t112 + t97 * t80 + t92 * t86, -qJD(3) * t117 + t97 * t79 - t92 * t85, t126 * t34 + t23 * t69, -t126 * t35 - t23 * t68 - t24 * t69 - t34 * t64, t34 * t102, -t35 * t102, 0, t74 * t24 + t70 * t35 + t115 * t102 + (t64 * t86 + t68 * t80) * pkin(3), t74 * t23 + t70 * t34 - t120 * t102 + (t126 * t86 + t69 * t80) * pkin(3), -t129 * t8 + t37 * t6, t114 * t37 + t129 * t9 - t27 * t8 - t36 * t6, t8 * t99, -t9 * t99, 0, t19 * t27 - t48 * t114 + t18 * t36 + t40 * t9 + (-t105 * t4 + t108 * t5 + (-t105 * t16 - t108 * t17) * qJD(5)) * t99, -t19 * t129 + t48 * t6 + t18 * t37 + t40 * t8 - (t105 * t5 + t108 * t4 + (-t105 * t17 + t108 * t16) * qJD(5)) * t99; 0, 0, 0, 0, 0, -t131, -qJ(2) * t131, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t84, t95 + (-t83 - t140) * qJD(3), 0, 0, 0, 0, 0, t24 + t156, t23 - t155, 0, 0, 0, 0, 0, -t114 - t160, t6 - t159; 0, 0, 0, 0, 0, 0, 0, t84 * t83, -t83 ^ 2 + t84 ^ 2, t95 + (t83 - t140) * qJD(3), 0, 0, -t92 * t84 - t118, t123 * t144 + t92 * t83, t180, t177, t183, t178, 0, -t64 * t163 - t169 - t134 * t102 + (t106 * t139 - t52) * qJD(4) + t137, -t126 * t163 + t181 + t157 * t102 + (qJD(4) * t139 - t38) * t109 - t136, -t179, t176, t175, t173, 0, -t46 * t27 - (-t105 * t15 + t108 * t14) * t99 + (-t105 * t109 - t150) * qJD(4) * t162 + ((-pkin(3) * t150 - t105 * t98) * t99 + t130) * qJD(5) + t172, t46 * t129 + (-t143 * t106 * t162 + t14 * t99 - t3) * t105 + (-qJD(5) * t10 - t2 + (-pkin(3) * t146 - qJD(5) * t98 + t15) * t99) * t108 + t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, t177, t183, t178, 0, -t102 * t128 + t116 - t169, t102 * t135 - t164 + t181, -t179, t176, t175, t173, 0, -t27 * t161 - (-t105 * t12 - t153) * t99 + (t105 * t142 - t153) * qJD(5) + t172, t129 * t161 + (qJD(5) * t142 + t12 * t99 - t2) * t108 + t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, t176, t175, t173, 0, t165 * t130 + t172, (-t165 * t10 - t2) * t108 + t174;];
tauc_reg = t1;
