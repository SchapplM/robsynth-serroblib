% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:49:21
% EndTime: 2019-03-09 01:49:25
% DurationCPUTime: 1.64s
% Computational Cost: add. (1483->259), mult. (3133->379), div. (0->0), fcn. (1985->6), ass. (0->146)
t103 = sin(qJ(4));
t159 = qJD(1) * t103;
t89 = qJD(6) + t159;
t189 = qJD(6) - t89;
t102 = sin(qJ(6));
t104 = cos(qJ(6));
t105 = cos(qJ(4));
t158 = qJD(1) * t105;
t98 = sin(pkin(9));
t145 = t98 * t158;
t99 = cos(pkin(9));
t162 = t99 * qJD(4);
t64 = t145 - t162;
t144 = t99 * t158;
t164 = qJD(4) * t98;
t66 = t144 + t164;
t115 = t102 * t64 - t104 * t66;
t188 = t115 * t89;
t100 = -pkin(7) + qJ(2);
t154 = qJD(4) * t105;
t157 = qJD(2) * t103;
t187 = t100 * t154 + t157;
t175 = t102 * t99;
t69 = t104 * t98 + t175;
t109 = t69 * qJD(1);
t110 = t69 * qJD(6);
t178 = t103 * t109 + t110;
t186 = t178 * t89;
t152 = qJD(6) * t104;
t153 = qJD(6) * t102;
t185 = t99 * t152 - t98 * t153;
t101 = pkin(1) + qJ(3);
t184 = qJD(1) * t101;
t149 = qJD(1) * qJD(4);
t135 = t103 * t149;
t128 = t104 * t135;
t129 = t102 * t135;
t8 = -qJD(6) * t115 - t98 * t128 - t99 * t129;
t23 = t102 * t66 + t104 * t64;
t183 = t23 * t89;
t150 = qJD(1) * qJD(2);
t155 = qJD(4) * t103;
t92 = qJD(1) * qJ(2) + qJD(3);
t85 = -pkin(7) * qJD(1) + t92;
t70 = t85 * t155;
t47 = -t105 * t150 + t70;
t182 = t47 * t98;
t181 = t47 * t99;
t180 = pkin(8) + qJ(5);
t118 = pkin(4) * t105 + qJ(5) * t103;
t51 = qJD(4) * t118 - t105 * qJD(5) + qJD(3);
t36 = t51 * qJD(1);
t173 = t105 * t85;
t39 = t103 * t150 + (qJD(5) + t173) * qJD(4);
t11 = t98 * t36 + t99 * t39;
t76 = pkin(4) * t103 - qJ(5) * t105 + t101;
t46 = qJD(1) * t76 - qJD(2);
t156 = qJD(4) * qJ(5);
t75 = t103 * t85;
t59 = t75 + t156;
t16 = t98 * t46 + t99 * t59;
t146 = t98 * t159;
t174 = t104 * t99;
t179 = -t102 * t146 + t159 * t174 + t185;
t171 = t105 * t99;
t71 = t118 * qJD(1);
t29 = t85 * t171 + t98 * t71;
t163 = t100 * t103;
t35 = t99 * t163 + t98 * t76;
t96 = t103 ^ 2;
t176 = -t105 ^ 2 + t96;
t172 = t105 * t98;
t107 = qJD(1) ^ 2;
t170 = t107 * t98;
t169 = t107 * t99;
t168 = t47 * t105;
t68 = t102 * t98 - t174;
t167 = qJD(4) * t68;
t166 = qJD(4) * t69;
t165 = qJD(4) * t89;
t86 = -qJD(2) + t184;
t161 = qJD(2) - t86;
t106 = qJD(4) ^ 2;
t160 = -t106 - t107;
t151 = t105 * qJD(2);
t148 = pkin(8) * t103 * t99;
t20 = t187 * t99 + t98 * t51;
t94 = 0.2e1 * t150;
t147 = 0.2e1 * qJD(3) * qJD(1);
t143 = t98 * t155;
t10 = t99 * t36 - t98 * t39;
t108 = (pkin(5) * t105 + t148) * qJD(1);
t4 = qJD(4) * t108 + t10;
t131 = pkin(8) * t143;
t5 = qJD(1) * t131 + t11;
t140 = -t102 * t5 + t104 * t4;
t139 = t179 * t89;
t137 = -t100 * t98 + pkin(5);
t136 = pkin(5) * t98 - t100;
t15 = t99 * t46 - t59 * t98;
t134 = t105 * t149;
t133 = qJD(4) * pkin(4) - qJD(5);
t132 = t161 * qJD(1);
t130 = -0.2e1 * t134;
t28 = -t85 * t172 + t99 * t71;
t53 = -t133 - t173;
t127 = t133 + t53;
t126 = -t10 * t99 - t11 * t98;
t125 = -t10 * t98 + t11 * t99;
t124 = t102 * t4 + t104 * t5;
t6 = pkin(5) * t159 - pkin(8) * t66 + t15;
t9 = -pkin(8) * t64 + t16;
t123 = t102 * t9 - t104 * t6;
t2 = t102 * t6 + t104 * t9;
t122 = -t15 * t99 - t16 * t98;
t121 = t15 * t98 - t16 * t99;
t120 = -t64 * t99 + t66 * t98;
t119 = qJD(2) + t86 + t184;
t61 = t99 * t76;
t22 = -pkin(8) * t171 + t137 * t103 + t61;
t26 = -pkin(8) * t172 + t35;
t117 = -t102 * t26 + t104 * t22;
t116 = t102 * t22 + t104 * t26;
t114 = -t100 * t106 + t147;
t83 = t180 * t99;
t113 = qJD(5) * t98 + qJD(6) * t83 + t108 + t28;
t82 = t180 * t98;
t112 = pkin(8) * t146 - qJD(5) * t99 + qJD(6) * t82 + t29;
t111 = qJD(1) * t68;
t7 = -t99 * t128 + t98 * t129 - t64 * t152 - t66 * t153;
t91 = -pkin(5) * t99 - pkin(4);
t63 = t136 * t105;
t50 = t68 * t105;
t49 = t69 * t105;
t45 = -pkin(5) * t146 + t75;
t42 = -t136 * t155 - t151;
t41 = t99 * t51;
t34 = -t98 * t163 + t61;
t30 = t70 + (-pkin(5) * t143 - t151) * qJD(1);
t27 = pkin(5) * t64 + t53;
t19 = -t187 * t98 + t41;
t18 = -t104 * t143 + t105 * t185 - t155 * t175;
t17 = -t105 * t110 + t68 * t155;
t13 = t131 + t20;
t12 = -t98 * t157 + t41 + (t137 * t105 + t148) * qJD(4);
t1 = [0, 0, 0, 0, t94, qJ(2) * t94, t94, t147, t92 * qJD(2) + t86 * qJD(3) + (qJ(2) * qJD(2) + qJD(3) * t101) * qJD(1), t103 * t130, 0.2e1 * t176 * t149, -t106 * t103, -t106 * t105, 0, t103 * t114 + t119 * t154, t105 * t114 - t119 * t155 (-qJD(2) * t64 + t182) * t105 + (qJD(1) * t19 + t10) * t103 + ((qJD(1) * t34 + t15) * t105 + (-t53 * t98 + (t64 + t145) * t100) * t103) * qJD(4) (-qJD(2) * t66 + t181) * t105 + (-qJD(1) * t20 - t11) * t103 + ((-qJD(1) * t35 - t16) * t105 + (-t53 * t99 + (t66 + t144) * t100) * t103) * qJD(4), -t19 * t66 - t20 * t64 + t126 * t105 + ((t34 * t99 + t35 * t98) * qJD(1) - t122) * t155, -t53 * t151 + t10 * t34 + t11 * t35 + t15 * t19 + t16 * t20 + (t53 * t155 - t168) * t100, -t115 * t17 - t50 * t7, t115 * t18 - t17 * t23 - t49 * t7 + t50 * t8, t7 * t103 + t17 * t89 + (-qJD(1) * t50 - t115) * t154, -t8 * t103 - t18 * t89 + (-qJD(1) * t49 - t23) * t154 (t89 + t159) * t154 (-t102 * t13 + t104 * t12) * t89 + t140 * t103 + t42 * t23 + t63 * t8 + t30 * t49 + t27 * t18 + (-t103 * t2 - t116 * t89) * qJD(6) + (qJD(1) * t117 - t123) * t154 -(t102 * t12 + t104 * t13) * t89 - t124 * t103 - t42 * t115 + t63 * t7 - t30 * t50 + t27 * t17 + (t103 * t123 - t117 * t89) * qJD(6) + (-qJD(1) * t116 - t2) * t154; 0, 0, 0, 0, -t107, -t107 * qJ(2), -t107, 0 (-qJD(3) - t92) * qJD(1), 0, 0, 0, 0, 0, t130, 0.2e1 * t135, t96 * t170 + (t64 - t162) * t158, t96 * t169 + (t66 + t164) * t158 ((-t98 ^ 2 - t99 ^ 2) * qJD(4) - t120) * t159 (t103 * t121 + t105 * t53) * qJD(1) + t126, 0, 0, 0, 0, 0, t186 + (t23 + t167) * t158, t139 + (-t115 + t166) * t158; 0, 0, 0, 0, 0, 0, 0, -t107, t132, 0, 0, 0, 0, 0, t160 * t103, t160 * t105 (-t169 + (t64 - t145) * qJD(4)) * t103 (t170 + (t66 - t144) * qJD(4)) * t103, t120 * t154 + (t64 * t98 + t66 * t99) * qJD(1), -t168 + t125 * t103 + t122 * qJD(1) + (t103 * t53 - t105 * t121) * qJD(4), 0, 0, 0, 0, 0, t89 * t111 + (-t69 * t165 - t8) * t105 + (t68 * t89 * qJD(6) + (-t158 * t69 + t23) * qJD(4)) * t103, t89 * t109 + (t68 * t165 - t7) * t105 + (t89 * t110 + (t105 * t111 - t115) * qJD(4)) * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, t105 * t107 * t103, -t176 * t107, 0, 0, 0, t105 * t132, -t161 * t159, -t64 * t75 - t181 + ((-t98 * t156 - t15) * t105 + (t127 * t98 - t28) * t103) * qJD(1), -t66 * t75 + t182 + ((-t99 * t156 + t16) * t105 + (t127 * t99 + t29) * t103) * qJD(1), t28 * t66 + t29 * t64 + (-qJD(5) * t64 - t15 * t159 + t11) * t99 + (qJD(5) * t66 - t16 * t159 - t10) * t98, -t47 * pkin(4) + t125 * qJ(5) - t121 * qJD(5) - t15 * t28 - t16 * t29 - t53 * t75, -t115 * t179 + t7 * t69, t115 * t178 - t179 * t23 - t7 * t68 - t69 * t8, t139 + (t115 + t166) * t158, -t186 + (t23 - t167) * t158, -t89 * t158, -t45 * t23 + t30 * t68 + t91 * t8 + (t102 * t112 - t104 * t113) * t89 + t178 * t27 + ((-t102 * t83 - t104 * t82) * qJD(4) + t123) * t158, t45 * t115 + t30 * t69 + t91 * t7 + (t102 * t113 + t104 * t112) * t89 + t179 * t27 + (-(-t102 * t82 + t104 * t83) * qJD(4) + t2) * t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t66 - t164) * t159 (-t64 - t162) * t159, -t64 ^ 2 - t66 ^ 2, t15 * t66 + t16 * t64 + t47, 0, 0, 0, 0, 0, t8 - t188, t7 - t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115 * t23, t115 ^ 2 - t23 ^ 2, t7 + t183, -t8 - t188, t134, t27 * t115 - t189 * t2 + t140, t189 * t123 + t27 * t23 - t124;];
tauc_reg  = t1;
