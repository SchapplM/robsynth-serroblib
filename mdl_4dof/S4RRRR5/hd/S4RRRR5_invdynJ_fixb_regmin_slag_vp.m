% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [4x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:15
% EndTime: 2019-12-31 17:28:21
% DurationCPUTime: 2.26s
% Computational Cost: add. (1621->294), mult. (3728->428), div. (0->0), fcn. (2620->10), ass. (0->157)
t115 = sin(qJ(2));
t116 = sin(qJ(1));
t120 = cos(qJ(1));
t141 = g(1) * t120 + g(2) * t116;
t119 = cos(qJ(2));
t199 = g(3) * t119;
t127 = t141 * t115 - t199;
t164 = t115 * qJDD(1);
t101 = pkin(5) * t164;
t165 = qJD(1) * qJD(2);
t153 = t119 * t165;
t53 = -qJDD(2) * pkin(2) + pkin(5) * t153 + t101;
t166 = t119 * qJD(1);
t96 = -qJD(3) + t166;
t214 = qJD(3) * pkin(6) * t96 + t127 - t53;
t113 = sin(qJ(4));
t117 = cos(qJ(4));
t114 = sin(qJ(3));
t175 = qJD(1) * t115;
t154 = t114 * t175;
t118 = cos(qJ(3));
t167 = t118 * qJD(2);
t71 = t154 - t167;
t168 = t114 * qJD(2);
t73 = t118 * t175 + t168;
t138 = t113 * t71 - t117 * t73;
t24 = t113 * t73 + t117 * t71;
t213 = t138 * t24;
t212 = -qJD(3) * t175 + qJDD(2);
t211 = t138 ^ 2 - t24 ^ 2;
t169 = qJD(4) * t117;
t170 = qJD(4) * t113;
t21 = qJD(3) * t167 + (t153 + t164) * t118 + t212 * t114;
t22 = ((qJD(3) + t166) * qJD(2) + t164) * t114 - t212 * t118;
t4 = -t113 * t22 + t117 * t21 - t71 * t169 - t73 * t170;
t91 = -qJD(4) + t96;
t210 = -t24 * t91 + t4;
t112 = qJ(3) + qJ(4);
t109 = cos(t112);
t81 = -t119 * pkin(2) - t115 * pkin(6) - pkin(1);
t65 = t81 * qJD(1);
t103 = pkin(5) * t166;
t85 = qJD(2) * pkin(6) + t103;
t33 = t114 * t65 + t118 * t85;
t16 = -t71 * pkin(7) + t33;
t14 = t16 * t170;
t200 = g(3) * t115;
t84 = -qJD(2) * pkin(2) + pkin(5) * t175;
t38 = t71 * pkin(3) + t84;
t108 = sin(t112);
t180 = t120 * t108;
t182 = t116 * t119;
t45 = -t109 * t182 + t180;
t179 = t120 * t109;
t47 = t116 * t108 + t119 * t179;
t209 = g(1) * t47 - g(2) * t45 + t109 * t200 + t38 * t24 + t14;
t142 = pkin(2) * t115 - pkin(6) * t119;
t79 = t142 * qJD(2);
t36 = qJD(1) * t79 + qJDD(1) * t81;
t29 = t118 * t36;
t106 = t119 * qJDD(1);
t151 = t115 * t165;
t205 = -t151 + t106;
t52 = t205 * pkin(5) + qJDD(2) * pkin(6);
t68 = qJDD(3) - t205;
t2 = t68 * pkin(3) - t21 * pkin(7) - t33 * qJD(3) - t114 * t52 + t29;
t171 = qJD(3) * t118;
t172 = qJD(3) * t114;
t131 = t114 * t36 + t118 * t52 + t65 * t171 - t85 * t172;
t3 = -t22 * pkin(7) + t131;
t158 = -t113 * t3 + t117 * t2;
t44 = t108 * t182 + t179;
t46 = t116 * t109 - t119 * t180;
t32 = -t114 * t85 + t118 * t65;
t15 = -t73 * pkin(7) + t32;
t13 = -t96 * pkin(3) + t15;
t190 = t117 * t16;
t7 = t113 * t13 + t190;
t208 = -g(1) * t46 + g(2) * t44 - t7 * qJD(4) + t108 * t200 + t38 * t138 + t158;
t124 = t138 * qJD(4) - t113 * t21 - t117 * t22;
t207 = t138 * t91 + t124;
t75 = t113 * t118 + t117 * t114;
t49 = t75 * t115;
t204 = qJD(3) + qJD(4);
t202 = pkin(6) + pkin(7);
t201 = pkin(3) * t114;
t198 = t71 * t96;
t197 = t73 * t96;
t74 = t113 * t114 - t117 * t118;
t130 = t74 * t119;
t196 = qJD(1) * t130 - t204 * t74;
t195 = (-t166 + t204) * t75;
t194 = t114 * t79 + t81 * t171;
t193 = t115 * pkin(5) * t168 + t118 * t79;
t76 = t142 * qJD(1);
t192 = pkin(5) * t154 + t118 * t76;
t181 = t118 * t119;
t97 = pkin(5) * t181;
t191 = t114 * t81 + t97;
t189 = t118 * t73;
t188 = t21 * t114;
t187 = qJD(2) * t71;
t186 = qJD(2) * t73;
t185 = t114 * t115;
t184 = t114 * t119;
t183 = t115 * t118;
t178 = t120 * t114;
t177 = t120 * t118;
t110 = t115 ^ 2;
t176 = -t119 ^ 2 + t110;
t174 = qJD(2) * t115;
t173 = qJD(2) * t119;
t162 = t96 * t167;
t161 = t96 * t172;
t160 = t96 * t171;
t157 = qJD(3) * t202;
t156 = t119 * t167;
t155 = t119 * t168;
t150 = qJD(4) * t13 + t3;
t148 = -qJD(3) * t65 - t52;
t146 = pkin(3) * t172 - t166 * t201 - t103;
t137 = pkin(3) * t115 - pkin(7) * t181;
t87 = t202 * t118;
t145 = qJD(1) * t137 + qJD(4) * t87 + t118 * t157 + t192;
t60 = t114 * t76;
t86 = t202 * t114;
t144 = -qJD(4) * t86 - t60 - (-pkin(5) * t183 - pkin(7) * t184) * qJD(1) - t114 * t157;
t143 = t85 * t171 - t29;
t140 = g(1) * t116 - g(2) * t120;
t139 = -pkin(6) * t68 + t84 * qJD(3);
t135 = t114 * t68 - t160;
t134 = t118 * t68 + t161;
t132 = -0.2e1 * pkin(1) * t165 - pkin(5) * qJDD(2);
t122 = qJD(1) ^ 2;
t129 = pkin(1) * t122 + t141;
t128 = t115 * t171 + t155;
t121 = qJD(2) ^ 2;
t125 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t121 + t140;
t100 = -t118 * pkin(3) - pkin(2);
t80 = (pkin(5) + t201) * t115;
t70 = t118 * t81;
t64 = qJDD(4) + t68;
t59 = t116 * t114 + t119 * t177;
t58 = t116 * t118 - t119 * t178;
t57 = -t116 * t181 + t178;
t56 = t114 * t182 + t177;
t50 = t74 * t115;
t39 = pkin(3) * t128 + pkin(5) * t173;
t37 = -pkin(7) * t185 + t191;
t31 = -pkin(7) * t183 + t70 + (-pkin(5) * t114 - pkin(3)) * t119;
t12 = t22 * pkin(3) + t53;
t11 = (t204 * t183 + t155) * t117 + (-t204 * t185 + t156) * t113;
t10 = -qJD(2) * t130 - t204 * t49;
t9 = -t128 * pkin(7) + (-t115 * t167 - t119 * t172) * pkin(5) + t194;
t8 = t137 * qJD(2) + (-t97 + (pkin(7) * t115 - t81) * t114) * qJD(3) + t193;
t6 = -t113 * t16 + t117 * t13;
t1 = [qJDD(1), t140, t141, t110 * qJDD(1) + 0.2e1 * t119 * t151, 0.2e1 * t115 * t106 - 0.2e1 * t176 * t165, qJDD(2) * t115 + t121 * t119, qJDD(2) * t119 - t121 * t115, 0, t115 * t132 + t119 * t125, -t115 * t125 + t119 * t132, t73 * t156 + (t21 * t118 - t73 * t172) * t115, (-t114 * t73 - t118 * t71) * t173 + (-t188 - t118 * t22 + (t114 * t71 - t189) * qJD(3)) * t115, (-t21 - t162) * t119 + (t134 + t186) * t115, (t96 * t168 + t22) * t119 + (-t135 - t187) * t115, -t68 * t119 - t96 * t174, -(-t81 * t172 + t193) * t96 + t70 * t68 - g(1) * t57 - g(2) * t59 + ((t160 + t187) * pkin(5) + (-pkin(5) * t68 + qJD(2) * t84 - t148) * t114 + t143) * t119 + (pkin(5) * t22 + t32 * qJD(2) + t53 * t114 + t84 * t171) * t115, t194 * t96 - t191 * t68 - g(1) * t56 - g(2) * t58 + (t84 * t167 + (-t161 + t186) * pkin(5) + t131) * t119 + (-t84 * t172 - t33 * qJD(2) + t53 * t118 + (t21 - t162) * pkin(5)) * t115, -t10 * t138 - t4 * t50, -t10 * t24 + t11 * t138 - t124 * t50 - t4 * t49, -t10 * t91 - t4 * t119 - t138 * t174 - t50 * t64, t11 * t91 - t119 * t124 - t24 * t174 - t49 * t64, -t64 * t119 - t91 * t174, -(-t113 * t9 + t117 * t8) * t91 + (-t113 * t37 + t117 * t31) * t64 - t158 * t119 + t6 * t174 + t39 * t24 - t80 * t124 + t12 * t49 + t38 * t11 - g(1) * t45 - g(2) * t47 + (-(-t113 * t31 - t117 * t37) * t91 + t7 * t119) * qJD(4), -t7 * t174 - g(1) * t44 - g(2) * t46 + t38 * t10 - t14 * t119 - t12 * t50 - t39 * t138 + t80 * t4 + ((-qJD(4) * t37 + t8) * t91 - t31 * t64 + t2 * t119) * t113 + ((qJD(4) * t31 + t9) * t91 - t37 * t64 + t150 * t119) * t117; 0, 0, 0, -t115 * t122 * t119, t176 * t122, t164, t106, qJDD(2), t115 * t129 - t101 - t199, t200 + (-pkin(5) * qJDD(1) + t129) * t119, -t96 * t189 + t188, (t21 + t198) * t118 + (-t22 + t197) * t114, (-t115 * t73 + t96 * t181) * qJD(1) + t135, (t115 * t71 - t96 * t184) * qJD(1) + t134, t96 * t175, -pkin(2) * t22 + t192 * t96 + t139 * t114 + (-t32 * t115 + (-pkin(5) * t71 - t114 * t84) * t119) * qJD(1) + t214 * t118, -pkin(2) * t21 - t60 * t96 + t139 * t118 + (-t84 * t181 + t33 * t115 + (-t119 * t73 + t96 * t183) * pkin(5)) * qJD(1) - t214 * t114, -t138 * t196 + t4 * t75, t124 * t75 + t138 * t195 - t196 * t24 - t4 * t74, t138 * t175 - t196 * t91 + t75 * t64, t24 * t175 + t195 * t91 - t74 * t64, t91 * t175, (-t113 * t87 - t117 * t86) * t64 - t100 * t124 + t12 * t74 - t6 * t175 + (t113 * t144 + t117 * t145) * t91 + t195 * t38 + t146 * t24 + t127 * t109, -(-t113 * t86 + t117 * t87) * t64 + t100 * t4 + t12 * t75 + t7 * t175 + (-t113 * t145 + t117 * t144) * t91 + t196 * t38 - t146 * t138 - t127 * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t71, -t71 ^ 2 + t73 ^ 2, t21 - t198, -t197 - t22, t68, -g(1) * t58 + g(2) * t56 - t33 * t96 - t84 * t73 + (t148 + t200) * t114 - t143, g(1) * t59 - g(2) * t57 + g(3) * t183 - t32 * t96 + t84 * t71 - t131, -t213, t211, t210, t207, t64, (-t113 * t15 - t190) * t91 + (t117 * t64 + t170 * t91 - t73 * t24) * pkin(3) + t208, (t16 * t91 - t2) * t113 + (-t15 * t91 - t150) * t117 + (-t113 * t64 + t138 * t73 + t169 * t91) * pkin(3) + t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213, t211, t210, t207, t64, -t7 * t91 + t208, -t113 * t2 - t117 * t150 - t6 * t91 + t209;];
tau_reg = t1;
