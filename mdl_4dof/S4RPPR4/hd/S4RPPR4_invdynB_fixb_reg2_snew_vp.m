% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPPR4
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPPR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:55
% EndTime: 2019-12-31 16:38:58
% DurationCPUTime: 1.86s
% Computational Cost: add. (2568->220), mult. (4776->300), div. (0->0), fcn. (2684->6), ass. (0->144)
t198 = sin(pkin(6));
t199 = cos(pkin(6));
t207 = qJD(1) ^ 2;
t172 = -qJDD(1) * t199 + t198 * t207;
t196 = g(3) - qJDD(2);
t152 = qJ(2) * t172 - t196 * t198;
t203 = sin(qJ(1));
t205 = cos(qJ(1));
t171 = qJDD(1) * t198 + t199 * t207;
t216 = t171 * t205 - t172 * t203;
t220 = -qJ(2) * t171 + t196 * t199;
t268 = -pkin(4) * t216 + t203 * t152 + t205 * t220;
t180 = g(1) * t203 - g(2) * t205;
t165 = qJDD(1) * pkin(1) + t180;
t181 = g(1) * t205 + g(2) * t203;
t166 = -pkin(1) * t207 - t181;
t126 = -t165 * t199 + t166 * t198;
t127 = t165 * t198 + t166 * t199;
t217 = t126 * t198 + t127 * t199;
t99 = t126 * t199 - t127 * t198;
t249 = t205 * t99;
t265 = -t203 * t217 + t249;
t250 = t203 * t99;
t78 = t205 * t217 + t250;
t214 = t171 * t203 + t172 * t205;
t253 = pkin(4) * t214 + t205 * t152 - t203 * t220;
t208 = (2 * qJD(3) * qJD(1)) + t127;
t226 = qJDD(1) * qJ(3);
t115 = -pkin(2) * t207 + t208 + t226;
t197 = qJDD(1) * pkin(2);
t116 = -qJ(3) * t207 + qJDD(3) + t126 - t197;
t218 = t115 * t199 + t116 * t198;
t93 = t115 * t198 - t116 * t199;
t72 = -t203 * t93 + t205 * t218;
t71 = t203 * t218 + t205 * t93;
t252 = pkin(2) + pkin(5);
t251 = pkin(1) * t171;
t202 = sin(qJ(4));
t194 = t202 ^ 2;
t247 = t194 * t207;
t204 = cos(qJ(4));
t195 = t204 ^ 2;
t246 = t195 * t207;
t228 = t194 + t195;
t173 = t228 * qJDD(1);
t245 = t198 * t173;
t244 = t199 * t173;
t109 = -pkin(5) * t207 + t115;
t242 = t202 * t109;
t223 = t202 * t207 * t204;
t178 = qJDD(4) + t223;
t241 = t202 * t178;
t179 = qJDD(4) - t223;
t240 = t202 * t179;
t235 = t204 * t109;
t234 = t204 * t178;
t233 = t204 * t179;
t227 = qJD(1) * qJD(4);
t225 = t202 * qJDD(1);
t224 = t204 * qJDD(1);
t222 = t202 * t227;
t221 = t204 * t227;
t114 = -qJDD(1) * pkin(5) + t116;
t103 = -t114 * t204 - t196 * t202;
t140 = -t180 * t203 - t181 * t205;
t212 = t198 * t223;
t211 = t199 * t223;
t175 = qJDD(1) * t205 - t203 * t207;
t210 = -pkin(4) * t175 - g(3) * t203;
t104 = t114 * t202 - t196 * t204;
t80 = -t204 * t103 + t202 * t104;
t81 = t103 * t202 + t104 * t204;
t139 = t180 * t205 - t181 * t203;
t209 = -pkin(1) * t172 - t126;
t206 = qJD(4) ^ 2;
t185 = -t206 - t246;
t184 = t206 - t246;
t183 = -t206 - t247;
t182 = -t206 + t247;
t177 = (-t194 + t195) * t207;
t176 = t228 * t207;
t174 = qJDD(1) * t203 + t205 * t207;
t170 = -0.2e1 * t222 + t224;
t169 = -t222 + t224;
t168 = -t221 - t225;
t167 = 0.2e1 * t221 + t225;
t163 = t228 * t227;
t157 = -pkin(4) * t174 + g(3) * t205;
t150 = -t169 * t202 - t195 * t227;
t149 = -t168 * t204 - t194 * t227;
t148 = qJDD(4) * t199 - t163 * t198;
t147 = qJDD(4) * t198 + t163 * t199;
t146 = -t185 * t202 - t234;
t145 = t183 * t204 - t240;
t144 = t185 * t204 - t241;
t143 = -t184 * t204 - t240;
t142 = t183 * t202 + t233;
t141 = -t182 * t202 - t234;
t138 = -t176 * t199 - t245;
t137 = -t176 * t198 + t244;
t128 = t167 * t202 - t170 * t204;
t124 = -t149 * t198 - t211;
t123 = -t150 * t198 + t211;
t122 = t149 * t199 - t212;
t121 = t150 * t199 + t212;
t120 = -t143 * t198 + t199 * t224;
t119 = -t141 * t198 - t199 * t225;
t118 = t143 * t199 + t198 * t224;
t117 = t141 * t199 - t198 * t225;
t113 = t144 * t198 + t170 * t199;
t112 = t142 * t198 + t167 * t199;
t111 = -t144 * t199 + t170 * t198;
t110 = -t142 * t199 + t167 * t198;
t106 = -t128 * t198 + t177 * t199;
t105 = t128 * t199 + t177 * t198;
t102 = -t137 * t203 + t138 * t205;
t101 = t137 * t205 + t138 * t203;
t96 = pkin(1) * t196 + qJ(2) * t217;
t91 = -t111 * t203 + t113 * t205;
t90 = -t110 * t203 + t112 * t205;
t89 = t111 * t205 + t113 * t203;
t88 = t110 * t205 + t112 * t203;
t87 = pkin(3) * t144 - qJ(3) * t146 - t104;
t86 = pkin(3) * t142 - qJ(3) * t145 - t103;
t85 = -qJ(2) * t93 + (-pkin(2) * t198 + qJ(3) * t199) * t196;
t84 = qJ(2) * t218 + (pkin(2) * t199 + qJ(3) * t198 + pkin(1)) * t196;
t83 = pkin(3) * t167 - t145 * t252 + t235;
t82 = pkin(3) * t170 - t146 * t252 - t242;
t79 = -pkin(3) * t176 - t81;
t76 = t109 * t199 + t198 * t80;
t75 = t109 * t198 - t199 * t80;
t74 = -pkin(3) * t244 - qJ(2) * t137 - t198 * t79;
t73 = -pkin(3) * t245 + qJ(2) * t138 + t199 * t79;
t70 = pkin(3) * t80 - qJ(3) * t81;
t69 = -qJ(2) * t111 - t198 * t82 + t199 * t87;
t68 = -qJ(2) * t110 - t198 * t83 + t199 * t86;
t67 = pkin(3) * t109 - t252 * t81;
t66 = -pkin(1) * t146 + qJ(2) * t113 + t198 * t87 + t199 * t82;
t65 = -pkin(1) * t145 + qJ(2) * t112 + t198 * t86 + t199 * t83;
t64 = -t203 * t75 + t205 * t76;
t63 = t203 * t76 + t205 * t75;
t62 = -qJ(2) * t75 - t198 * t67 + t199 * t70;
t61 = -pkin(1) * t81 + qJ(2) * t76 + t198 * t70 + t199 * t67;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t174, -t175, 0, t140, 0, 0, 0, 0, 0, 0, -t216, t214, 0, t78, 0, 0, 0, 0, 0, 0, 0, t216, -t214, t72, 0, 0, 0, 0, 0, 0, t90, t91, t102, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t175, -t174, 0, t139, 0, 0, 0, 0, 0, 0, -t214, -t216, 0, -t265, 0, 0, 0, 0, 0, 0, 0, t214, t216, t71, 0, 0, 0, 0, 0, 0, t88, t89, t101, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t196, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t196, 0, 0, 0, 0, 0, 0, t145, t146, 0, t81; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t175, 0, -t174, 0, t210, -t157, -t139, -pkin(4) * t139, 0, 0, -t214, 0, -t216, 0, t253, -t268, t265, pkin(4) * t265 + qJ(2) * t249 - t203 * t96, 0, t214, t216, 0, 0, 0, -t71, -t253, t268, -pkin(4) * t71 - t203 * t84 + t205 * t85, -t121 * t203 + t123 * t205, -t105 * t203 + t106 * t205, -t118 * t203 + t120 * t205, -t122 * t203 + t124 * t205, -t117 * t203 + t119 * t205, -t147 * t203 + t148 * t205, -pkin(4) * t88 - t203 * t65 + t205 * t68, -pkin(4) * t89 - t203 * t66 + t205 * t69, -pkin(4) * t101 - t203 * t73 + t205 * t74, -pkin(4) * t63 - t203 * t61 + t205 * t62; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t174, 0, t175, 0, t157, t210, t140, pkin(4) * t140, 0, 0, t216, 0, -t214, 0, t268, t253, t78, pkin(4) * t78 + qJ(2) * t250 + t205 * t96, 0, -t216, t214, 0, 0, 0, t72, -t268, -t253, pkin(4) * t72 + t203 * t85 + t205 * t84, t121 * t205 + t123 * t203, t105 * t205 + t106 * t203, t118 * t205 + t120 * t203, t122 * t205 + t124 * t203, t117 * t205 + t119 * t203, t147 * t205 + t148 * t203, pkin(4) * t90 + t203 * t68 + t205 * t65, pkin(4) * t91 + t203 * t69 + t205 * t66, pkin(4) * t102 + t203 * t74 + t205 * t73, pkin(4) * t64 + t203 * t62 + t205 * t61; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t180, t181, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t209, -t127 - t251, 0, -pkin(1) * t99, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - 0.2e1 * t197 - t209, t208 + 0.2e1 * t226 + t251, pkin(1) * t93 - pkin(2) * t116 + qJ(3) * t115, (t169 - t222) * t204, -t167 * t204 - t170 * t202, -t184 * t202 + t233, (-t168 + t221) * t202, t182 * t204 - t241, 0, pkin(1) * t110 + qJ(3) * t167 - t142 * t252 + t242, pkin(1) * t111 + qJ(3) * t170 - t144 * t252 + t235, pkin(1) * t137 - qJ(3) * t176 + t173 * t252 - t80, pkin(1) * t75 + qJ(3) * t109 - t252 * t80;];
tauB_reg = t1;
