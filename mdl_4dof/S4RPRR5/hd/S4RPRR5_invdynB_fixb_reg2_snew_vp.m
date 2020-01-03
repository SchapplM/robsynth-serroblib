% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPRR5
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPRR5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR5_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:41
% EndTime: 2019-12-31 16:51:44
% DurationCPUTime: 1.44s
% Computational Cost: add. (3353->215), mult. (5021->291), div. (0->0), fcn. (2300->6), ass. (0->139)
t211 = sin(qJ(3));
t206 = -qJD(1) + qJD(3);
t204 = t206 ^ 2;
t205 = qJDD(1) - qJDD(3);
t214 = cos(qJ(3));
t231 = t214 * t205;
t222 = t211 * t204 + t231;
t158 = pkin(5) * t222 + t211 * g(3);
t212 = sin(qJ(1));
t215 = cos(qJ(1));
t235 = t211 * t205;
t226 = -t214 * t204 + t235;
t244 = t212 * t226 + t215 * t222;
t245 = -pkin(5) * t226 + t214 * g(3);
t250 = -pkin(4) * t244 + t158 * t215 - t212 * t245;
t136 = t212 * t222 - t215 * t226;
t249 = -pkin(4) * t136 + t158 * t212 + t215 * t245;
t217 = qJD(1) ^ 2;
t207 = qJDD(1) * qJ(2);
t190 = t215 * g(1) + t212 * g(2);
t221 = (2 * qJD(2) * qJD(1)) - t190;
t220 = t207 + t221;
t243 = pkin(1) + pkin(2);
t157 = -t243 * t217 + t220;
t189 = t212 * g(1) - t215 * g(2);
t223 = -qJDD(2) + t189;
t219 = -t217 * qJ(2) - t223;
t218 = -t243 * qJDD(1) + t219;
t122 = t211 * t157 - t214 * t218;
t123 = t214 * t157 + t211 * t218;
t103 = t214 * t122 - t211 * t123;
t104 = t211 * t122 + t214 * t123;
t86 = t103 * t215 + t212 * t104;
t246 = t103 * t212 - t215 * t104;
t241 = qJDD(1) * pkin(1);
t210 = sin(qJ(4));
t208 = t210 ^ 2;
t240 = t208 * t204;
t116 = t205 * pkin(3) - t204 * pkin(6) + t122;
t239 = t210 * t116;
t213 = cos(qJ(4));
t188 = t213 * t204 * t210;
t180 = qJDD(4) + t188;
t238 = t210 * t180;
t181 = qJDD(4) - t188;
t237 = t210 * t181;
t236 = t210 * t205;
t234 = t213 * t116;
t233 = t213 * t180;
t232 = t213 * t181;
t192 = t213 * t205;
t117 = -t204 * pkin(3) - t205 * pkin(6) + t123;
t112 = t210 * g(3) + t213 * t117;
t209 = t213 ^ 2;
t230 = t208 + t209;
t229 = qJD(4) * t206;
t228 = t210 * t229;
t227 = t213 * t229;
t111 = -t213 * g(3) + t210 * t117;
t162 = -t217 * pkin(1) + t220;
t164 = -t219 + t241;
t129 = t215 * t162 - t212 * t164;
t152 = -t212 * t189 - t215 * t190;
t225 = t211 * t188;
t224 = t214 * t188;
t182 = t212 * qJDD(1) + t215 * t217;
t167 = -pkin(4) * t182 + t215 * g(3);
t183 = t215 * qJDD(1) - t212 * t217;
t166 = pkin(4) * t183 + t212 * g(3);
t93 = t213 * t111 - t210 * t112;
t94 = t210 * t111 + t213 * t112;
t128 = t212 * t162 + t215 * t164;
t151 = t215 * t189 - t212 * t190;
t216 = qJD(4) ^ 2;
t194 = t209 * t204;
t187 = -t194 - t216;
t186 = t194 - t216;
t185 = -t216 - t240;
t184 = t216 - t240;
t178 = t194 - t240;
t177 = t194 + t240;
t172 = t230 * t205;
t171 = -t192 - 0.2e1 * t228;
t170 = -t192 - t228;
t169 = t227 - t236;
t168 = 0.2e1 * t227 - t236;
t165 = t230 * t229;
t150 = t211 * qJDD(4) + t214 * t165;
t149 = t214 * qJDD(4) - t211 * t165;
t148 = t213 * t169 - t208 * t229;
t147 = -t210 * t170 - t209 * t229;
t146 = -t210 * t185 - t232;
t145 = -t210 * t184 + t233;
t144 = t213 * t187 - t238;
t143 = t213 * t186 - t237;
t142 = t213 * t185 - t237;
t141 = t210 * t187 + t233;
t138 = -t214 * t172 - t211 * t177;
t135 = -t211 * t172 + t214 * t177;
t134 = -t210 * t168 + t213 * t171;
t133 = t214 * t145 - t210 * t235;
t132 = t214 * t143 - t211 * t192;
t131 = -t211 * t145 - t210 * t231;
t130 = -t211 * t143 - t213 * t231;
t127 = t214 * t148 - t225;
t126 = t214 * t147 + t225;
t125 = -t211 * t148 - t224;
t124 = -t211 * t147 + t224;
t121 = t214 * t146 + t211 * t168;
t120 = t214 * t144 - t211 * t171;
t119 = t211 * t146 - t214 * t168;
t118 = t211 * t144 + t214 * t171;
t115 = t214 * t134 - t211 * t178;
t114 = -t211 * t134 - t214 * t178;
t110 = t212 * t135 + t215 * t138;
t109 = -t215 * t135 + t212 * t138;
t108 = -pkin(6) * t142 + t234;
t107 = -pkin(6) * t141 + t239;
t106 = -pkin(3) * t142 + t112;
t105 = -pkin(3) * t141 + t111;
t100 = t212 * t119 + t215 * t121;
t99 = t212 * t118 + t215 * t120;
t98 = -t215 * t119 + t212 * t121;
t97 = -t215 * t118 + t212 * t120;
t96 = pkin(5) * t103 + qJ(2) * g(3);
t95 = -pkin(5) * t104 + t243 * g(3);
t91 = -pkin(5) * t135 + t214 * t93;
t90 = -pkin(5) * t138 - t211 * t93;
t89 = t211 * t116 + t214 * t94;
t88 = -t214 * t116 + t211 * t94;
t85 = -pkin(5) * t119 + qJ(2) * t142 - t211 * t106 + t214 * t108;
t84 = -pkin(5) * t118 + qJ(2) * t141 - t211 * t105 + t214 * t107;
t83 = -pkin(5) * t121 - t214 * t106 - t211 * t108 + t243 * t142;
t82 = -pkin(5) * t120 - t214 * t105 - t211 * t107 + t243 * t141;
t81 = t212 * t88 + t215 * t89;
t80 = t212 * t89 - t215 * t88;
t79 = -pkin(5) * t88 - (pkin(3) * t211 - pkin(6) * t214 + qJ(2)) * t93;
t78 = -pkin(5) * t89 - (pkin(3) * t214 + pkin(6) * t211 + t243) * t93;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t182, -t183, 0, t152, 0, 0, 0, 0, 0, 0, -t182, 0, t183, t129, 0, 0, 0, 0, 0, 0, -t136, t244, 0, -t246, 0, 0, 0, 0, 0, 0, t99, t100, t110, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t183, -t182, 0, t151, 0, 0, 0, 0, 0, 0, t183, 0, t182, t128, 0, 0, 0, 0, 0, 0, t244, t136, 0, t86, 0, 0, 0, 0, 0, 0, t97, t98, t109, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t141, -t142, 0, t93; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t183, 0, -t182, 0, -t166, -t167, -t151, -pkin(4) * t151, 0, t183, 0, 0, t182, 0, -t166, -t128, t167, -pkin(4) * t128 + (-pkin(1) * t212 + qJ(2) * t215) * g(3), 0, 0, -t244, 0, -t136, 0, t250, t249, t86, -pkin(4) * t86 - t212 * t95 + t215 * t96, -t125 * t212 + t127 * t215, -t114 * t212 + t115 * t215, -t131 * t212 + t133 * t215, -t124 * t212 + t126 * t215, -t130 * t212 + t132 * t215, -t149 * t212 + t150 * t215, -pkin(4) * t97 - t212 * t82 + t215 * t84, -pkin(4) * t98 - t212 * t83 + t215 * t85, -pkin(4) * t109 - t212 * t90 + t215 * t91, -pkin(4) * t80 - t212 * t78 + t215 * t79; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t182, 0, t183, 0, t167, -t166, t152, pkin(4) * t152, 0, t182, 0, 0, -t183, 0, t167, t129, t166, pkin(4) * t129 + (pkin(1) * t215 + qJ(2) * t212) * g(3), 0, 0, -t136, 0, t244, 0, t249, -t250, t246, -pkin(4) * t246 + t212 * t96 + t215 * t95, t125 * t215 + t127 * t212, t114 * t215 + t115 * t212, t131 * t215 + t133 * t212, t124 * t215 + t126 * t212, t130 * t215 + t132 * t212, t149 * t215 + t150 * t212, pkin(4) * t99 + t212 * t84 + t215 * t82, pkin(4) * t100 + t212 * t85 + t215 * t83, pkin(4) * t110 + t212 * t91 + t215 * t90, pkin(4) * t81 + t212 * t79 + t215 * t78; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t189, t190, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t223 + 0.2e1 * t241, 0, 0.2e1 * t207 + t221, pkin(1) * t164 + qJ(2) * t162, 0, 0, 0, 0, 0, t205, qJ(2) * t226 + t222 * t243 + t122, qJ(2) * t222 - t226 * t243 + t123, 0, qJ(2) * t104 + t103 * t243, (-t169 - t227) * t210, -t168 * t213 - t171 * t210, -t184 * t213 - t238, (-t170 + t228) * t213, -t186 * t210 - t232, 0, -pkin(3) * t171 - pkin(6) * t144 + qJ(2) * t120 - t243 * t118 + t234, pkin(3) * t168 - pkin(6) * t146 + qJ(2) * t121 - t243 * t119 - t239, -pkin(3) * t177 + pkin(6) * t172 + qJ(2) * t138 - t243 * t135 - t94, pkin(3) * t116 - pkin(6) * t94 + qJ(2) * t89 - t243 * t88;];
tauB_reg = t1;
