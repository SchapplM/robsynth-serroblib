% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PRPR7
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PRPR7_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR7_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:56
% EndTime: 2019-12-31 16:25:58
% DurationCPUTime: 1.31s
% Computational Cost: add. (1922->213), mult. (3553->298), div. (0->0), fcn. (2202->6), ass. (0->148)
t188 = cos(qJ(2));
t190 = qJD(2) ^ 2;
t186 = sin(qJ(2));
t205 = t186 * qJDD(2);
t159 = t188 * t190 + t205;
t181 = sin(pkin(6));
t182 = cos(pkin(6));
t163 = t181 * g(1) - t182 * g(2);
t120 = -pkin(4) * t159 + t188 * t163;
t238 = t181 * t120;
t237 = t182 * t120;
t203 = t188 * qJDD(2);
t160 = t186 * t190 - t203;
t119 = pkin(4) * t160 - t186 * t163;
t236 = t181 * t119;
t235 = t182 * t119;
t164 = t182 * g(1) + t181 * g(2);
t179 = g(3) - qJDD(1);
t136 = -t186 * t164 + t188 * t179;
t195 = -pkin(1) * t160 + qJ(1) * t159 - t136;
t146 = t182 * t163;
t116 = -t181 * t164 + t146;
t234 = pkin(1) * t159 + qJ(1) * t160;
t233 = pkin(2) + pkin(5);
t185 = sin(qJ(4));
t177 = t185 ^ 2;
t228 = t177 * t190;
t187 = cos(qJ(4));
t178 = t187 ^ 2;
t227 = t178 * t190;
t226 = t181 * t159;
t225 = t181 * t160;
t224 = t181 * t163;
t222 = t181 * t179;
t221 = t182 * t159;
t220 = t182 * t160;
t219 = t182 * t179;
t137 = -t188 * t164 - t186 * t179;
t191 = (2 * qJD(3) * qJD(2)) + t137;
t208 = qJDD(2) * qJ(3);
t112 = -t190 * pkin(2) + t191 + t208;
t106 = -t190 * pkin(5) + t112;
t218 = t185 * t106;
t202 = t187 * t190 * t185;
t165 = qJDD(4) + t202;
t217 = t185 * t165;
t210 = t177 + t178;
t158 = t210 * qJDD(2);
t216 = t186 * t158;
t215 = t187 * t106;
t214 = t187 * t165;
t166 = qJDD(4) - t202;
t213 = t187 * t166;
t212 = t188 * t158;
t209 = qJD(2) * qJD(4);
t207 = t182 * qJDD(2);
t206 = t185 * qJDD(2);
t204 = t187 * qJDD(2);
t201 = t185 * t209;
t200 = t187 * t209;
t161 = t210 * t190;
t114 = -t186 * t161 + t212;
t115 = -t188 * t161 - t216;
t180 = qJDD(2) * pkin(2);
t113 = -t190 * qJ(3) + qJDD(3) + t136 - t180;
t109 = -qJDD(2) * pkin(5) + t113;
t89 = -t187 * t109 - t185 * t163;
t90 = t185 * t109 - t187 * t163;
t70 = t185 * t90 - t187 * t89;
t199 = -pkin(1) * t114 + qJ(1) * t115 + qJ(3) * t161 - t233 * t158 + t70;
t198 = t191 + 0.2e1 * t208 + t234;
t197 = -qJDD(3) + 0.2e1 * t180 + t195;
t196 = t137 + t234;
t85 = t188 * t112 + t186 * t113;
t94 = t186 * t136 + t188 * t137;
t117 = -t182 * t164 - t224;
t189 = qJD(4) ^ 2;
t194 = -t189 - t227;
t193 = t186 * t202;
t192 = t188 * t202;
t71 = t185 * t89 + t187 * t90;
t83 = t186 * t112 - t188 * t113;
t93 = t188 * t136 - t186 * t137;
t172 = t181 * qJDD(2);
t169 = t189 - t227;
t168 = -t189 - t228;
t167 = -t189 + t228;
t162 = (-t177 + t178) * t190;
t156 = pkin(1) * t163;
t155 = -0.2e1 * t201 + t204;
t154 = -t201 + t204;
t153 = -t200 - t206;
t152 = 0.2e1 * t200 + t206;
t151 = t185 * t166;
t150 = t210 * t209;
t135 = t188 * qJDD(4) - t186 * t150;
t134 = -t185 * t154 - t178 * t209;
t133 = -t187 * t153 - t177 * t209;
t131 = -t185 * t194 - t214;
t130 = t185 * t169 - t213;
t129 = (-t154 + t201) * t187;
t128 = t187 * t168 - t151;
t127 = -t187 * t167 + t217;
t126 = t187 * t194 - t217;
t125 = -t187 * t169 - t151;
t124 = t185 * t168 + t213;
t123 = -t185 * t167 - t214;
t122 = (t153 - t200) * t185;
t111 = t187 * t152 + t185 * t155;
t110 = t185 * t152 - t187 * t155;
t104 = -t186 * t133 - t192;
t103 = -t186 * t134 + t192;
t102 = -t186 * t125 + t187 * t203;
t101 = -t186 * t123 - t185 * t203;
t98 = t186 * t126 + t188 * t155;
t97 = t186 * t124 + t188 * t152;
t96 = -t188 * t126 + t186 * t155;
t95 = -t188 * t124 + t186 * t152;
t91 = -t186 * t110 + t188 * t162;
t87 = t182 * t94 - t224;
t86 = t181 * t94 + t146;
t82 = t181 * t131 + t182 * t98;
t81 = t181 * t128 + t182 * t97;
t80 = -t182 * t131 + t181 * t98;
t79 = -t182 * t128 + t181 * t97;
t78 = t182 * t85 - t224;
t77 = t181 * t85 + t146;
t76 = -pkin(4) * t83 + (-pkin(2) * t186 + qJ(3) * t188) * t163;
t75 = pkin(3) * t126 - qJ(3) * t131 - t90;
t74 = pkin(3) * t124 - qJ(3) * t128 - t89;
t73 = pkin(3) * t152 - t233 * t128 + t215;
t72 = pkin(3) * t155 - t233 * t131 - t218;
t69 = -pkin(3) * t161 - t71;
t68 = -pkin(1) * t83 + pkin(2) * t113 - qJ(3) * t112;
t67 = t188 * t106 + t186 * t70;
t66 = t186 * t106 - t188 * t70;
t65 = -pkin(1) * t96 - qJ(3) * t155 + t233 * t126 - t215;
t64 = -pkin(1) * t95 - qJ(3) * t152 + t233 * t124 - t218;
t63 = -pkin(3) * t212 - pkin(4) * t114 - t186 * t69;
t61 = pkin(3) * t70 - qJ(3) * t71;
t60 = pkin(3) * t106 - t233 * t71;
t59 = -pkin(4) * t96 - t186 * t72 + t188 * t75;
t58 = -pkin(4) * t95 - t186 * t73 + t188 * t74;
t57 = t181 * t71 + t182 * t67;
t56 = t181 * t67 - t182 * t71;
t55 = -pkin(1) * t66 - qJ(3) * t106 + t233 * t70;
t54 = -pkin(4) * t66 - t186 * t60 + t188 * t61;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, 0, 0, 0, 0, 0, 0, -t221, t220, 0, t87, 0, 0, 0, 0, 0, 0, 0, t221, -t220, t78, 0, 0, 0, 0, 0, 0, t81, t82, t182 * t115, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, 0, 0, 0, 0, 0, 0, -t226, t225, 0, t86, 0, 0, 0, 0, 0, 0, 0, t226, -t225, t77, 0, 0, 0, 0, 0, 0, t79, t80, t181 * t115, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, 0, 0, 0, 0, 0, 0, -t160, -t159, 0, -t93, 0, 0, 0, 0, 0, 0, 0, t160, t159, t83, 0, 0, 0, 0, 0, 0, t95, t96, t114, t66; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t222, -t219, -t116, -qJ(1) * t116, 0, 0, -t220, 0, -t221, t172, t195 * t181 + t235, -t196 * t181 - t237, t182 * t93, -qJ(1) * t86 - (pkin(1) * t181 - pkin(4) * t182) * t93, t172, t220, t221, 0, 0, 0, -t182 * t83, -t197 * t181 - t235, t198 * t181 + t237, -qJ(1) * t77 - t181 * t68 + t182 * t76, t182 * t103 - t181 * t129, -t181 * t111 + t182 * t91, t182 * t102 - t181 * t130, t182 * t104 - t181 * t122, t182 * t101 - t181 * t127, t182 * t135, -qJ(1) * t79 - t181 * t64 + t182 * t58, -qJ(1) * t80 - t181 * t65 + t182 * t59, -t199 * t181 + t182 * t63, -qJ(1) * t56 - t181 * t55 + t182 * t54; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t219, -t222, t117, qJ(1) * t117, 0, 0, -t225, 0, -t226, -t207, -t195 * t182 + t236, t196 * t182 - t238, t181 * t93, qJ(1) * t87 - (-pkin(1) * t182 - pkin(4) * t181) * t93, -t207, t225, t226, 0, 0, 0, -t181 * t83, t197 * t182 - t236, -t198 * t182 + t238, qJ(1) * t78 + t181 * t76 + t182 * t68, t181 * t103 + t182 * t129, t182 * t111 + t181 * t91, t181 * t102 + t182 * t130, t181 * t104 + t182 * t122, t181 * t101 + t182 * t127, t181 * t135, qJ(1) * t81 + t181 * t58 + t182 * t64, qJ(1) * t82 + t181 * t59 + t182 * t65, t181 * t63 + t199 * t182, qJ(1) * t57 + t181 * t54 + t182 * t55; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t163, t164, 0, 0, 0, 0, t159, 0, -t160, 0, t120, t119, t94, pkin(4) * t94 + t156, 0, -t159, t160, 0, 0, 0, t85, -t120, -t119, pkin(4) * t85 + t156 + (pkin(2) * t188 + qJ(3) * t186) * t163, t188 * t134 + t193, t188 * t110 + t186 * t162, t188 * t125 + t186 * t204, t188 * t133 - t193, t188 * t123 - t185 * t205, t186 * qJDD(4) + t188 * t150, -pkin(1) * t128 + pkin(4) * t97 + t186 * t74 + t188 * t73, -pkin(1) * t131 + pkin(4) * t98 + t186 * t75 + t188 * t72, -pkin(3) * t216 + pkin(4) * t115 + t188 * t69, -pkin(1) * t71 + pkin(4) * t67 + t186 * t61 + t188 * t60;];
tauB_reg = t1;
