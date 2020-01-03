% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PRPR4
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
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PRPR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:01
% EndTime: 2019-12-31 16:22:05
% DurationCPUTime: 1.82s
% Computational Cost: add. (2181->212), mult. (4172->292), div. (0->0), fcn. (2660->6), ass. (0->142)
t189 = sin(qJ(2));
t193 = qJD(2) ^ 2;
t191 = cos(qJ(2));
t210 = t191 * qJDD(2);
t161 = -t189 * t193 + t210;
t182 = g(3) - qJDD(1);
t143 = pkin(4) * t161 + t189 * t182;
t184 = sin(pkin(6));
t185 = cos(pkin(6));
t212 = t189 * qJDD(2);
t162 = t191 * t193 + t212;
t200 = t184 * t161 + t185 * t162;
t206 = -pkin(4) * t162 + t191 * t182;
t258 = -qJ(1) * t200 - t184 * t143 + t185 * t206;
t165 = t184 * g(1) - t185 * g(2);
t166 = t185 * g(1) + t184 * g(2);
t132 = t189 * t165 - t191 * t166;
t198 = t191 * t165 + t189 * t166;
t203 = t191 * t132 - t189 * t198;
t93 = -t189 * t132 - t191 * t198;
t238 = t185 * t93;
t255 = -t184 * t203 + t238;
t239 = t184 * t93;
t70 = t185 * t203 + t239;
t202 = t185 * t161 - t184 * t162;
t243 = qJ(1) * t202 + t185 * t143 + t184 * t206;
t194 = (2 * qJD(3) * qJD(2)) + t132;
t214 = qJDD(2) * qJ(3);
t107 = -t193 * pkin(2) + t194 + t214;
t183 = qJDD(2) * pkin(2);
t116 = -t193 * qJ(3) + qJDD(3) - t183 - t198;
t204 = t191 * t107 + t189 * t116;
t85 = t189 * t107 - t191 * t116;
t64 = -t184 * t85 + t185 * t204;
t63 = t184 * t204 + t185 * t85;
t242 = pkin(5) + pkin(2);
t241 = pkin(1) * t162;
t188 = sin(qJ(4));
t180 = t188 ^ 2;
t237 = t180 * t193;
t190 = cos(qJ(4));
t181 = t190 ^ 2;
t236 = t181 * t193;
t231 = t184 * t182;
t226 = t185 * t182;
t105 = -t193 * pkin(5) + t107;
t225 = t188 * t105;
t209 = t190 * t193 * t188;
t167 = qJDD(4) + t209;
t224 = t188 * t167;
t168 = qJDD(4) - t209;
t223 = t188 * t168;
t216 = t180 + t181;
t160 = t216 * qJDD(2);
t222 = t189 * t160;
t221 = t190 * t105;
t220 = t190 * t167;
t219 = t190 * t168;
t218 = t191 * t160;
t215 = qJD(2) * qJD(4);
t213 = t188 * qJDD(2);
t211 = t190 * qJDD(2);
t208 = t188 * t215;
t207 = t190 * t215;
t106 = -qJDD(2) * pkin(5) + t116;
t95 = -t190 * t106 - t188 * t182;
t130 = -t184 * t165 - t185 * t166;
t197 = t189 * t209;
t196 = t191 * t209;
t96 = t188 * t106 - t190 * t182;
t74 = t188 * t96 - t190 * t95;
t75 = t188 * t95 + t190 * t96;
t129 = t185 * t165 - t184 * t166;
t195 = pkin(1) * t161 + t198;
t192 = qJD(4) ^ 2;
t172 = -t192 - t236;
t171 = t192 - t236;
t170 = -t192 - t237;
t169 = -t192 + t237;
t164 = (-t180 + t181) * t193;
t163 = t216 * t193;
t159 = -0.2e1 * t208 + t211;
t158 = -t208 + t211;
t157 = -t207 - t213;
t156 = 0.2e1 * t207 + t213;
t155 = t216 * t215;
t142 = t191 * qJDD(4) - t189 * t155;
t141 = t189 * qJDD(4) + t191 * t155;
t140 = -t188 * t158 - t181 * t215;
t139 = -t190 * t157 - t180 * t215;
t138 = -t188 * t172 - t220;
t137 = t190 * t170 - t223;
t136 = t190 * t172 - t224;
t135 = -t190 * t171 - t223;
t134 = t188 * t170 + t219;
t133 = -t188 * t169 - t220;
t128 = -t191 * t163 - t222;
t127 = -t189 * t163 + t218;
t117 = t188 * t156 - t190 * t159;
t115 = -t189 * t139 - t196;
t114 = -t189 * t140 + t196;
t113 = t191 * t139 - t197;
t112 = t191 * t140 + t197;
t111 = -t189 * t135 + t190 * t210;
t110 = -t189 * t133 - t188 * t210;
t109 = t191 * t135 + t189 * t211;
t108 = t191 * t133 - t188 * t212;
t102 = t189 * t136 + t191 * t159;
t101 = t189 * t134 + t191 * t156;
t100 = -t191 * t136 + t189 * t159;
t99 = -t191 * t134 + t189 * t156;
t98 = -t189 * t117 + t191 * t164;
t97 = t191 * t117 + t189 * t164;
t90 = -t184 * t127 + t185 * t128;
t89 = t185 * t127 + t184 * t128;
t88 = pkin(1) * t182 + pkin(4) * t203;
t83 = -t184 * t100 + t185 * t102;
t82 = t185 * t101 - t184 * t99;
t81 = t185 * t100 + t184 * t102;
t80 = t184 * t101 + t185 * t99;
t79 = -pkin(4) * t85 + (-pkin(2) * t189 + qJ(3) * t191) * t182;
t78 = pkin(3) * t136 - qJ(3) * t138 - t96;
t77 = pkin(3) * t134 - qJ(3) * t137 - t95;
t76 = pkin(4) * t204 + (pkin(2) * t191 + qJ(3) * t189 + pkin(1)) * t182;
t73 = pkin(3) * t156 - t242 * t137 + t221;
t72 = pkin(3) * t159 - t242 * t138 - t225;
t71 = -pkin(3) * t163 - t75;
t68 = t191 * t105 + t189 * t74;
t67 = t189 * t105 - t191 * t74;
t66 = -pkin(3) * t218 - pkin(4) * t127 - t189 * t71;
t65 = -pkin(3) * t222 + pkin(4) * t128 + t191 * t71;
t62 = pkin(3) * t74 - qJ(3) * t75;
t61 = pkin(3) * t105 - t242 * t75;
t60 = -pkin(4) * t100 - t189 * t72 + t191 * t78;
t59 = -pkin(4) * t99 - t189 * t73 + t191 * t77;
t58 = -pkin(1) * t138 + pkin(4) * t102 + t189 * t78 + t191 * t72;
t57 = -pkin(1) * t137 + pkin(4) * t101 + t189 * t77 + t191 * t73;
t56 = -t184 * t67 + t185 * t68;
t55 = t184 * t68 + t185 * t67;
t54 = -pkin(4) * t67 - t189 * t61 + t191 * t62;
t53 = -pkin(1) * t75 + pkin(4) * t68 + t189 * t62 + t191 * t61;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, 0, 0, 0, 0, 0, 0, -t200, -t202, 0, t70, 0, 0, 0, 0, 0, 0, 0, t200, t202, t64, 0, 0, 0, 0, 0, 0, t82, t83, t90, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, 0, 0, 0, 0, 0, 0, t202, -t200, 0, -t255, 0, 0, 0, 0, 0, 0, 0, -t202, t200, t63, 0, 0, 0, 0, 0, 0, t80, t81, t89, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, 0, 0, 0, 0, 0, 0, t137, t138, 0, t75; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t231, -t226, -t129, -qJ(1) * t129, 0, 0, t202, 0, -t200, 0, -t243, -t258, t255, pkin(4) * t238 + qJ(1) * t255 - t184 * t88, 0, -t202, t200, 0, 0, 0, -t63, t243, t258, -qJ(1) * t63 - t184 * t76 + t185 * t79, -t184 * t112 + t185 * t114, -t184 * t97 + t185 * t98, -t184 * t109 + t185 * t111, -t184 * t113 + t185 * t115, -t184 * t108 + t185 * t110, -t184 * t141 + t185 * t142, -qJ(1) * t80 - t184 * t57 + t185 * t59, -qJ(1) * t81 - t184 * t58 + t185 * t60, -qJ(1) * t89 - t184 * t65 + t185 * t66, -qJ(1) * t55 - t184 * t53 + t185 * t54; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t226, -t231, t130, qJ(1) * t130, 0, 0, t200, 0, t202, 0, t258, -t243, t70, pkin(4) * t239 + qJ(1) * t70 + t185 * t88, 0, -t200, -t202, 0, 0, 0, t64, -t258, t243, qJ(1) * t64 + t184 * t79 + t185 * t76, t185 * t112 + t184 * t114, t184 * t98 + t185 * t97, t185 * t109 + t184 * t111, t185 * t113 + t184 * t115, t185 * t108 + t184 * t110, t185 * t141 + t184 * t142, qJ(1) * t82 + t184 * t59 + t185 * t57, qJ(1) * t83 + t184 * t60 + t185 * t58, qJ(1) * t90 + t184 * t66 + t185 * t65, qJ(1) * t56 + t184 * t54 + t185 * t53; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t165, t166, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t195, -t132 - t241, 0, -pkin(1) * t93, qJDD(2), 0, 0, 0, 0, 0, 0, qJDD(3) - 0.2e1 * t183 - t195, t194 + 0.2e1 * t214 + t241, pkin(1) * t85 - pkin(2) * t116 + qJ(3) * t107, (t158 - t208) * t190, -t190 * t156 - t188 * t159, -t188 * t171 + t219, (-t157 + t207) * t188, t190 * t169 - t224, 0, pkin(1) * t99 + qJ(3) * t156 - t242 * t134 + t225, pkin(1) * t100 + qJ(3) * t159 - t242 * t136 + t221, pkin(1) * t127 - qJ(3) * t163 + t242 * t160 - t74, pkin(1) * t67 + qJ(3) * t105 - t242 * t74;];
tauB_reg = t1;
