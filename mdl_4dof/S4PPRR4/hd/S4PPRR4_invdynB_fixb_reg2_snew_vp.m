% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PPRR4
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
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PPRR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:41
% EndTime: 2019-12-31 16:18:45
% DurationCPUTime: 2.08s
% Computational Cost: add. (3837->245), mult. (6945->378), div. (0->0), fcn. (5272->8), ass. (0->171)
t213 = sin(pkin(6));
t220 = cos(qJ(3));
t222 = qJD(3) ^ 2;
t218 = sin(qJ(3));
t235 = t218 * qJDD(3);
t190 = t220 * t222 + t235;
t234 = t220 * qJDD(3);
t191 = -t218 * t222 + t234;
t212 = sin(pkin(7));
t214 = cos(pkin(7));
t227 = -t212 * t190 + t214 * t191;
t261 = t213 * t227;
t215 = cos(pkin(6));
t260 = t215 * t227;
t195 = t215 * g(1) + t213 * g(2);
t210 = g(3) - qJDD(1);
t174 = -t212 * t195 + t214 * t210;
t175 = -t214 * t195 - t212 * t210;
t123 = t220 * t174 + t218 * t175;
t124 = -t218 * t174 + t220 * t175;
t228 = t218 * t123 + t220 * t124;
t90 = t220 * t123 - t218 * t124;
t253 = t212 * t90;
t75 = t214 * t228 + t253;
t252 = t214 * t90;
t74 = -t212 * t228 + t252;
t194 = t213 * g(1) - t215 * g(2);
t188 = -qJDD(2) + t194;
t155 = pkin(4) * t190 - t220 * t188;
t223 = -pkin(4) * t191 - t218 * t188;
t99 = -qJ(2) * t227 + t212 * t155 + t214 * t223;
t255 = t214 * t190 + t212 * t191;
t100 = qJ(2) * t255 + t214 * t155 - t212 * t223;
t217 = sin(qJ(4));
t208 = t217 ^ 2;
t251 = t208 * t222;
t248 = t213 * t188;
t247 = t213 * t210;
t246 = t214 * t188;
t176 = t215 * t188;
t245 = t215 * t210;
t116 = -qJDD(3) * pkin(3) - t222 * pkin(5) + t123;
t244 = t217 * t116;
t219 = cos(qJ(4));
t202 = t217 * t222 * t219;
t196 = qJDD(4) + t202;
t243 = t217 * t196;
t197 = qJDD(4) - t202;
t242 = t217 * t197;
t241 = t219 * t116;
t240 = t219 * t196;
t239 = t219 * t197;
t117 = -t222 * pkin(3) + qJDD(3) * pkin(5) + t124;
t110 = t219 * t117 - t217 * t188;
t209 = t219 ^ 2;
t238 = t208 + t209;
t237 = qJD(3) * qJD(4);
t236 = t217 * qJDD(3);
t204 = t219 * qJDD(3);
t233 = t217 * t237;
t232 = t219 * t237;
t189 = t238 * qJDD(3);
t206 = t209 * t222;
t192 = t206 + t251;
t148 = t218 * t189 + t220 * t192;
t149 = t220 * t189 - t218 * t192;
t113 = t214 * t148 + t212 * t149;
t114 = -t212 * t148 + t214 * t149;
t109 = t217 * t117 + t219 * t188;
t82 = t217 * t109 + t219 * t110;
t231 = -pkin(1) * t113 - pkin(2) * t148 - pkin(3) * t192 - pkin(5) * t189 + qJ(1) * t114 - t82;
t230 = pkin(1) * t255 + pkin(2) * t190 - qJ(1) * t227 + t124;
t229 = -pkin(1) * t227 - pkin(2) * t191 - qJ(1) * t255 + t123;
t122 = t212 * t174 + t214 * t175;
t151 = -t213 * t194 - t215 * t195;
t225 = t218 * t202;
t224 = t220 * t202;
t81 = t219 * t109 - t217 * t110;
t121 = t214 * t174 - t212 * t175;
t150 = t215 * t194 - t213 * t195;
t221 = qJD(4) ^ 2;
t201 = -t206 - t221;
t200 = t206 - t221;
t199 = -t221 - t251;
t198 = t221 - t251;
t193 = t206 - t251;
t187 = t204 - 0.2e1 * t233;
t186 = t204 - t233;
t185 = t232 + t236;
t184 = 0.2e1 * t232 + t236;
t183 = pkin(1) * t188;
t182 = t238 * t237;
t173 = t218 * qJDD(4) + t220 * t182;
t172 = t219 * t185 - t208 * t237;
t171 = -t220 * qJDD(4) + t218 * t182;
t170 = -t217 * t186 - t209 * t237;
t166 = -t217 * t199 - t239;
t165 = -t217 * t198 + t240;
t164 = t219 * t201 - t243;
t163 = t219 * t200 - t242;
t162 = t219 * t199 - t242;
t161 = -t219 * t198 - t243;
t160 = t217 * t201 + t240;
t159 = -t217 * t200 - t239;
t158 = (-t185 - t232) * t217;
t157 = (-t186 + t233) * t219;
t142 = -t217 * t184 + t219 * t187;
t141 = -t219 * t184 - t217 * t187;
t140 = t215 * t255;
t139 = t213 * t255;
t138 = t220 * t172 - t225;
t137 = t220 * t170 + t225;
t136 = t218 * t172 + t224;
t135 = t218 * t170 - t224;
t134 = t220 * t165 + t217 * t235;
t133 = t220 * t163 + t218 * t204;
t132 = t218 * t165 - t217 * t234;
t131 = t218 * t163 - t219 * t234;
t130 = t220 * t166 + t218 * t184;
t129 = t220 * t164 - t218 * t187;
t128 = t218 * t166 - t220 * t184;
t127 = t218 * t164 + t220 * t187;
t126 = t220 * t142 - t218 * t193;
t125 = t218 * t142 + t220 * t193;
t119 = -t212 * t171 + t214 * t173;
t112 = t215 * t122 - t248;
t111 = t213 * t122 + t176;
t108 = -t212 * t136 + t214 * t138;
t107 = -t212 * t135 + t214 * t137;
t106 = -t212 * t132 + t214 * t134;
t105 = -t212 * t131 + t214 * t133;
t104 = -pkin(5) * t162 + t241;
t103 = -pkin(5) * t160 + t244;
t98 = -t212 * t128 + t214 * t130;
t97 = -t212 * t127 + t214 * t129;
t96 = t214 * t128 + t212 * t130;
t95 = t214 * t127 + t212 * t129;
t94 = -pkin(3) * t162 + t110;
t93 = -pkin(3) * t160 + t109;
t92 = -t212 * t125 + t214 * t126;
t87 = pkin(2) * t188 + pkin(4) * t228;
t86 = t213 * t162 + t215 * t98;
t85 = t213 * t160 + t215 * t97;
t84 = -t215 * t162 + t213 * t98;
t83 = -t215 * t160 + t213 * t97;
t79 = -pkin(4) * t148 + t220 * t81;
t78 = pkin(4) * t149 + t218 * t81;
t77 = t218 * t116 + t220 * t82;
t76 = -t220 * t116 + t218 * t82;
t72 = -pkin(4) * t128 + t220 * t104 - t218 * t94;
t71 = -pkin(4) * t127 + t220 * t103 - t218 * t93;
t70 = t215 * t75 - t248;
t69 = t213 * t75 + t176;
t68 = -pkin(1) * t96 - pkin(2) * t128 + pkin(3) * t184 - pkin(5) * t166 - t244;
t67 = -pkin(1) * t95 - pkin(2) * t127 - pkin(3) * t187 - pkin(5) * t164 + t241;
t66 = -pkin(2) * t162 + pkin(4) * t130 + t218 * t104 + t220 * t94;
t65 = -pkin(2) * t160 + pkin(4) * t129 + t218 * t103 + t220 * t93;
t63 = pkin(1) * t74 + pkin(2) * t90;
t62 = -t212 * t76 + t214 * t77;
t61 = t212 * t77 + t214 * t76;
t60 = -qJ(2) * t113 - t212 * t78 + t214 * t79;
t59 = pkin(4) * t252 + qJ(2) * t74 - t212 * t87;
t58 = -pkin(4) * t76 - (pkin(3) * t218 - pkin(5) * t220) * t81;
t57 = -qJ(2) * t96 - t212 * t66 + t214 * t72;
t56 = -qJ(2) * t95 - t212 * t65 + t214 * t71;
t55 = -t213 * t81 + t215 * t62;
t54 = t213 * t62 + t215 * t81;
t53 = pkin(4) * t77 - (-pkin(3) * t220 - pkin(5) * t218 - pkin(2)) * t81;
t52 = -pkin(1) * t61 - pkin(2) * t76 + pkin(3) * t116 - pkin(5) * t82;
t51 = -qJ(2) * t61 - t212 * t53 + t214 * t58;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, 0, 0, 0, 0, 0, 0, -t140, -t260, 0, t70, 0, 0, 0, 0, 0, 0, t85, t86, t215 * t114, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, 0, 0, 0, 0, 0, 0, -t139, -t261, 0, t69, 0, 0, 0, 0, 0, 0, t83, t84, t213 * t114, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, 0, 0, 0, 0, 0, 0, t227, -t255, 0, -t74, 0, 0, 0, 0, 0, 0, t95, t96, t113, t61; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t247, -t245, -t150, -qJ(1) * t150, 0, 0, 0, 0, 0, 0, -t213 * t174 - t212 * t176, -t213 * t175 - t214 * t176, t215 * t121, -qJ(1) * t111 - (pkin(1) * t213 - qJ(2) * t215) * t121, 0, 0, t260, 0, -t140, t213 * qJDD(3), -t229 * t213 + t215 * t99, t215 * t100 - t230 * t213, t215 * t74, -qJ(1) * t69 - t213 * t63 + t215 * t59, t215 * t108 - t213 * t158, -t213 * t141 + t215 * t92, t215 * t106 - t213 * t161, t215 * t107 - t213 * t157, t215 * t105 - t213 * t159, t215 * t119, -qJ(1) * t83 - t213 * t67 + t215 * t56, -qJ(1) * t84 - t213 * t68 + t215 * t57, -t213 * t231 + t215 * t60, -qJ(1) * t54 - t213 * t52 + t215 * t51; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t245, -t247, t151, qJ(1) * t151, 0, 0, 0, 0, 0, 0, t215 * t174 - t212 * t248, t215 * t175 - t213 * t246, t213 * t121, qJ(1) * t112 - (-pkin(1) * t215 - qJ(2) * t213) * t121, 0, 0, t261, 0, -t139, -t215 * qJDD(3), t213 * t99 + t229 * t215, t213 * t100 + t230 * t215, t213 * t74, qJ(1) * t70 + t213 * t59 + t215 * t63, t213 * t108 + t215 * t158, t215 * t141 + t213 * t92, t213 * t106 + t215 * t161, t213 * t107 + t215 * t157, t213 * t105 + t215 * t159, t213 * t119, qJ(1) * t85 + t213 * t56 + t215 * t67, qJ(1) * t86 + t213 * t57 + t215 * t68, t213 * t60 + t215 * t231, qJ(1) * t55 + t213 * t51 + t215 * t52; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t194, t195, 0, 0, 0, 0, 0, 0, 0, 0, t246, -t212 * t188, t122, qJ(2) * t122 + t183, 0, 0, t255, 0, t227, 0, -t100, t99, t75, pkin(4) * t253 + qJ(2) * t75 + t214 * t87 + t183, t214 * t136 + t212 * t138, t214 * t125 + t212 * t126, t214 * t132 + t212 * t134, t214 * t135 + t212 * t137, t214 * t131 + t212 * t133, t214 * t171 + t212 * t173, -pkin(1) * t160 + qJ(2) * t97 + t212 * t71 + t214 * t65, -pkin(1) * t162 + qJ(2) * t98 + t212 * t72 + t214 * t66, qJ(2) * t114 + t212 * t79 + t214 * t78, pkin(1) * t81 + qJ(2) * t62 + t212 * t58 + t214 * t53;];
tauB_reg = t1;
