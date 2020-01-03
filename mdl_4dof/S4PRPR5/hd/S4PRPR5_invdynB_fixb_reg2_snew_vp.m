% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PRPR5
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PRPR5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR5_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:19
% EndTime: 2019-12-31 16:23:23
% DurationCPUTime: 2.29s
% Computational Cost: add. (4359->263), mult. (7789->395), div. (0->0), fcn. (5338->8), ass. (0->181)
t234 = sin(pkin(6));
t233 = sin(pkin(7));
t235 = cos(pkin(7));
t243 = qJD(2) ^ 2;
t205 = t233 * qJDD(2) + t235 * t243;
t206 = t235 * qJDD(2) - t233 * t243;
t239 = sin(qJ(2));
t241 = cos(qJ(2));
t249 = -t239 * t205 + t241 * t206;
t286 = t234 * t249;
t236 = cos(pkin(6));
t285 = t236 * t249;
t213 = t236 * g(1) + t234 * g(2);
t231 = g(3) - qJDD(1);
t189 = -t241 * t213 - t239 * t231;
t183 = -t243 * pkin(2) + t189;
t188 = -t239 * t213 + t241 * t231;
t244 = qJDD(2) * pkin(2) - t188;
t128 = t233 * t183 - t235 * t244;
t129 = t235 * t183 + t233 * t244;
t250 = t233 * t128 + t235 * t129;
t99 = t235 * t128 - t233 * t129;
t278 = t239 * t99;
t84 = t241 * t250 + t278;
t277 = t241 * t99;
t83 = -t239 * t250 + t277;
t212 = t234 * g(1) - t236 * g(2);
t204 = -qJDD(3) + t212;
t163 = qJ(3) * t205 - t235 * t204;
t245 = -qJ(3) * t206 - t233 * t204;
t110 = -pkin(4) * t249 + t239 * t163 + t241 * t245;
t280 = t241 * t205 + t239 * t206;
t111 = pkin(4) * t280 + t241 * t163 - t239 * t245;
t196 = t236 * t212;
t165 = -t234 * t213 + t196;
t238 = sin(qJ(4));
t229 = t238 ^ 2;
t276 = t229 * t243;
t209 = t241 * qJDD(2) - t239 * t243;
t275 = t234 * t209;
t274 = t234 * t212;
t272 = t234 * t231;
t271 = t236 * t209;
t270 = t236 * t231;
t125 = -qJDD(2) * pkin(3) - t243 * pkin(5) + t128;
t269 = t238 * t125;
t240 = cos(qJ(4));
t220 = t238 * t243 * t240;
t214 = qJDD(4) + t220;
t268 = t238 * t214;
t215 = qJDD(4) - t220;
t267 = t238 * t215;
t264 = t240 * t125;
t263 = t240 * t214;
t262 = t240 * t215;
t126 = -t243 * pkin(3) + qJDD(2) * pkin(5) + t129;
t119 = t240 * t126 - t238 * t204;
t230 = t240 ^ 2;
t261 = t229 + t230;
t260 = qJD(2) * qJD(4);
t259 = t236 * qJDD(2);
t258 = t238 * qJDD(2);
t225 = t240 * qJDD(2);
t257 = t238 * t260;
t256 = t240 * t260;
t207 = t261 * qJDD(2);
t227 = t230 * t243;
t210 = t227 + t276;
t161 = t233 * t207 + t235 * t210;
t164 = t235 * t207 - t233 * t210;
t122 = t241 * t161 + t239 * t164;
t123 = -t239 * t161 + t241 * t164;
t118 = t238 * t126 + t240 * t204;
t91 = t238 * t118 + t240 * t119;
t255 = -pkin(1) * t122 - pkin(2) * t161 - pkin(3) * t210 - pkin(5) * t207 + qJ(1) * t123 - t91;
t254 = pkin(1) * t280 + pkin(2) * t205 - qJ(1) * t249 + t129;
t253 = -pkin(1) * t249 - pkin(2) * t206 - qJ(1) * t280 + t128;
t208 = t239 * qJDD(2) + t241 * t243;
t252 = -pkin(1) * t208 + qJ(1) * t209 - t189;
t251 = pkin(1) * t209 + qJ(1) * t208 - t188;
t135 = t239 * t188 + t241 * t189;
t166 = -t236 * t213 - t274;
t247 = t233 * t220;
t246 = t235 * t220;
t171 = pkin(4) * t208 - t241 * t212;
t170 = -pkin(4) * t209 - t239 * t212;
t90 = t240 * t118 - t238 * t119;
t134 = t241 * t188 - t239 * t189;
t242 = qJD(4) ^ 2;
t223 = t234 * qJDD(2);
t219 = -t227 - t242;
t218 = t227 - t242;
t217 = -t242 - t276;
t216 = t242 - t276;
t211 = t227 - t276;
t203 = t225 - 0.2e1 * t257;
t202 = t225 - t257;
t201 = t256 + t258;
t200 = 0.2e1 * t256 + t258;
t199 = t261 * t260;
t195 = t236 * t208;
t194 = t234 * t208;
t187 = t240 * t201 - t229 * t260;
t186 = -t238 * t202 - t230 * t260;
t185 = t233 * qJDD(4) + t235 * t199;
t184 = -t235 * qJDD(4) + t233 * t199;
t181 = -t238 * t217 - t262;
t180 = -t238 * t216 + t263;
t179 = t240 * t219 - t268;
t178 = t240 * t218 - t267;
t177 = t240 * t217 - t267;
t176 = -t240 * t216 - t268;
t175 = t238 * t219 + t263;
t174 = -t238 * t218 - t262;
t173 = (-t201 - t256) * t238;
t172 = (-t202 + t257) * t240;
t153 = -t238 * t200 + t240 * t203;
t152 = -t240 * t200 - t238 * t203;
t149 = t236 * t280;
t148 = t234 * t280;
t147 = t235 * t187 - t247;
t146 = t235 * t186 + t247;
t145 = t233 * t187 + t246;
t144 = t233 * t186 - t246;
t143 = t235 * t180 + t233 * t258;
t142 = t235 * t178 + t225 * t233;
t141 = t233 * t180 - t235 * t258;
t140 = t233 * t178 - t225 * t235;
t139 = t235 * t181 + t233 * t200;
t138 = t235 * t179 - t233 * t203;
t137 = t233 * t181 - t235 * t200;
t136 = t233 * t179 + t235 * t203;
t132 = t235 * t153 - t233 * t211;
t131 = t233 * t153 + t235 * t211;
t130 = -t239 * t184 + t241 * t185;
t121 = t236 * t135 - t274;
t120 = t234 * t135 + t196;
t117 = -t239 * t145 + t241 * t147;
t116 = -t239 * t144 + t241 * t146;
t115 = -t239 * t141 + t241 * t143;
t114 = -t239 * t140 + t241 * t142;
t113 = -pkin(5) * t177 + t264;
t112 = -pkin(5) * t175 + t269;
t109 = -t239 * t137 + t241 * t139;
t108 = -t239 * t136 + t241 * t138;
t107 = t241 * t137 + t239 * t139;
t106 = t241 * t136 + t239 * t138;
t103 = -t239 * t131 + t241 * t132;
t102 = -pkin(3) * t177 + t119;
t101 = -pkin(3) * t175 + t118;
t96 = t236 * t109 + t234 * t177;
t95 = t236 * t108 + t234 * t175;
t94 = t234 * t109 - t236 * t177;
t93 = t234 * t108 - t236 * t175;
t92 = pkin(2) * t204 + qJ(3) * t250;
t88 = -qJ(3) * t161 + t235 * t90;
t87 = qJ(3) * t164 + t233 * t90;
t86 = t233 * t125 + t235 * t91;
t85 = -t235 * t125 + t233 * t91;
t81 = -pkin(1) * t107 - pkin(2) * t137 + pkin(3) * t200 - pkin(5) * t181 - t269;
t80 = -pkin(1) * t106 - pkin(2) * t136 - pkin(3) * t203 - pkin(5) * t179 + t264;
t79 = -qJ(3) * t137 - t233 * t102 + t235 * t113;
t78 = -qJ(3) * t136 - t233 * t101 + t235 * t112;
t77 = -t234 * t204 + t236 * t84;
t76 = t236 * t204 + t234 * t84;
t75 = -pkin(2) * t177 + qJ(3) * t139 + t235 * t102 + t233 * t113;
t74 = -pkin(2) * t175 + qJ(3) * t138 + t235 * t101 + t233 * t112;
t72 = pkin(1) * t83 + pkin(2) * t99;
t71 = -pkin(4) * t122 - t239 * t87 + t241 * t88;
t70 = -t239 * t85 + t241 * t86;
t69 = t239 * t86 + t241 * t85;
t68 = pkin(4) * t83 + qJ(3) * t277 - t239 * t92;
t67 = -qJ(3) * t85 - (pkin(3) * t233 - pkin(5) * t235) * t90;
t66 = -pkin(4) * t107 - t239 * t75 + t241 * t79;
t65 = -pkin(4) * t106 - t239 * t74 + t241 * t78;
t64 = -t234 * t90 + t236 * t70;
t63 = t234 * t70 + t236 * t90;
t62 = qJ(3) * t86 - (-pkin(3) * t235 - pkin(5) * t233 - pkin(2)) * t90;
t61 = -pkin(1) * t69 - pkin(2) * t85 + pkin(3) * t125 - pkin(5) * t91;
t60 = -pkin(4) * t69 - t239 * t62 + t241 * t67;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, 0, 0, 0, 0, 0, 0, -t195, -t271, 0, t121, 0, 0, 0, 0, 0, 0, -t149, -t285, 0, t77, 0, 0, 0, 0, 0, 0, t95, t96, t236 * t123, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, 0, 0, 0, 0, 0, 0, -t194, -t275, 0, t120, 0, 0, 0, 0, 0, 0, -t148, -t286, 0, t76, 0, 0, 0, 0, 0, 0, t93, t94, t234 * t123, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t231, 0, 0, 0, 0, 0, 0, t209, -t208, 0, -t134, 0, 0, 0, 0, 0, 0, t249, -t280, 0, -t83, 0, 0, 0, 0, 0, 0, t106, t107, t122, t69; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t272, -t270, -t165, -qJ(1) * t165, 0, 0, t271, 0, -t195, t223, t236 * t170 + t234 * t251, t236 * t171 + t234 * t252, t236 * t134, -qJ(1) * t120 - (pkin(1) * t234 - pkin(4) * t236) * t134, 0, 0, t285, 0, -t149, t223, t236 * t110 - t234 * t253, t236 * t111 - t234 * t254, t236 * t83, -qJ(1) * t76 - t234 * t72 + t236 * t68, t236 * t117 - t234 * t173, t236 * t103 - t234 * t152, t236 * t115 - t234 * t176, t236 * t116 - t234 * t172, t236 * t114 - t234 * t174, t236 * t130, -qJ(1) * t93 - t234 * t80 + t236 * t65, -qJ(1) * t94 - t234 * t81 + t236 * t66, -t234 * t255 + t236 * t71, -qJ(1) * t63 - t234 * t61 + t236 * t60; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t270, -t272, t166, qJ(1) * t166, 0, 0, t275, 0, -t194, -t259, t234 * t170 - t236 * t251, t234 * t171 - t236 * t252, t234 * t134, qJ(1) * t121 - (-pkin(1) * t236 - pkin(4) * t234) * t134, 0, 0, t286, 0, -t148, -t259, t234 * t110 + t236 * t253, t234 * t111 + t236 * t254, t234 * t83, qJ(1) * t77 + t234 * t68 + t236 * t72, t234 * t117 + t236 * t173, t234 * t103 + t236 * t152, t234 * t115 + t236 * t176, t234 * t116 + t236 * t172, t234 * t114 + t236 * t174, t234 * t130, qJ(1) * t95 + t234 * t65 + t236 * t80, qJ(1) * t96 + t234 * t66 + t236 * t81, t234 * t71 + t236 * t255, qJ(1) * t64 + t234 * t60 + t236 * t61; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t212, t213, 0, 0, 0, 0, t208, 0, t209, 0, -t171, t170, t135, pkin(1) * t212 + pkin(4) * t135, 0, 0, t280, 0, t249, 0, -t111, t110, t84, pkin(1) * t204 + pkin(4) * t84 + qJ(3) * t278 + t241 * t92, t241 * t145 + t239 * t147, t241 * t131 + t239 * t132, t241 * t141 + t239 * t143, t241 * t144 + t239 * t146, t241 * t140 + t239 * t142, t241 * t184 + t239 * t185, -pkin(1) * t175 + pkin(4) * t108 + t239 * t78 + t241 * t74, -pkin(1) * t177 + pkin(4) * t109 + t239 * t79 + t241 * t75, pkin(4) * t123 + t239 * t88 + t241 * t87, pkin(1) * t90 + pkin(4) * t70 + t239 * t67 + t241 * t62;];
tauB_reg = t1;
