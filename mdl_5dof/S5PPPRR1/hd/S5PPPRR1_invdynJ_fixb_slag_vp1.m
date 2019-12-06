% Calculate vector of inverse dynamics joint torques for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPPRR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:50
% EndTime: 2019-12-05 14:58:06
% DurationCPUTime: 8.99s
% Computational Cost: add. (14278->655), mult. (21391->1009), div. (0->0), fcn. (23276->8), ass. (0->242)
t210 = sin(qJ(5));
t211 = cos(qJ(5));
t284 = rSges(6,1) * t211;
t297 = rSges(6,2) * t210 - t284;
t205 = pkin(9) + qJ(4);
t203 = sin(t205);
t207 = sin(pkin(7));
t208 = cos(pkin(8));
t262 = t207 * t208;
t204 = cos(t205);
t209 = cos(pkin(7));
t267 = t204 * t209;
t171 = t203 * t262 + t267;
t206 = sin(pkin(8));
t250 = qJD(4) * t206;
t147 = qJD(5) * t171 + t207 * t250;
t261 = t209 * t203;
t173 = -t207 * t204 + t208 * t261;
t248 = qJD(4) * t209;
t148 = qJD(5) * t173 + t206 * t248;
t247 = qJD(5) * t206;
t249 = qJD(4) * t208;
t185 = t203 * t247 - t249;
t264 = t206 * t210;
t183 = -t204 * t264 - t208 * t211;
t263 = t206 * t211;
t184 = t204 * t263 - t208 * t210;
t269 = t203 * t206;
t241 = t204 * t262;
t172 = t241 - t261;
t141 = -t172 * t210 + t207 * t263;
t142 = t172 * t211 + t207 * t264;
t63 = Icges(6,5) * t142 + Icges(6,6) * t141 + Icges(6,3) * t171;
t272 = Icges(6,4) * t142;
t65 = Icges(6,2) * t141 + Icges(6,6) * t171 + t272;
t135 = Icges(6,4) * t141;
t67 = Icges(6,1) * t142 + Icges(6,5) * t171 + t135;
t28 = t183 * t65 + t184 * t67 + t269 * t63;
t174 = t203 * t207 + t208 * t267;
t143 = -t174 * t210 + t209 * t263;
t144 = t174 * t211 + t209 * t264;
t64 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t173;
t271 = Icges(6,4) * t144;
t66 = Icges(6,2) * t143 + Icges(6,6) * t173 + t271;
t136 = Icges(6,4) * t143;
t68 = Icges(6,1) * t144 + Icges(6,5) * t173 + t136;
t29 = t183 * t66 + t184 * t68 + t269 * t64;
t115 = Icges(6,5) * t184 + Icges(6,6) * t183 + Icges(6,3) * t269;
t270 = Icges(6,4) * t184;
t116 = Icges(6,2) * t183 + Icges(6,6) * t269 + t270;
t179 = Icges(6,4) * t183;
t117 = Icges(6,1) * t184 + Icges(6,5) * t269 + t179;
t43 = t115 * t269 + t116 * t183 + t117 * t184;
t296 = (t147 * t28 + t148 * t29 + t185 * t43) * t204;
t176 = (-Icges(5,5) * t203 - Icges(5,6) * t204) * t206;
t160 = qJD(4) * t176;
t157 = qJD(4) * t241 - t203 * t248;
t244 = qJDD(4) * t206;
t95 = qJD(5) * t157 + qJDD(5) * t171 + t207 * t244;
t295 = t95 / 0.2e1;
t159 = t174 * qJD(4);
t96 = qJD(5) * t159 + qJDD(5) * t173 + t209 * t244;
t294 = t96 / 0.2e1;
t293 = -t147 / 0.2e1;
t292 = t147 / 0.2e1;
t291 = -t148 / 0.2e1;
t290 = t148 / 0.2e1;
t243 = qJDD(4) * t208;
t251 = qJD(4) * t204;
t149 = -t243 + (qJD(5) * t251 + qJDD(5) * t203) * t206;
t289 = t149 / 0.2e1;
t288 = -t185 / 0.2e1;
t287 = t185 / 0.2e1;
t156 = t171 * qJD(4);
t113 = -pkin(4) * t156 + pkin(6) * t157;
t91 = -qJD(5) * t142 + t156 * t210;
t92 = qJD(5) * t141 - t156 * t211;
t55 = rSges(6,1) * t92 + rSges(6,2) * t91 + rSges(6,3) * t157;
t282 = t113 + t55;
t158 = t173 * qJD(4);
t114 = -pkin(4) * t158 + pkin(6) * t159;
t93 = -qJD(5) * t144 + t158 * t210;
t94 = qJD(5) * t143 - t158 * t211;
t56 = rSges(6,1) * t94 + rSges(6,2) * t93 + rSges(6,3) * t159;
t281 = -t114 - t56;
t276 = Icges(5,4) * t172;
t266 = t206 * t207;
t99 = -Icges(5,2) * t171 + Icges(5,6) * t266 + t276;
t280 = -Icges(5,1) * t171 - t276 - t99;
t128 = pkin(4) * t172 + pkin(6) * t171;
t69 = rSges(6,1) * t142 + rSges(6,2) * t141 + rSges(6,3) * t171;
t279 = t128 + t69;
t130 = pkin(4) * t174 + pkin(6) * t173;
t70 = rSges(6,1) * t144 + rSges(6,2) * t143 + rSges(6,3) * t173;
t278 = -t130 - t70;
t168 = (-pkin(4) * t203 + pkin(6) * t204) * t250;
t240 = t203 * t250;
t145 = -qJD(5) * t184 + t210 * t240;
t146 = qJD(5) * t183 - t211 * t240;
t239 = t204 * t250;
t74 = rSges(6,1) * t146 + rSges(6,2) * t145 + rSges(6,3) * t239;
t277 = t168 + t74;
t275 = Icges(5,4) * t174;
t274 = Icges(5,4) * t203;
t273 = Icges(5,4) * t204;
t268 = t204 * t206;
t265 = t206 * t209;
t100 = -Icges(5,2) * t173 + Icges(5,6) * t265 + t275;
t260 = -Icges(5,1) * t173 - t100 - t275;
t164 = Icges(5,4) * t171;
t101 = Icges(5,1) * t172 + Icges(5,5) * t266 - t164;
t259 = Icges(5,2) * t172 - t101 + t164;
t165 = Icges(5,4) * t173;
t102 = Icges(5,1) * t174 + Icges(5,5) * t265 - t165;
t258 = Icges(5,2) * t174 - t102 + t165;
t118 = rSges(6,1) * t184 + rSges(6,2) * t183 + rSges(6,3) * t269;
t182 = (pkin(4) * t204 + pkin(6) * t203) * t206;
t257 = t118 + t182;
t151 = -Icges(5,6) * t208 + (-Icges(5,2) * t203 + t273) * t206;
t178 = (-Icges(5,1) * t203 - t273) * t206;
t256 = t151 - t178;
t152 = -Icges(5,5) * t208 + (Icges(5,1) * t204 - t274) * t206;
t177 = (-Icges(5,2) * t204 - t274) * t206;
t255 = t152 + t177;
t254 = t203 * rSges(6,2) * t264 + rSges(6,3) * t268;
t252 = qJD(3) * t206;
t253 = qJD(2) * t207 + t209 * t252;
t245 = qJDD(3) * t206;
t186 = qJDD(2) * t207 + t209 * t245;
t246 = -m(4) - m(5) - m(6);
t242 = -m(3) + t246;
t236 = t250 / 0.2e1;
t127 = -t171 * pkin(4) + pkin(6) * t172;
t129 = -t173 * pkin(4) + pkin(6) * t174;
t233 = -qJD(2) * t209 + t207 * t252;
t232 = -qJD(3) * t208 + qJD(1);
t199 = -qJDD(3) * t208 + qJDD(1);
t105 = -Icges(5,5) * t156 - Icges(5,6) * t157;
t106 = -Icges(5,5) * t158 - Icges(5,6) * t159;
t107 = -Icges(5,4) * t156 - Icges(5,2) * t157;
t108 = -Icges(5,4) * t158 - Icges(5,2) * t159;
t109 = -Icges(5,1) * t156 - Icges(5,4) * t157;
t110 = -Icges(5,1) * t158 - Icges(5,4) * t159;
t231 = t207 * (-t101 * t156 + t105 * t266 - t107 * t171 + t109 * t172 - t157 * t99) + t209 * (-t100 * t157 - t102 * t156 + t106 * t266 - t108 * t171 + t110 * t172);
t230 = t207 * (-t101 * t158 + t105 * t265 - t107 * t173 + t109 * t174 - t159 * t99) + t209 * (-t100 * t159 - t102 * t158 + t106 * t265 - t108 * t173 + t110 * t174);
t229 = t207 * (-t105 * t208 + (-t107 * t203 + t109 * t204 + (-t101 * t203 - t204 * t99) * qJD(4)) * t206) + t209 * (-t106 * t208 + (-t108 * t203 + t110 * t204 + (-t100 * t204 - t102 * t203) * qJD(4)) * t206);
t97 = Icges(5,5) * t172 - Icges(5,6) * t171 + Icges(5,3) * t266;
t98 = Icges(5,5) * t174 - Icges(5,6) * t173 + Icges(5,3) * t265;
t228 = t207 * (t101 * t172 - t171 * t99 + t266 * t97) + t209 * (-t100 * t171 + t102 * t172 + t266 * t98);
t227 = t207 * (t101 * t174 - t173 * t99 + t265 * t97) + t209 * (-t100 * t173 + t102 * t174 + t265 * t98);
t111 = -rSges(5,1) * t156 - rSges(5,2) * t157;
t180 = (-rSges(5,1) * t203 - rSges(5,2) * t204) * t206;
t163 = qJD(4) * t180;
t103 = rSges(5,1) * t172 - rSges(5,2) * t171 + rSges(5,3) * t266;
t155 = -rSges(5,3) * t208 + (rSges(5,1) * t204 - rSges(5,2) * t203) * t206;
t217 = t103 * t208 + t155 * t266;
t44 = t217 * qJDD(4) + (t111 * t208 + t163 * t266) * qJD(4) + t186;
t104 = rSges(5,1) * t174 - rSges(5,2) * t173 + rSges(5,3) * t265;
t112 = -rSges(5,1) * t158 - rSges(5,2) * t159;
t195 = t207 * t245;
t45 = t195 + (-qJD(4) * t112 - qJDD(4) * t104) * t208 + (-qJDD(2) + (-qJD(4) * t163 - qJDD(4) * t155) * t206) * t209;
t226 = t207 * t44 - t209 * t45;
t225 = t207 * (-t208 * t97 + (t101 * t204 - t203 * t99) * t206) + t209 * (-t208 * t98 + (-t100 * t203 + t102 * t204) * t206);
t60 = qJD(4) * t217 + t253;
t61 = (-t104 * t208 - t155 * t265) * qJD(4) + t233;
t224 = t207 * t60 - t209 * t61;
t223 = -Icges(6,1) * t211 + Icges(6,4) * t210;
t222 = -Icges(6,4) * t211 + Icges(6,2) * t210;
t221 = -Icges(6,5) * t211 + Icges(6,6) * t210;
t220 = t103 * t209 - t104 * t207;
t219 = t111 * t209 - t112 * t207;
t218 = t128 * t209 - t130 * t207;
t81 = t172 * rSges(6,3) + t297 * t171;
t82 = t174 * rSges(6,3) + t297 * t173;
t216 = t128 * t208 + t182 * t266;
t215 = (Icges(6,5) * t183 - Icges(6,6) * t184) * t185 + t147 * (Icges(6,5) * t141 - Icges(6,6) * t142) + t148 * (Icges(6,5) * t143 - Icges(6,6) * t144);
t214 = (Icges(6,1) * t143 - t271 - t66) * t148 + (Icges(6,1) * t141 - t272 - t65) * t147 + (Icges(6,1) * t183 - t116 - t270) * t185;
t213 = (-Icges(6,2) * t144 + t136 + t68) * t148 + (-Icges(6,2) * t142 + t135 + t67) * t147 + (-Icges(6,2) * t184 + t117 + t179) * t185;
t212 = (Icges(6,3) * t174 + t173 * t221 + t210 * t66 - t211 * t68) * t148 + (Icges(6,3) * t172 + t171 * t221 + t210 * t65 - t211 * t67) * t147 + (t116 * t210 - t117 * t211 + (Icges(6,3) * t204 + t203 * t221) * t206) * t185;
t192 = pkin(6) * t268;
t187 = -qJDD(2) * t209 + t195;
t162 = qJD(4) * t178;
t161 = qJD(4) * t177;
t150 = -Icges(5,3) * t208 + (Icges(5,5) * t204 - Icges(5,6) * t203) * t206;
t140 = -rSges(6,1) * t203 * t263 + t254;
t139 = (Icges(6,5) * t204 + t203 * t223) * t206;
t138 = (Icges(6,6) * t204 + t203 * t222) * t206;
t134 = rSges(6,1) * t183 - rSges(6,2) * t184;
t126 = -rSges(5,1) * t173 - rSges(5,2) * t174;
t125 = -rSges(5,1) * t171 - rSges(5,2) * t172;
t120 = -Icges(5,5) * t173 - Icges(5,6) * t174;
t119 = -Icges(5,5) * t171 - Icges(5,6) * t172;
t90 = rSges(6,1) * t143 - rSges(6,2) * t144;
t89 = rSges(6,1) * t141 - rSges(6,2) * t142;
t80 = Icges(6,5) * t174 + t173 * t223;
t79 = Icges(6,5) * t172 + t171 * t223;
t78 = Icges(6,6) * t174 + t173 * t222;
t77 = Icges(6,6) * t172 + t171 * t222;
t73 = Icges(6,1) * t146 + Icges(6,4) * t145 + Icges(6,5) * t239;
t72 = Icges(6,4) * t146 + Icges(6,2) * t145 + Icges(6,6) * t239;
t71 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t239;
t62 = -t150 * t208 + (-t151 * t203 + t152 * t204) * t206;
t59 = t150 * t265 - t151 * t173 + t152 * t174;
t58 = t150 * t266 - t151 * t171 + t152 * t172;
t57 = t220 * t250 + t232;
t54 = Icges(6,1) * t94 + Icges(6,4) * t93 + Icges(6,5) * t159;
t53 = Icges(6,1) * t92 + Icges(6,4) * t91 + Icges(6,5) * t157;
t52 = Icges(6,4) * t94 + Icges(6,2) * t93 + Icges(6,6) * t159;
t51 = Icges(6,4) * t92 + Icges(6,2) * t91 + Icges(6,6) * t157;
t50 = Icges(6,5) * t94 + Icges(6,6) * t93 + Icges(6,3) * t159;
t49 = Icges(6,5) * t92 + Icges(6,6) * t91 + Icges(6,3) * t157;
t48 = -t160 * t208 + (-t161 * t203 + t162 * t204 + (-t151 * t204 - t152 * t203) * qJD(4)) * t206;
t38 = -t151 * t159 - t152 * t158 + t160 * t265 - t161 * t173 + t162 * t174;
t37 = -t151 * t157 - t152 * t156 + t160 * t266 - t161 * t171 + t162 * t172;
t36 = t115 * t173 + t116 * t143 + t117 * t144;
t35 = t115 * t171 + t116 * t141 + t117 * t142;
t34 = (qJD(4) * t219 + qJDD(4) * t220) * t206 + t199;
t33 = -t118 * t148 + t185 * t70 + (-t130 * t208 - t182 * t265) * qJD(4) + t233;
t32 = qJD(4) * t216 + t118 * t147 - t185 * t69 + t253;
t23 = -t147 * t70 + t148 * t69 + t218 * t250 + t232;
t22 = t143 * t66 + t144 * t68 + t173 * t64;
t21 = t143 * t65 + t144 * t67 + t173 * t63;
t20 = t141 * t66 + t142 * t68 + t171 * t64;
t19 = t141 * t65 + t142 * t67 + t171 * t63;
t18 = t116 * t145 + t117 * t146 + t183 * t72 + t184 * t73 + (t115 * t251 + t203 * t71) * t206;
t17 = t115 * t159 + t116 * t93 + t117 * t94 + t143 * t72 + t144 * t73 + t173 * t71;
t16 = t115 * t157 + t116 * t91 + t117 * t92 + t141 * t72 + t142 * t73 + t171 * t71;
t15 = -t118 * t96 - t148 * t74 + t149 * t70 + t185 * t56 + t195 + (-qJD(4) * t114 - qJDD(4) * t130) * t208 + (-qJDD(2) + (-qJD(4) * t168 - qJDD(4) * t182) * t206) * t209;
t14 = t118 * t95 + t147 * t74 - t149 * t69 - t185 * t55 + t216 * qJDD(4) + (t113 * t208 + t168 * t266) * qJD(4) + t186;
t13 = t145 * t66 + t146 * t68 + t183 * t52 + t184 * t54 + (t203 * t50 + t251 * t64) * t206;
t12 = t145 * t65 + t146 * t67 + t183 * t51 + t184 * t53 + (t203 * t49 + t251 * t63) * t206;
t11 = -t147 * t56 + t148 * t55 + t69 * t96 - t70 * t95 + (t218 * qJDD(4) + (t113 * t209 - t114 * t207) * qJD(4)) * t206 + t199;
t9 = t143 * t52 + t144 * t54 + t159 * t64 + t173 * t50 + t66 * t93 + t68 * t94;
t8 = t143 * t51 + t144 * t53 + t159 * t63 + t173 * t49 + t65 * t93 + t67 * t94;
t7 = t141 * t52 + t142 * t54 + t157 * t64 + t171 * t50 + t66 * t91 + t68 * t92;
t6 = t141 * t51 + t142 * t53 + t157 * t63 + t171 * t49 + t65 * t91 + t67 * t92;
t5 = t147 * t21 + t148 * t22 + t185 * t36;
t4 = t147 * t19 + t148 * t20 + t185 * t35;
t3 = t12 * t147 + t13 * t148 + t149 * t43 + t185 * t18 + t28 * t95 + t29 * t96;
t2 = t147 * t8 + t148 * t9 + t149 * t36 + t17 * t185 + t21 * t95 + t96 * t22;
t1 = t147 * t6 + t148 * t7 + t149 * t35 + t16 * t185 + t95 * t19 + t20 * t96;
t10 = [(m(2) + m(3)) * qJDD(1) + m(4) * t199 + m(5) * t34 + m(6) * t11 + (-m(2) + t242) * g(3); t242 * (g(1) * t207 - g(2) * t209) + m(4) * (t186 * t207 - t187 * t209) + m(5) * t226 + m(6) * (t14 * t207 - t15 * t209) + m(3) * (t207 ^ 2 + t209 ^ 2) * qJDD(2); t246 * (-g(3) * t208 + (g(1) * t209 + g(2) * t207) * t206) + m(4) * (-t199 * t208 + (t186 * t209 + t187 * t207) * t206) + m(5) * (-t208 * t34 + (t207 * t45 + t209 * t44) * t206) + m(6) * (-t11 * t208 + (t14 * t209 + t15 * t207) * t206); (t208 ^ 2 * t160 + (((t207 * t280 + t209 * t260) * t204 + (t207 * t259 + t209 * t258) * t203) * t206 + (-t119 * t207 - t120 * t209 + t203 * t255 + t204 * t256) * t208) * t250) * t249 / 0.2e1 - (t206 * t229 - t208 * t48) * t249 / 0.2e1 - (t206 * t225 - t208 * t62) * t243 / 0.2e1 + (-t18 * t208 + (t12 * t207 + t13 * t209) * t206) * t287 + ((t183 * t78 + t184 * t80) * t148 + (t183 * t77 + t184 * t79) * t147 + (t183 * t138 + t184 * t139) * t185 + (t172 * t28 + t174 * t29) * qJD(5) + ((t43 * qJD(5) + t115 * t185 + t63 * t147 + t64 * t148) * t204 + t212 * t203) * t206) * t288 + (-t208 * t43 + (t207 * t28 + t209 * t29) * t206) * t289 + (-t17 * t208 + (t207 * t8 + t209 * t9) * t206) * t290 + ((t143 * t78 + t144 * t80 + t174 * t64) * t148 + (t143 * t77 + t144 * t79 + t174 * t63) * t147 + (t174 * t115 + t143 * t138 + t144 * t139) * t185 + (t172 * t21 + t174 * t22 + t268 * t36) * qJD(5) + t212 * t173) * t291 + (-t16 * t208 + (t207 * t6 + t209 * t7) * t206) * t292 + ((t141 * t78 + t142 * t80 + t172 * t64) * t148 + (t141 * t77 + t142 * t79 + t172 * t63) * t147 + (t172 * t115 + t141 * t138 + t142 * t139) * t185 + (t172 * t19 + t174 * t20 + t268 * t35) * qJD(5) + t212 * t171) * t293 + (-t208 * t36 + (t207 * t21 + t209 * t22) * t206) * t294 + (-t208 * t35 + (t19 * t207 + t20 * t209) * t206) * t295 - t247 * t296 / 0.2e1 - ((-qJD(4) * t48 - qJDD(4) * t62) * t208 + (qJD(4) * t229 + qJDD(4) * t225) * t206 + t3) * t208 / 0.2e1 - (t172 * t4 + t174 * t5) * qJD(5) / 0.2e1 + ((-qJD(4) * t37 - qJDD(4) * t58) * t208 + (qJD(4) * t231 + qJDD(4) * t228) * t206 + t1) * t266 / 0.2e1 + ((-qJD(4) * t38 - qJDD(4) * t59) * t208 + (qJD(4) * t230 + qJDD(4) * t227) * t206 + t2) * t265 / 0.2e1 + (t209 * (t230 * t206 - t208 * t38) + t207 * (t231 * t206 - t208 * t37)) * t236 - (t209 * ((t120 * t265 + t173 * t258 + t174 * t260) * t265 + (t119 * t265 + t173 * t259 + t174 * t280) * t266 - (-t173 * t255 - t174 * t256 + t176 * t265) * t208) + t207 * ((t120 * t266 + t171 * t258 + t172 * t260) * t265 + (t119 * t266 + t171 * t259 + t172 * t280) * t266 - (-t171 * t255 - t172 * t256 + t176 * t266) * t208)) * t206 * qJD(4) ^ 2 / 0.2e1 + (t209 * (t227 * t206 - t208 * t59) + t207 * (t228 * t206 - t208 * t58)) * t244 / 0.2e1 + ((t14 * t279 + t15 * t278 + t281 * t33 + t282 * t32) * t208 + ((t11 * t279 - t15 * t257 + t23 * t282 - t277 * t33) * t209 + (t11 * t278 + t14 * t257 + t23 * t281 + t277 * t32) * t207) * t206 - t32 * (t140 * t147 - t185 * t81) - t33 * (-t140 * t148 + t185 * t82) - t23 * (-t147 * t82 + t148 * t81) - (t32 * (t118 * t172 - t268 * t69) + t33 * (-t118 * t174 + t268 * t70) + t23 * (-t172 * t70 + t174 * t69)) * qJD(5) - ((t32 * t127 - t33 * t129) * t208 + (t23 * (t127 * t209 - t129 * t207) + (t207 * t32 - t209 * t33) * (-pkin(4) * t269 + t192)) * t206) * qJD(4) - g(1) * (t129 + t82) - g(2) * (t127 + t81) - g(3) * (t192 + (-pkin(4) - t284) * t269 + t254)) * m(6) + (-g(1) * t126 - g(2) * t125 - g(3) * t180 - ((t60 * t125 - t61 * t126) * t208 + (t57 * (t125 * t209 - t126 * t207) + t224 * t180) * t206) * qJD(4) + (t103 * t44 - t104 * t45 + t111 * t60 - t112 * t61) * t208 + (t155 * t226 + t163 * t224 + t219 * t57 + t220 * t34) * t206) * m(5); t159 * t5 / 0.2e1 + t173 * t2 / 0.2e1 + (t171 * t21 + t173 * t22 + t269 * t36) * t294 + (t157 * t21 + t159 * t22 + t171 * t8 + t173 * t9 + (t17 * t203 + t251 * t36) * t206) * t290 + t157 * t4 / 0.2e1 + t171 * t1 / 0.2e1 + (t171 * t19 + t173 * t20 + t269 * t35) * t295 + (t157 * t19 + t159 * t20 + t171 * t6 + t173 * t7 + (t16 * t203 + t251 * t35) * t206) * t292 + t236 * t296 + t3 * t269 / 0.2e1 + (t171 * t28 + t173 * t29 + t269 * t43) * t289 + (t12 * t171 + t13 * t173 + t157 * t28 + t159 * t29 + (t18 * t203 + t251 * t43) * t206) * t287 + (t143 * t213 + t144 * t214 + t173 * t215) * t291 + (t141 * t213 + t142 * t214 + t171 * t215) * t293 + (t183 * t213 + t184 * t214 + t215 * t269) * t288 + (t11 * (-t171 * t70 + t173 * t69) + (t171 * t32 - t173 * t33) * t74 + (t14 * t171 - t15 * t173 + t157 * t32 - t159 * t33) * t118 + ((-t32 * t69 + t33 * t70) * t251 + (-t14 * t69 + t15 * t70 - t32 * t55 + t33 * t56) * t203) * t206 - t32 * (t134 * t147 - t185 * t89) - t33 * (-t134 * t148 + t185 * t90) - g(1) * t90 - g(2) * t89 - g(3) * t134 + (t147 * t90 - t148 * t89 - t157 * t70 + t159 * t69 - t171 * t56 + t173 * t55) * t23) * m(6);];
tau = t10;
