% Calculate vector of inverse dynamics joint torques for
% S5PRRPR1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:34
% EndTime: 2019-12-05 16:15:41
% DurationCPUTime: 3.12s
% Computational Cost: add. (9248->359), mult. (5423->448), div. (0->0), fcn. (4094->8), ass. (0->218)
t182 = pkin(8) + qJ(2);
t176 = sin(t182);
t277 = pkin(2) * qJD(2);
t239 = t176 * t277;
t179 = qJ(3) + t182;
t171 = sin(t179);
t172 = cos(t179);
t117 = rSges(4,1) * t171 + rSges(4,2) * t172;
t183 = qJD(2) + qJD(3);
t267 = t117 * t183;
t91 = -t239 - t267;
t246 = qJD(5) * t183;
t108 = -qJDD(5) * t172 + t171 * t246;
t181 = pkin(9) + qJ(5);
t175 = sin(t181);
t177 = cos(t181);
t135 = rSges(6,1) * t177 - rSges(6,2) * t175;
t113 = t135 * qJD(5);
t278 = rSges(6,2) * t177;
t133 = rSges(6,1) * t175 + t278;
t186 = -pkin(7) - qJ(4);
t260 = t171 * t186;
t137 = t183 * t260;
t185 = cos(pkin(9));
t173 = pkin(4) * t185 + pkin(3);
t142 = t172 * t173;
t180 = qJDD(2) + qJDD(3);
t178 = cos(t182);
t187 = qJD(2) ^ 2;
t203 = (-qJDD(2) * t176 - t178 * t187) * pkin(2);
t250 = qJD(4) * t183;
t200 = qJDD(4) * t171 + t172 * t250 + t203;
t159 = t172 * qJ(4);
t288 = pkin(3) * t171;
t116 = -t159 + t288;
t302 = -t171 * t173 - t172 * t186;
t74 = t116 + t302;
t263 = t171 * t175;
t143 = rSges(6,2) * t263;
t262 = t171 * t177;
t242 = rSges(6,1) * t262;
t85 = -t172 * rSges(6,3) - t143 + t242;
t245 = -t116 + t74 - t85;
t248 = qJD(5) * t172;
t300 = t172 * pkin(3) + t171 * qJ(4);
t258 = t172 * t177;
t208 = rSges(6,1) * t258 + t171 * rSges(6,3);
t236 = qJD(5) * t278;
t259 = t172 * t175;
t240 = rSges(6,2) * t259;
t247 = qJD(5) * t175;
t235 = -t183 * t240 + (-rSges(6,1) * t247 - t236) * t171;
t54 = t208 * t183 + t235;
t156 = qJD(4) * t172;
t88 = t183 * t300 - t156;
t9 = -t113 * t248 + t108 * t133 + t245 * t180 + (-t88 + t137 - t54 + (t300 - t142) * t183) * t183 + t200;
t309 = -g(1) + t9;
t184 = sin(pkin(9));
t279 = rSges(5,2) * t184;
t241 = t172 * t279;
t130 = t183 * t241;
t280 = rSges(5,1) * t185;
t243 = t171 * t280;
t151 = t171 * t279;
t251 = t172 * rSges(5,3) + t151;
t89 = t243 - t251;
t271 = -t116 - t89;
t301 = -t171 * rSges(5,3) - t172 * t280;
t23 = t271 * t180 + (t301 * t183 + t130 - t88) * t183 + t200;
t308 = -g(1) + t23;
t107 = qJDD(5) * t171 + t172 * t246;
t257 = t172 * t183;
t147 = qJ(4) * t257;
t170 = pkin(2) * t178;
t289 = pkin(2) * t176;
t225 = qJDD(2) * t170 - t187 * t289;
t155 = qJD(4) * t171;
t252 = t147 + t155;
t261 = t171 * t183;
t195 = -qJDD(4) * t172 + t180 * t300 + t171 * t250 + t225 + t183 * (-pkin(3) * t261 + t252);
t249 = qJD(5) * t171;
t86 = -t240 + t208;
t62 = t142 - t260 + t86;
t285 = -t300 + t62;
t207 = rSges(6,3) * t257 + t183 * t143 - t172 * t236;
t234 = t172 * t247;
t53 = (-t177 * t261 - t234) * rSges(6,1) + t207;
t10 = -t113 * t249 - t107 * t133 + t285 * t180 + (-t147 + t53 + (t302 + t288) * t183) * t183 + t195;
t307 = -g(2) + t10;
t254 = rSges(5,3) * t257 + t183 * t151;
t90 = -t241 - t301;
t24 = t180 * t90 + t183 * (-t183 * t243 + t254) + t195;
t306 = -g(2) + t24;
t102 = rSges(4,1) * t257 - rSges(4,2) * t261;
t305 = -t102 * t183 - t117 * t180 - g(1) + t203;
t119 = t172 * rSges(4,1) - rSges(4,2) * t171;
t304 = t119 * t180 - t183 * t267 - g(2) + t225;
t69 = t300 + t90;
t303 = t69 * t183;
t169 = Icges(6,4) * t177;
t212 = -Icges(6,2) * t175 + t169;
t127 = Icges(6,1) * t175 + t169;
t109 = t183 * t116;
t299 = t183 * t89 + t109 + t254;
t124 = Icges(6,5) * t177 - Icges(6,6) * t175;
t123 = Icges(6,5) * t175 + Icges(6,6) * t177;
t196 = Icges(6,3) * t183 - t123 * qJD(5);
t205 = t212 * t172;
t81 = Icges(6,6) * t171 + t205;
t274 = t175 * t81;
t269 = Icges(6,4) * t175;
t128 = Icges(6,1) * t177 - t269;
t206 = t128 * t172;
t83 = Icges(6,5) * t171 + t206;
t215 = -t177 * t83 + t274;
t298 = -t124 * t261 + t196 * t172 + t215 * t183;
t204 = t124 * t172;
t80 = Icges(6,4) * t262 - Icges(6,2) * t263 - Icges(6,6) * t172;
t275 = t175 * t80;
t140 = Icges(6,4) * t263;
t82 = Icges(6,1) * t262 - Icges(6,5) * t172 - t140;
t216 = -t177 * t82 + t275;
t297 = t196 * t171 + (t204 + t216) * t183;
t125 = Icges(6,2) * t177 + t269;
t210 = t175 * t125 - t177 * t127;
t296 = t124 * qJD(5) + t210 * t183;
t78 = Icges(6,5) * t262 - Icges(6,6) * t263 - Icges(6,3) * t172;
t30 = -t216 * t171 - t172 * t78;
t76 = t183 * t85;
t295 = -rSges(6,1) * t234 - t183 * t74 + t109 + t155 + t207 + t76;
t282 = -Icges(6,2) * t262 - t140 + t82;
t284 = t127 * t171 + t80;
t294 = -t282 * t175 - t284 * t177;
t293 = t107 / 0.2e1;
t292 = t108 / 0.2e1;
t291 = t171 / 0.2e1;
t290 = -t172 / 0.2e1;
t287 = -t171 * t78 - t82 * t258;
t79 = Icges(6,3) * t171 + t204;
t286 = t171 * t79 + t83 * t258;
t283 = -t127 * t172 - t81;
t281 = -t125 * t172 + t83;
t222 = -t133 * t248 + t155;
t201 = t222 - t239;
t28 = t245 * t183 + t201;
t273 = t183 * t28;
t265 = t123 * t172;
t45 = -t210 * t171 - t265;
t272 = t45 * t183;
t266 = t123 * t171;
t264 = t124 * t183;
t256 = -t125 + t128;
t255 = t127 + t212;
t244 = t300 + t285;
t238 = t178 * t277;
t237 = m(2) + m(3) + m(4) + m(5);
t233 = -pkin(3) - t280;
t232 = -t249 / 0.2e1;
t231 = t249 / 0.2e1;
t230 = -t248 / 0.2e1;
t229 = t248 / 0.2e1;
t65 = t83 * t262;
t228 = t172 * t79 - t65;
t227 = -t78 + t274;
t223 = -t156 + t238;
t136 = rSges(3,1) * t178 - rSges(3,2) * t176;
t134 = rSges(3,1) * t176 + rSges(3,2) * t178;
t105 = t133 * t249;
t29 = t244 * t183 - t105 + t223;
t220 = -t29 * t171 - t28 * t172;
t31 = -t81 * t263 - t228;
t219 = t31 * t171 - t30 * t172;
t32 = -t80 * t259 - t287;
t33 = -t81 * t259 + t286;
t218 = t33 * t171 - t32 * t172;
t217 = t171 * t85 + t172 * t86;
t42 = t175 * t82 + t177 * t80;
t43 = t175 * t83 + t177 * t81;
t214 = t137 + t156 - t235;
t211 = t177 * t125 + t127 * t175;
t61 = -t85 + t302;
t202 = -t281 * t175 + t283 * t177;
t68 = t233 * t171 + t159 + t251;
t199 = (-t255 * t175 + t256 * t177) * t183;
t198 = Icges(6,5) * t183 - qJD(5) * t127;
t197 = Icges(6,6) * t183 - t125 * qJD(5);
t11 = t219 * qJD(5) + t272;
t111 = t212 * qJD(5);
t112 = t128 * qJD(5);
t46 = -t210 * t172 + t266;
t44 = t46 * t183;
t12 = t218 * qJD(5) + t44;
t50 = t197 * t171 + t183 * t205;
t52 = t198 * t171 + t183 * t206;
t16 = -t216 * qJD(5) + t175 * t52 + t177 * t50;
t49 = t197 * t172 - t212 * t261;
t51 = -t128 * t261 + t198 * t172;
t17 = -t215 * qJD(5) + t175 * t51 + t177 * t49;
t190 = -t211 * qJD(5) - t111 * t175 + t112 * t177 + t123 * t183;
t20 = t296 * t171 + t190 * t172;
t21 = t190 * t171 - t296 * t172;
t193 = (t44 + ((t31 - t65 + (t79 + t275) * t172 + t287) * t172 + t286 * t171) * qJD(5)) * t229 + (-t210 * qJD(5) + t111 * t177 + t112 * t175) * t183 + (t43 + t46) * t293 + (t42 + t45) * t292 + (-t272 + ((t227 * t172 - t286 + t33) * t172 + (t227 * t171 + t228 + t32) * t171) * qJD(5) + t11) * t232 + (t17 + t20) * t231 + (t16 + t21 + t12) * t230 + (Icges(4,3) + t211 + Icges(5,2) * t185 ^ 2 + (Icges(5,1) * t184 + 0.2e1 * Icges(5,4) * t185) * t184) * t180;
t192 = -t43 * qJD(5) - t175 * t49 + t177 * t51 + t183 * t79;
t191 = -t42 * qJD(5) - t175 * t50 + t177 * t52 + t183 * t78;
t55 = t271 * t183 + t155 - t239;
t56 = t223 + t303;
t189 = (t55 * t233 * t172 + (t55 * (-rSges(5,3) - qJ(4)) + t56 * t233) * t171) * t183;
t188 = (t28 * (-t208 - t142) + t29 * (t302 - t242)) * t183;
t100 = t133 * t172;
t99 = t133 * t171;
t92 = t119 * t183 + t238;
t39 = t217 * qJD(5) + qJD(1);
t13 = t107 * t85 - t108 * t86 + qJDD(1) + (t171 * t54 + t172 * t53) * qJD(5);
t6 = t192 * t171 - t298 * t172;
t5 = t191 * t171 - t297 * t172;
t4 = t298 * t171 + t192 * t172;
t3 = t297 * t171 + t191 * t172;
t1 = [m(6) * t13 + t237 * qJDD(1) + (-m(6) - t237) * g(3); Icges(3,3) * qJDD(2) + t193 + (t304 * (t119 + t170) + t305 * (-t117 - t289) + (-t102 - t238 + t92) * t91) * m(4) + ((t134 ^ 2 + t136 ^ 2) * qJDD(2) + g(1) * t134 - g(2) * t136) * m(3) + (t28 * (t214 - t238) + t188 + t307 * (t170 + t62) + t309 * (t61 - t289) + (-t201 + t28 - t239 + t295) * t29) * m(6) + (t55 * (t130 - t223) + t189 + t306 * (t170 + t69) + t308 * (t68 - t289) + (t55 + t147 + t299) * t56) * m(5); t193 + (t244 * t273 + t188 + t307 * t62 + t309 * t61 + (-t222 + t295) * t29 + (-t105 - t156 + t214) * t28) * m(6) + (t189 + t306 * t69 + t308 * t68 + (-t155 + t252 + t299) * t56 + (t130 + t303) * t55) * m(5) + (-t267 * t92 - t102 * t91 + (t91 * t183 + t304) * t119 + (t92 * t183 - t305) * t117) * m(4); (-m(5) - m(6)) * (g(1) * t171 - g(2) * t172) + 0.2e1 * (t10 * t290 + t9 * t291) * m(6) + 0.2e1 * (t23 * t291 + t24 * t290) * m(5); t12 * t257 / 0.2e1 + (t33 * t107 + t32 * t108 + t46 * t180 + t20 * t183 + (t171 * t4 - t3 * t172) * qJD(5)) * t291 + t218 * t293 + ((t183 * t33 - t3) * t172 + (t183 * t32 + t4) * t171) * t231 + t11 * t261 / 0.2e1 + (t31 * t107 + t30 * t108 + t45 * t180 + t21 * t183 + (t6 * t171 - t172 * t5) * qJD(5)) * t290 + t219 * t292 + ((t183 * t31 - t5) * t172 + (t183 * t30 + t6) * t171) * t230 + t180 * (t43 * t171 - t42 * t172) / 0.2e1 + t183 * ((t183 * t43 - t16) * t172 + (t183 * t42 + t17) * t171) / 0.2e1 + ((-t249 * t265 + t264) * t171 + (t199 + (-t294 * t172 + (t266 + t202) * t171) * qJD(5)) * t172) * t232 + ((-t248 * t266 - t264) * t172 + (t199 + (t202 * t171 + (-t294 + t265) * t172) * qJD(5)) * t171) * t229 - t183 * ((t256 * t175 + t255 * t177) * t183 + ((t281 * t171 - t282 * t172) * t177 + (t283 * t171 + t284 * t172) * t175) * qJD(5)) / 0.2e1 + (t13 * t217 + t39 * ((t53 + t76) * t172 + (-t183 * t86 + t54) * t171) + t220 * t113 + ((-t183 * t29 - t9) * t172 + (-t10 + t273) * t171) * t133 - (-t100 * t29 + t28 * t99) * t183 - (t39 * (-t100 * t172 - t171 * t99) + t220 * t135) * qJD(5) + g(1) * t100 + g(2) * t99 - g(3) * t135) * m(6);];
tau = t1;
