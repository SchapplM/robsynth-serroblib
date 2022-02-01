% Calculate vector of inverse dynamics joint torques for
% S5RPRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:36
% EndTime: 2022-01-23 09:20:45
% DurationCPUTime: 6.54s
% Computational Cost: add. (11923->435), mult. (9980->590), div. (0->0), fcn. (9392->10), ass. (0->216)
t200 = qJ(1) + pkin(8);
t196 = qJ(3) + t200;
t190 = sin(t196);
t191 = cos(t196);
t150 = t191 * pkin(3) + t190 * qJ(4);
t180 = qJD(4) * t191;
t199 = qJD(1) + qJD(3);
t109 = t150 * t199 - t180;
t201 = sin(pkin(9));
t259 = qJD(5) * t199;
t116 = (qJDD(5) * t190 + t191 * t259) * t201;
t202 = cos(pkin(9));
t203 = sin(qJ(5));
t205 = cos(qJ(5));
t141 = -rSges(6,3) * t202 + (rSges(6,1) * t205 - rSges(6,2) * t203) * t201;
t156 = (-rSges(6,1) * t203 - rSges(6,2) * t205) * t201;
t147 = qJD(5) * t156;
t198 = qJDD(1) + qJDD(3);
t176 = -qJDD(5) * t202 + t198;
t177 = -qJD(5) * t202 + t199;
t194 = sin(t200);
t195 = cos(t200);
t207 = qJD(1) ^ 2;
t204 = sin(qJ(1));
t206 = cos(qJ(1));
t223 = (-qJDD(1) * t204 - t206 * t207) * pkin(1);
t211 = (-qJDD(1) * t194 - t195 * t207) * pkin(2) + t223;
t260 = qJD(4) * t199;
t210 = qJDD(4) * t190 + t191 * t260 + t211;
t296 = pkin(4) * t202;
t238 = pkin(7) * t201 + t296;
t258 = qJD(5) * t201;
t250 = t190 * t258;
t134 = t238 * t190;
t183 = t191 * qJ(4);
t148 = pkin(3) * t190 - t183;
t266 = -t148 - t134;
t278 = t191 * t199;
t272 = t202 * t205;
t129 = t190 * t272 - t191 * t203;
t273 = t202 * t203;
t130 = t190 * t205 - t191 * t273;
t92 = -qJD(5) * t129 + t130 * t199;
t128 = t190 * t273 + t191 * t205;
t131 = t190 * t203 + t191 * t272;
t93 = -qJD(5) * t128 + t131 * t199;
t241 = -rSges(6,1) * t93 - rSges(6,2) * t92;
t274 = t199 * t201;
t252 = t191 * t274;
t51 = rSges(6,3) * t252 - t241;
t271 = -t129 * rSges(6,1) + t128 * rSges(6,2);
t279 = t190 * t201;
t83 = rSges(6,3) * t279 - t271;
t16 = t147 * t250 + t116 * t141 - t176 * t83 - t177 * t51 + t266 * t198 + (-t238 * t278 - t109) * t199 + t210;
t323 = t16 - g(1);
t322 = t141 * t250 - t177 * t83;
t117 = (qJDD(5) * t191 - t190 * t259) * t201;
t276 = t191 * t202;
t277 = t191 * t201;
t135 = pkin(4) * t276 + pkin(7) * t277;
t189 = pkin(2) * t195;
t197 = t206 * pkin(1);
t192 = qJDD(1) * t197;
t297 = pkin(1) * t204;
t240 = -pkin(2) * t194 - t297;
t216 = qJDD(1) * t189 + t207 * t240 + t192;
t179 = qJD(4) * t190;
t264 = qJ(4) * t278 + t179;
t280 = t190 * t199;
t213 = t199 * (-pkin(3) * t280 + t264) + t198 * t150 + t190 * t260 + t216;
t253 = t190 * t274;
t90 = -qJD(5) * t131 + t128 * t199;
t91 = qJD(5) * t130 - t129 * t199;
t292 = t91 * rSges(6,1) + t90 * rSges(6,2);
t50 = -rSges(6,3) * t253 + t292;
t85 = t131 * rSges(6,1) + t130 * rSges(6,2) + rSges(6,3) * t277;
t17 = -t117 * t141 + t198 * t135 + t176 * t85 + t177 * t50 + (-t147 * t258 - qJDD(4)) * t191 - t199 ^ 2 * t134 + t213;
t321 = t17 - g(2);
t158 = rSges(5,2) * t252;
t291 = rSges(5,1) * t202;
t256 = t190 * t291;
t263 = rSges(5,2) * t279 + t191 * rSges(5,3);
t112 = t256 - t263;
t267 = -t148 - t112;
t307 = -rSges(5,1) * t276 - t190 * rSges(5,3);
t32 = t267 * t198 + (t199 * t307 - t109 + t158) * t199 + t210;
t320 = t32 - g(1);
t113 = -rSges(5,2) * t277 - t307;
t265 = rSges(5,2) * t253 + rSges(5,3) * t278;
t33 = -qJDD(4) * t191 + t198 * t113 + t199 * (-t199 * t256 + t265) + t213;
t319 = t33 - g(2);
t166 = rSges(4,2) * t280;
t127 = rSges(4,1) * t278 - t166;
t149 = rSges(4,1) * t190 + rSges(4,2) * t191;
t318 = -t127 * t199 - t149 * t198 - g(1) + t211;
t186 = t191 * rSges(4,1);
t151 = -rSges(4,2) * t190 + t186;
t275 = t199 * t149;
t317 = t151 * t198 - t199 * t275 - g(2) + t216;
t222 = -t296 - pkin(3) + (-rSges(6,3) - pkin(7)) * t201;
t228 = t240 * qJD(1);
t219 = t179 + t228;
t36 = t199 * t266 + t219 + t322;
t306 = t189 + t197;
t227 = t306 * qJD(1);
t218 = t227 - t180;
t270 = t135 + t150;
t303 = t191 * t141 * t258 - t177 * t85 - t199 * t270;
t37 = t218 - t303;
t315 = t36 * (t180 + t241) + t37 * (t264 + t292) + (t36 * t222 * t191 + (-t36 * qJ(4) + t37 * (-rSges(6,3) * t201 - pkin(3) - t238)) * t190) * t199 + t323 * t222 * t190;
t74 = Icges(6,5) * t129 - Icges(6,6) * t128 + Icges(6,3) * t279;
t120 = Icges(6,4) * t129;
t77 = -Icges(6,2) * t128 + Icges(6,6) * t279 + t120;
t119 = Icges(6,4) * t128;
t81 = -Icges(6,1) * t129 - Icges(6,5) * t279 + t119;
t30 = -(t203 * t77 + t205 * t81) * t201 - t202 * t74;
t314 = t183 + t271;
t313 = -t128 * t77 - t129 * t81;
t312 = t130 * t77 - t131 * t81;
t309 = t74 * t277;
t26 = t309 + t312;
t308 = t26 - t309;
t160 = t195 * rSges(3,1) - rSges(3,2) * t194;
t146 = t160 + t197;
t76 = Icges(6,5) * t131 + Icges(6,6) * t130 + Icges(6,3) * t277;
t284 = Icges(6,4) * t131;
t79 = Icges(6,2) * t130 + Icges(6,6) * t277 + t284;
t121 = Icges(6,4) * t130;
t82 = Icges(6,1) * t131 + Icges(6,5) * t277 + t121;
t25 = -t128 * t79 + t129 * t82 + t76 * t279;
t304 = -t199 * t134 + t322;
t302 = t190 * (-Icges(6,2) * t129 - t119 - t81) + t191 * (-Icges(6,2) * t131 + t121 + t82);
t301 = t116 / 0.2e1;
t300 = t117 / 0.2e1;
t299 = t190 / 0.2e1;
t298 = -t191 / 0.2e1;
t282 = Icges(6,4) * t205;
t138 = -Icges(6,6) * t202 + (-Icges(6,2) * t203 + t282) * t201;
t283 = Icges(6,4) * t203;
t139 = -Icges(6,5) * t202 + (Icges(6,1) * t205 - t283) * t201;
t153 = (-Icges(6,5) * t203 - Icges(6,6) * t205) * t201;
t142 = qJD(5) * t153;
t154 = (-Icges(6,2) * t205 - t283) * t201;
t143 = qJD(5) * t154;
t155 = (-Icges(6,1) * t203 - t282) * t201;
t144 = qJD(5) * t155;
t42 = -t142 * t202 + (-t143 * t203 + t144 * t205 + (-t138 * t205 - t139 * t203) * qJD(5)) * t201;
t137 = -Icges(6,3) * t202 + (Icges(6,5) * t205 - Icges(6,6) * t203) * t201;
t65 = -t137 * t202 + (-t138 * t203 + t139 * t205) * t201;
t295 = t65 * t176 + t42 * t177;
t52 = -t128 * t138 + t129 * t139 + t137 * t279;
t289 = t177 * t52;
t24 = t279 * t74 + t313;
t287 = t190 * t24;
t286 = t30 * t116;
t31 = -t202 * t76 + (-t203 * t79 + t205 * t82) * t201;
t285 = t31 * t117;
t106 = t151 * t199 + t227;
t281 = t106 * t149;
t103 = t113 + t150;
t269 = -t138 + t155;
t268 = t139 + t154;
t140 = t199 * t148;
t262 = t179 - t140;
t257 = m(3) + m(4) + m(5);
t27 = t130 * t79 + t131 * t82 + t76 * t277;
t248 = -pkin(3) - t291;
t247 = -t258 / 0.2e1;
t246 = t258 / 0.2e1;
t245 = t190 * t247;
t244 = t190 * t246;
t243 = t191 * t247;
t242 = t191 * t246;
t174 = rSges(2,1) * t206 - rSges(2,2) * t204;
t173 = rSges(2,1) * t204 + rSges(2,2) * t206;
t159 = rSges(3,1) * t194 + rSges(3,2) * t195;
t235 = t190 * t36 - t191 * t37;
t234 = -t190 * t50 + t191 * t51;
t233 = -t190 * t85 + t191 * t83;
t232 = t190 * (-Icges(6,5) * t128 - Icges(6,6) * t129) + t191 * (Icges(6,5) * t130 - Icges(6,6) * t131);
t226 = t201 * t232;
t225 = (t191 * t25 + t287) * t201;
t224 = (t190 * t26 + t191 * t27) * t201;
t55 = t85 + t270;
t221 = (Icges(6,1) * t130 - t284 - t79) * t191 + (-Icges(6,1) * t128 - t120 - t77) * t190;
t102 = t190 * t248 + t183 + t263;
t215 = -t140 + t219;
t105 = t228 - t275;
t53 = t130 * t138 + t131 * t139 + t137 * t277;
t43 = t53 * t177;
t10 = qJD(5) * t224 + t43;
t45 = Icges(6,5) * t93 + Icges(6,6) * t92 + Icges(6,3) * t252;
t47 = Icges(6,4) * t93 + Icges(6,2) * t92 + Icges(6,6) * t252;
t49 = Icges(6,1) * t93 + Icges(6,4) * t92 + Icges(6,5) * t252;
t13 = -t202 * t45 + (-t203 * t47 + t205 * t49 + (t203 * t81 - t205 * t77) * qJD(5)) * t201;
t44 = Icges(6,5) * t91 + Icges(6,6) * t90 - Icges(6,3) * t253;
t46 = Icges(6,4) * t91 + Icges(6,2) * t90 - Icges(6,6) * t253;
t48 = Icges(6,1) * t91 + Icges(6,4) * t90 - Icges(6,5) * t253;
t14 = -t202 * t44 + (-t203 * t46 + t205 * t48 + (-t203 * t82 - t205 * t79) * qJD(5)) * t201;
t21 = t130 * t143 + t131 * t144 + t138 * t90 + t139 * t91 + (-t137 * t280 + t142 * t191) * t201;
t22 = -t128 * t143 + t129 * t144 + t138 * t92 + t139 * t93 + (t137 * t278 + t142 * t190) * t201;
t9 = qJD(5) * t225 + t289;
t212 = (t43 + ((t24 + t27 - t313) * t191 + t308 * t190) * t258) * t245 + t286 / 0.2e1 + t285 / 0.2e1 + t52 * t301 + t53 * t300 + t295 + (t14 + t21) * t242 + (Icges(5,2) * t202 ^ 2 + (Icges(5,1) * t201 + 0.2e1 * Icges(5,4) * t202) * t201 + Icges(4,3)) * t198 + (t13 + t22 + t10) * t244 + (-t289 + ((-t25 + t308 - t312) * t191 - t287) * t258 + t9) * t243;
t66 = t199 * t267 + t219;
t67 = t103 * t199 + t218;
t209 = t66 * (t158 + t180) + t67 * (t264 + t265) + (t66 * t248 * t191 + (t66 * (-rSges(5,3) - qJ(4)) + t67 * t248) * t190) * t199;
t107 = t199 * t112;
t101 = rSges(6,1) * t130 - rSges(6,2) * t131;
t100 = -rSges(6,1) * t128 - rSges(6,2) * t129;
t40 = t233 * t258 + qJD(2);
t15 = -t116 * t85 + t117 * t83 + t234 * t258 + qJDD(2);
t6 = -t128 * t46 + t129 * t48 + t79 * t92 + t82 * t93 + (t190 * t44 + t278 * t76) * t201;
t5 = -t128 * t47 + t129 * t49 + t77 * t92 - t81 * t93 + (t190 * t45 + t278 * t74) * t201;
t4 = t130 * t46 + t131 * t48 + t79 * t90 + t82 * t91 + (t191 * t44 - t280 * t76) * t201;
t3 = t130 * t47 + t131 * t49 + t77 * t90 - t81 * t91 + (t191 * t45 - t280 * t74) * t201;
t1 = [t212 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + (t105 * t166 + (-t105 * t186 - t281) * t199 + (-t105 * t306 + t106 * t240) * qJD(1) + t317 * (t151 + t306) + t318 * (-t149 + t240)) * m(4) + ((qJDD(1) * t160 - g(2) + t192) * t146 + (-qJDD(1) * t159 + t223 + (-0.2e1 * t160 + 0.2e1 * t146 - t197) * t207 - g(1)) * (-t159 - t297)) * m(3) + ((t173 ^ 2 + t174 ^ 2) * qJDD(1) + g(1) * t173 - g(2) * t174) * m(2) + ((t240 * t37 - t306 * t36) * qJD(1) - (-t36 + t215 + t304) * t37 + t321 * (t55 + t306) + t315 + t323 * (t314 + t240)) * m(6) + (-(-t107 - t66 + t215) * t67 + (t240 * t67 - t306 * t66) * qJD(1) + t209 + t319 * (t103 + t306) + t320 * (t102 + t240)) * m(5); m(6) * t15 + t257 * qJDD(2) + (-m(6) - t257) * g(3); t212 + (-t36 * (t180 + t303) - t37 * (t262 + t304) + t321 * t55 + t315 + t323 * t314) * m(6) + (t209 - t66 * t180 - t67 * (-t107 + t262) + (t199 * t66 + t319) * t103 + t320 * t102) * m(5) + (-t105 * t127 - t106 * t275 + t281 * t199 + (t105 * t199 + t317) * t151 - t318 * t149) * m(4); (-m(5) - m(6)) * (g(1) * t190 - g(2) * t191) + 0.2e1 * (t16 * t299 + t17 * t298) * m(6) + 0.2e1 * (t298 * t33 + t299 * t32) * m(5); -t10 * t253 / 0.2e1 + (-t202 * t53 + t224) * t300 + (-t202 * t21 + ((t199 * t26 + t4) * t191 + (-t199 * t27 + t3) * t190) * t201) * t242 + (t116 * t24 + t117 * t25 + t176 * t52 + t177 * t22 + (t190 * t5 + t191 * t6) * t258) * t279 / 0.2e1 + (-t202 * t52 + t225) * t301 + (-t202 * t22 + ((t199 * t24 + t6) * t191 + (-t199 * t25 + t5) * t190) * t201) * t244 - t202 * (t286 + t285 + (t13 * t190 + t14 * t191) * t258 + t295) / 0.2e1 + t176 * (-t202 * t65 + (t190 * t30 + t191 * t31) * t201) / 0.2e1 + t177 * (-t202 * t42 + ((t199 * t30 + t14) * t191 + (-t199 * t31 + t13) * t190) * t201) / 0.2e1 + ((t130 * t268 + t131 * t269 + t153 * t277) * t177 + (t130 * t302 + t221 * t131 + t191 * t226) * t258) * t243 + ((-t128 * t268 + t129 * t269 + t153 * t279) * t177 + (-t128 * t302 + t129 * t221 + t190 * t226) * t258) * t245 - t177 * (-t202 * t153 * t177 + ((-t203 * t268 + t205 * t269) * t177 + ((-t203 * t302 + t205 * t221) * t201 - t232 * t202) * qJD(5)) * t201) / 0.2e1 + (t116 * t26 + t117 * t27 + t176 * t53 + t177 * t21 + (t190 * t3 + t191 * t4) * t258 + t199 * t9) * t277 / 0.2e1 + ((t16 * t83 - t17 * t85 + t36 * t51 - t37 * t50) * t202 + (t15 * t233 + t40 * (-t278 * t85 - t280 * t83 + t234) + t235 * t147 + ((t199 * t36 - t17) * t191 + (t199 * t37 + t16) * t190) * t141) * t201 - (-t100 * t36 + t101 * t37) * t177 - (t40 * (t100 * t191 - t101 * t190) + t235 * t156) * t258 - g(1) * t101 - g(2) * t100 - g(3) * t156) * m(6);];
tau = t1;
