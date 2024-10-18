% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR14_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR14_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR14_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:24
% EndTime: 2024-09-27 18:42:30
% DurationCPUTime: 4.05s
% Computational Cost: add. (9910->417), mult. (20353->596), div. (0->0), fcn. (14694->10), ass. (0->245)
t258 = sin(qJ(2));
t254 = cos(pkin(5));
t261 = cos(qJ(3));
t315 = t254 * t261;
t293 = qJD(3) * t315;
t262 = cos(qJ(2));
t312 = t261 * t262;
t257 = sin(qJ(3));
t316 = t254 * t257;
t329 = pkin(1) * qJD(1);
t382 = pkin(2) * t293 - (-t258 * t316 + t312) * t329;
t253 = sin(pkin(5));
t298 = t253 * (-pkin(8) - pkin(9));
t286 = t257 * t298;
t381 = qJD(3) * t286 + t382;
t269 = pkin(1) * (-t257 * t262 - t258 * t315);
t206 = qJD(1) * t269;
t244 = pkin(2) * t316;
t380 = (t261 * t298 - t244) * qJD(3) - t206;
t245 = pkin(2) * t315;
t250 = t254 * pkin(3);
t190 = t245 + t250 + t286;
t317 = t253 * t261;
t221 = pkin(8) * t317 + t244;
t242 = pkin(9) * t317;
t205 = t242 + t221;
t256 = sin(qJ(4));
t260 = cos(qJ(4));
t133 = t190 * t256 + t205 * t260;
t358 = -qJD(4) * t133 - t256 * t381 + t260 * t380;
t305 = qJD(4) * t260;
t306 = qJD(4) * t256;
t357 = t190 * t305 - t205 * t306 + t256 * t380 + t260 * t381;
t210 = (-t256 * t257 + t260 * t261) * t253;
t353 = qJD(3) + qJD(4);
t154 = t353 * t210;
t340 = pkin(10) * t154;
t379 = -t340 + t358;
t211 = (t256 * t261 + t257 * t260) * t253;
t155 = t353 * t211;
t151 = t155 * pkin(10);
t378 = t151 - t357;
t252 = qJD(1) + qJD(2);
t194 = t252 * t210;
t195 = t252 * t211;
t255 = sin(qJ(5));
t259 = cos(qJ(5));
t289 = t194 * t259 - t195 * t255;
t119 = Ifges(6,4) * t289;
t126 = t194 * t255 + t195 * t259;
t231 = pkin(2) * t252 + t262 * t329;
t342 = pkin(3) * t261;
t193 = (-t252 * t342 - t231) * t253;
t138 = -pkin(4) * t194 + t193;
t237 = t252 * t254 + qJD(3);
t233 = qJD(4) + t237;
t227 = qJD(5) + t233;
t330 = Ifges(6,4) * t126;
t134 = t252 * t154;
t135 = t252 * t155;
t50 = -qJD(5) * t126 - t134 * t255 - t135 * t259;
t375 = Ifges(6,6) * t50;
t49 = qJD(5) * t289 + t134 * t259 - t135 * t255;
t376 = Ifges(6,5) * t49;
t337 = t375 + t376;
t328 = pkin(1) * qJD(2);
t299 = t258 * t328;
t285 = qJD(1) * t299;
t276 = t254 * t285;
t300 = t258 * t329;
t319 = t252 * t253;
t222 = pkin(8) * t319 + t300;
t280 = pkin(9) * t319 + t222;
t241 = t312 * t328;
t310 = qJD(1) * t241 + t231 * t293;
t108 = (-qJD(3) * t280 - t276) * t257 + t310;
t297 = t231 * t316;
t148 = t261 * t280 + t297;
t268 = qJD(2) * t269;
t267 = qJD(1) * t268;
t109 = -qJD(3) * t148 + t267;
t216 = t231 * t315;
t147 = -t257 * t280 + t216;
t129 = pkin(3) * t237 + t147;
t33 = t108 * t260 + t109 * t256 + t129 * t305 - t148 * t306;
t19 = -pkin(10) * t135 + t33;
t142 = t260 * t148;
t82 = t129 * t256 + t142;
t34 = -qJD(4) * t82 - t108 * t256 + t109 * t260;
t20 = -pkin(10) * t134 + t34;
t339 = pkin(10) * t194;
t60 = t82 + t339;
t325 = t255 * t60;
t189 = t195 * pkin(10);
t140 = t256 * t148;
t81 = t129 * t260 - t140;
t59 = -t189 + t81;
t57 = pkin(4) * t233 + t59;
t21 = t259 * t57 - t325;
t4 = qJD(5) * t21 + t19 * t259 + t20 * t255;
t324 = t259 * t60;
t22 = t255 * t57 + t324;
t5 = -qJD(5) * t22 - t19 * t255 + t20 * t259;
t370 = t5 * mrSges(6,1) - t4 * mrSges(6,2);
t67 = Ifges(6,1) * t126 + Ifges(6,5) * t227 + t119;
t377 = t337 + t370 - (Ifges(6,5) * t289 - Ifges(6,6) * t126) * t227 / 0.2e1 - (-Ifges(6,2) * t126 + t119 + t67) * t289 / 0.2e1 - t138 * (mrSges(6,1) * t126 + mrSges(6,2) * t289) - (Ifges(6,1) * t289 - t330) * t126 / 0.2e1;
t374 = Ifges(5,5) * t134;
t373 = Ifges(5,6) * t135;
t327 = t126 * t22;
t247 = pkin(3) * t260 + pkin(4);
t303 = qJD(5) * t259;
t304 = qJD(5) * t255;
t314 = t255 * t256;
t89 = -t147 * t256 - t142;
t68 = t89 - t339;
t90 = t147 * t260 - t140;
t69 = -t189 + t90;
t372 = -t255 * t68 - t259 * t69 + t247 * t303 + (-t256 * t304 + (t259 * t260 - t314) * qJD(4)) * pkin(3);
t313 = t256 * t259;
t371 = t255 * t69 - t259 * t68 - t247 * t304 + (-t256 * t303 + (-t255 * t260 - t313) * qJD(4)) * pkin(3);
t369 = t34 * mrSges(5,1) - t33 * mrSges(5,2);
t116 = (-qJD(3) * t222 - t276) * t257 + t310;
t165 = t222 * t261 + t297;
t117 = -qJD(3) * t165 + t267;
t368 = t117 * mrSges(4,1) - t116 * mrSges(4,2);
t66 = Ifges(6,2) * t289 + Ifges(6,6) * t227 + t330;
t365 = t66 / 0.2e1;
t208 = t210 * pkin(10);
t100 = t208 + t133;
t132 = t190 * t260 - t205 * t256;
t290 = pkin(4) * t254 - pkin(10) * t211;
t96 = t132 + t290;
t53 = t100 * t259 + t255 * t96;
t364 = -qJD(5) * t53 + t255 * t378 + t259 * t379;
t52 = -t100 * t255 + t259 * t96;
t363 = qJD(5) * t52 + t255 * t379 - t259 * t378;
t359 = t21 * t289;
t248 = pkin(1) * t262 + pkin(2);
t230 = t248 * t315;
t236 = pkin(1) * t258 + pkin(8) * t253;
t291 = -pkin(9) * t253 - t236;
t163 = t257 * t291 + t230 + t250;
t229 = t248 * t316;
t197 = t236 * t261 + t229;
t173 = t242 + t197;
t99 = t163 * t256 + t173 * t260;
t307 = qJD(3) * t257;
t294 = t253 * t307;
t356 = -pkin(8) * t294 + t382;
t355 = -qJD(3) * t221 - t206;
t354 = mrSges(3,1) * t258 + mrSges(3,2) * t262;
t251 = t253 ^ 2;
t350 = t126 / 0.2e1;
t348 = t195 / 0.2e1;
t345 = t257 / 0.2e1;
t343 = pkin(3) * t257;
t341 = pkin(4) * t210;
t149 = t210 * t259 - t211 * t255;
t61 = qJD(5) * t149 + t154 * t259 - t155 * t255;
t338 = t21 * t61;
t334 = mrSges(5,3) * t194;
t333 = Ifges(4,4) * t257;
t332 = Ifges(4,4) * t261;
t331 = Ifges(5,4) * t195;
t326 = t195 * mrSges(5,3);
t323 = t261 * Ifges(4,2);
t322 = t81 * t154;
t277 = mrSges(4,1) * t257 + mrSges(4,2) * t261;
t308 = qJD(3) * t253;
t212 = t277 * t308;
t320 = t212 * t231;
t318 = t253 * t257;
t311 = -t373 + t374;
t309 = t248 * t293 + t241;
t284 = t252 * t294;
t201 = pkin(3) * t284 + t253 * t285;
t296 = t252 * t318;
t295 = t252 * t317;
t292 = t308 / 0.2e1;
t238 = pkin(3) * t294;
t139 = pkin(4) * t155 + t238;
t98 = t163 * t260 - t173 * t256;
t288 = t253 * t300;
t287 = t254 * t299;
t282 = -t294 / 0.2e1;
t281 = t261 * t292;
t225 = (-pkin(2) - t342) * t253;
t217 = (-t248 - t342) * t253;
t278 = -mrSges(4,1) * t261 + mrSges(4,2) * t257;
t87 = t290 + t98;
t91 = t208 + t99;
t39 = -t255 * t91 + t259 * t87;
t40 = t255 * t87 + t259 * t91;
t164 = -t222 * t257 + t216;
t275 = t164 * t261 + t165 * t257;
t150 = t210 * t255 + t211 * t259;
t272 = t237 * (Ifges(4,5) * t261 - Ifges(4,6) * t257);
t271 = t257 * (Ifges(4,1) * t261 - t333);
t270 = (t323 + t333) * t253;
t130 = (qJD(3) * t291 - t287) * t257 + t309;
t131 = t268 + (t261 * t291 - t229) * qJD(3);
t41 = t130 * t260 + t131 * t256 + t163 * t305 - t173 * t306;
t42 = -qJD(4) * t99 - t130 * t256 + t131 * t260;
t111 = Ifges(5,2) * t194 + Ifges(5,6) * t233 + t331;
t188 = Ifges(5,4) * t194;
t112 = Ifges(5,1) * t195 + Ifges(5,5) * t233 + t188;
t264 = mrSges(6,3) * t359 + t126 * t365 - t193 * (mrSges(5,1) * t195 + mrSges(5,2) * t194) - t233 * (Ifges(5,5) * t194 - Ifges(5,6) * t195) / 0.2e1 + t311 + t81 * t334 + t111 * t348 - t195 * (Ifges(5,1) * t194 - t331) / 0.2e1 - (-Ifges(5,2) * t195 + t112 + t188) * t194 / 0.2e1 + t369 + t377;
t169 = Ifges(4,6) * t237 + t252 * t270;
t232 = Ifges(4,4) * t295;
t170 = Ifges(4,1) * t296 + Ifges(4,5) * t237 + t232;
t226 = Ifges(4,5) * qJD(3) * t295;
t62 = -qJD(5) * t150 - t154 * t255 - t155 * t259;
t97 = pkin(4) * t135 + t201;
t263 = t227 * (Ifges(6,5) * t61 + Ifges(6,6) * t62) / 0.2e1 + t272 * t292 + (t270 * t282 + (Ifges(4,1) * t257 + t332) * t253 * t281) * t252 + t154 * t112 / 0.2e1 - t155 * t111 / 0.2e1 + t138 * (-mrSges(6,1) * t62 + mrSges(6,2) * t61) + t61 * t67 / 0.2e1 + (-mrSges(6,1) * t97 + mrSges(6,3) * t4 + Ifges(6,4) * t49 + Ifges(6,2) * t50) * t149 + t289 * (Ifges(6,4) * t61 + Ifges(6,2) * t62) / 0.2e1 + t62 * t365 + (-Ifges(4,6) * t284 + t226 + t311 + t337) * t254 / 0.2e1 + t233 * (Ifges(5,5) * t154 - Ifges(5,6) * t155) / 0.2e1 + t193 * (mrSges(5,1) * t155 + mrSges(5,2) * t154) + t194 * (Ifges(5,4) * t154 - Ifges(5,2) * t155) / 0.2e1 + (Ifges(5,1) * t154 - Ifges(5,4) * t155) * t348 + (-mrSges(5,1) * t201 + mrSges(5,3) * t33 + Ifges(5,4) * t134 - Ifges(5,2) * t135) * t210 + (mrSges(5,2) * t201 - mrSges(5,3) * t34 + Ifges(5,1) * t134 - Ifges(5,4) * t135) * t211 + (-t373 / 0.2e1 + t374 / 0.2e1 + t375 / 0.2e1 + t376 / 0.2e1 + (Ifges(4,5) * t281 + Ifges(4,6) * t282) * t252 + t368 + t369 + t370) * t254 + (mrSges(6,2) * t97 - mrSges(6,3) * t5 + Ifges(6,1) * t49 + Ifges(6,4) * t50) * t150 + (t278 * t285 + (t261 * (-Ifges(4,2) * t257 + t332) + t271) * qJD(3) * t252) * t251 + t170 * t281 + t169 * t282 + (Ifges(6,1) * t61 + Ifges(6,4) * t62) * t350 - t82 * t155 * mrSges(5,3) + t22 * t62 * mrSges(6,3) + (t116 * t317 - t117 * t318) * mrSges(4,3);
t239 = t253 * t299;
t220 = -pkin(8) * t318 + t245;
t219 = pkin(3) * t313 + t247 * t255;
t218 = -pkin(3) * t314 + t247 * t259;
t213 = t239 + t238;
t204 = t278 * t319;
t203 = -mrSges(4,2) * t237 + mrSges(4,3) * t295;
t202 = mrSges(4,1) * t237 - mrSges(4,3) * t296;
t200 = t252 * t212;
t196 = -t236 * t257 + t230;
t176 = t225 - t341;
t162 = t217 - t341;
t161 = pkin(3) * t296 + pkin(4) * t195;
t160 = mrSges(5,1) * t233 - t326;
t159 = -mrSges(5,2) * t233 + t334;
t144 = -qJD(3) * t197 + t268;
t143 = (-qJD(3) * t236 - t287) * t257 + t309;
t137 = -mrSges(5,1) * t194 + mrSges(5,2) * t195;
t115 = t139 + t239;
t105 = mrSges(6,1) * t227 - mrSges(6,3) * t126;
t104 = -mrSges(6,2) * t227 + mrSges(6,3) * t289;
t75 = mrSges(5,1) * t135 + mrSges(5,2) * t134;
t72 = -mrSges(6,1) * t289 + mrSges(6,2) * t126;
t38 = t42 - t340;
t37 = -t151 + t41;
t27 = t259 * t59 - t325;
t26 = -t255 * t59 - t324;
t14 = -mrSges(6,1) * t50 + mrSges(6,2) * t49;
t7 = -qJD(5) * t40 - t255 * t37 + t259 * t38;
t6 = qJD(5) * t39 + t255 * t38 + t259 * t37;
t1 = [t213 * t137 + t217 * t75 + (t204 * t299 - t200 * t248 - t320 + ((-t196 * t261 - t197 * t257) * t252 - t275) * qJD(3) * mrSges(4,3)) * t253 + t144 * t202 + t143 * t203 + t41 * t159 + t42 * t160 + t162 * t14 + t115 * t72 + t6 * t104 + t7 * t105 + m(6) * (t115 * t138 + t162 * t97 + t21 * t7 + t22 * t6 + t39 * t5 + t4 * t40) + m(5) * (t193 * t213 + t201 * t217 + t33 * t99 + t34 * t98 + t41 * t82 + t42 * t81) + (-t39 * t49 + t40 * t50 - t338) * mrSges(6,3) + (-t134 * t98 - t135 * t99 - t322) * mrSges(5,3) + t263 + m(4) * (t116 * t197 + t117 * t196 + t143 * t165 + t144 * t164 + (-qJD(1) * t248 - t231) * t251 * t299) + t354 * t328 * (-qJD(1) - t252); t225 * t75 + (-pkin(2) * t200 - t320 + (-t137 - t204 - t72) * t300 + (t137 * t343 + ((-t220 * t261 - t221 * t257) * t252 - t275) * mrSges(4,3)) * qJD(3)) * t253 + t176 * t14 + t139 * t72 + t364 * t105 + t363 * t104 + t356 * t203 + t355 * t202 + t263 + t358 * t160 + (-t49 * t52 + t50 * t53 - t338) * mrSges(6,3) + (-t132 * t134 - t133 * t135 - t322) * mrSges(5,3) + t357 * t159 + t354 * t329 * (-qJD(2) + t252) + (t176 * t97 + t4 * t53 + t5 * t52 + t363 * t22 + t364 * t21 + (t139 - t288) * t138) * m(6) + (t132 * t34 + t133 * t33 + t201 * t225 + t357 * t82 + t358 * t81 + (t238 - t288) * t193) * m(5) + (t116 * t221 + t117 * t220 + (-pkin(2) * t285 + t231 * t300) * t251 + t356 * t165 + t355 * t164) * m(4); t226 + t165 * t202 - t164 * t203 - t90 * t159 - t89 * t160 - t161 * t72 + t371 * t105 + t372 * t104 + (-t218 * t49 + t219 * t50 + t327) * mrSges(6,3) - m(5) * (t81 * t89 + t82 * t90) + t264 + (-t272 / 0.2e1 + t169 * t345 - Ifges(4,6) * t307 + (t231 * t277 + (t323 * t345 - t271 / 0.2e1) * t252) * t253 + (-m(5) * t193 - t137) * t343 + t275 * mrSges(4,3) - (t232 + t170) * t261 / 0.2e1) * t319 + (m(5) * (t256 * t33 + t260 * t34 + t305 * t82 - t306 * t81) + t159 * t305 - t160 * t306 + (-t134 * t260 - t135 * t256) * mrSges(5,3)) * pkin(3) + t82 * t326 + (-t138 * t161 + t21 * t371 + t5 * t218 + t4 * t219 + t22 * t372) * m(6) + t368; -t81 * t159 - t27 * t104 - t26 * t105 + t264 - m(6) * (t21 * t26 + t22 * t27) + (-t195 * t72 + (t104 * t259 - t105 * t255) * qJD(5) + (t255 * t50 - t259 * t49) * mrSges(6,3) + (-t138 * t195 - t21 * t304 + t22 * t303 + t255 * t4 + t259 * t5) * m(6)) * pkin(4) + mrSges(6,3) * t327 + (t160 + t326) * t82; t66 * t350 - t21 * t104 + t22 * t105 + (t359 + t327) * mrSges(6,3) + t377;];
tauc = t1(:);
