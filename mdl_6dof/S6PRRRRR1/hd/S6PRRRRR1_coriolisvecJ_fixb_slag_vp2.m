% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:32
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:32:01
% EndTime: 2018-11-23 15:32:11
% DurationCPUTime: 9.89s
% Computational Cost: add. (14017->574), mult. (35711->812), div. (0->0), fcn. (27413->12), ass. (0->276)
t248 = sin(qJ(4));
t249 = sin(qJ(3));
t253 = cos(qJ(4));
t254 = cos(qJ(3));
t217 = t248 * t254 + t249 * t253;
t255 = cos(qJ(2));
t244 = sin(pkin(6));
t321 = qJD(1) * t244;
t302 = t255 * t321;
t184 = t217 * t302;
t376 = -pkin(9) - pkin(8);
t228 = t376 * t249;
t229 = t376 * t254;
t187 = t248 * t228 - t253 * t229;
t304 = qJD(3) * t376;
t222 = t249 * t304;
t223 = t254 * t304;
t388 = -qJD(4) * t187 - t222 * t248 + t253 * t223 + t184;
t216 = -t248 * t249 + t253 * t254;
t312 = qJD(4) * t253;
t313 = qJD(4) * t248;
t387 = -t216 * t302 + t253 * t222 + t248 * t223 + t228 * t312 + t229 * t313;
t243 = qJD(3) + qJD(4);
t178 = t243 * t217;
t404 = pkin(10) * t178 - t387;
t177 = t243 * t216;
t403 = -pkin(10) * t177 + t388;
t247 = sin(qJ(5));
t252 = cos(qJ(5));
t186 = t253 * t228 + t229 * t248;
t146 = -pkin(10) * t217 + t186;
t147 = pkin(10) * t216 + t187;
t281 = t252 * t146 - t147 * t247;
t402 = -qJD(5) * t281 - t403 * t247 + t252 * t404;
t314 = qJD(3) * t249;
t241 = pkin(3) * t314;
t156 = pkin(4) * t178 + t241;
t250 = sin(qJ(2));
t303 = t250 * t321;
t277 = t252 * t216 - t217 * t247;
t84 = qJD(5) * t277 + t177 * t252 - t178 * t247;
t172 = t216 * t247 + t217 * t252;
t85 = qJD(5) * t172 + t177 * t247 + t252 * t178;
t401 = pkin(5) * t85 - pkin(11) * t84 + t156 - t303;
t242 = qJD(5) + t243;
t246 = sin(qJ(6));
t251 = cos(qJ(6));
t209 = t216 * qJD(2);
t210 = t217 * qJD(2);
t278 = t209 * t247 + t252 * t210;
t134 = t242 * t251 - t246 * t278;
t374 = -t134 / 0.2e1;
t135 = t242 * t246 + t251 * t278;
t373 = -t135 / 0.2e1;
t294 = t252 * t209 - t210 * t247;
t148 = qJD(6) - t294;
t371 = -t148 / 0.2e1;
t164 = t177 * qJD(2);
t165 = t178 * qJD(2);
t76 = qJD(5) * t294 + t164 * t252 - t165 * t247;
t55 = qJD(6) * t134 + t251 * t76;
t56 = -qJD(6) * t135 - t246 * t76;
t16 = -mrSges(7,1) * t56 + mrSges(7,2) * t55;
t245 = cos(pkin(6));
t319 = qJD(1) * t254;
t232 = t245 * t319;
t224 = qJD(2) * pkin(8) + t303;
t295 = pkin(9) * qJD(2) + t224;
t182 = -t295 * t249 + t232;
t176 = qJD(3) * pkin(3) + t182;
t320 = qJD(1) * t249;
t301 = t245 * t320;
t183 = t254 * t295 + t301;
t318 = qJD(2) * t244;
t296 = qJD(1) * t318;
t293 = t255 * t296;
t275 = t254 * t293;
t276 = t249 * t293;
t61 = -t183 * t313 + t248 * (-qJD(3) * t183 - t276) + t253 * (qJD(3) * t182 + t275) + t176 * t312;
t36 = -pkin(10) * t165 + t61;
t360 = t164 * pkin(10);
t175 = t253 * t183;
t118 = t176 * t248 + t175;
t362 = pkin(10) * t209;
t107 = t118 + t362;
t322 = t252 * t107;
t173 = t248 * t183;
t117 = t253 * t176 - t173;
t200 = t210 * pkin(10);
t106 = t117 - t200;
t93 = pkin(4) * t243 + t106;
t48 = t247 * t93 + t322;
t189 = -t224 * t249 + t232;
t327 = t224 * t254;
t190 = t301 + t327;
t315 = qJD(2) * t254;
t317 = qJD(2) * t249;
t260 = (-t248 * (-pkin(9) * t317 + t189) + t253 * (-pkin(9) * t315 - t190)) * qJD(3);
t62 = -qJD(2) * t184 - qJD(4) * t118 + t260;
t9 = t247 * t36 - t252 * (t62 - t360) + t48 * qJD(5);
t400 = m(7) * t9 + t16;
t290 = mrSges(7,1) * t246 + mrSges(7,2) * t251;
t332 = t107 * t247;
t47 = t252 * t93 - t332;
t40 = -pkin(5) * t242 - t47;
t272 = t40 * t290;
t285 = Ifges(7,5) * t251 - Ifges(7,6) * t246;
t351 = Ifges(7,4) * t251;
t287 = -Ifges(7,2) * t246 + t351;
t352 = Ifges(7,4) * t246;
t289 = Ifges(7,1) * t251 - t352;
t365 = t251 / 0.2e1;
t366 = -t246 / 0.2e1;
t372 = t135 / 0.2e1;
t348 = t135 * Ifges(7,4);
t66 = t134 * Ifges(7,2) + t148 * Ifges(7,6) + t348;
t133 = Ifges(7,4) * t134;
t67 = t135 * Ifges(7,1) + t148 * Ifges(7,5) + t133;
t399 = t134 * t287 / 0.2e1 + t289 * t372 + t148 * t285 / 0.2e1 + t272 + t66 * t366 + t67 * t365;
t398 = -Ifges(4,1) / 0.2e1;
t397 = -Ifges(4,4) * t315 / 0.2e1;
t100 = t146 * t247 + t147 * t252;
t238 = -pkin(3) * t254 - pkin(2);
t193 = -pkin(4) * t216 + t238;
t94 = -pkin(5) * t277 - pkin(11) * t172 + t193;
t46 = t100 * t251 + t246 * t94;
t395 = -qJD(6) * t46 + t246 * t402 + t401 * t251;
t45 = -t100 * t246 + t251 * t94;
t394 = qJD(6) * t45 + t401 * t246 - t251 * t402;
t353 = Ifges(6,4) * t278;
t145 = Ifges(6,4) * t294;
t393 = -qJD(5) * t100 + t247 * t404 + t403 * t252;
t337 = mrSges(6,1) * t242 + mrSges(7,1) * t134 - mrSges(7,2) * t135 - mrSges(6,3) * t278;
t237 = pkin(3) * t253 + pkin(4);
t310 = qJD(5) * t252;
t311 = qJD(5) * t247;
t324 = t247 * t248;
t162 = t237 * t310 + (-t248 * t311 + (t252 * t253 - t324) * qJD(4)) * pkin(3);
t120 = t253 * t182 - t173;
t110 = -t200 + t120;
t119 = -t182 * t248 - t175;
t274 = t119 - t362;
t58 = t252 * t110 + t247 * t274;
t391 = t162 - t58;
t323 = t248 * t252;
t390 = -t110 * t247 + t252 * t274 + t237 * t311 + (t248 * t310 + (t247 * t253 + t323) * qJD(4)) * pkin(3);
t167 = -mrSges(5,1) * t209 + mrSges(5,2) * t210;
t199 = qJD(2) * t238 - t302;
t389 = -m(5) * t199 - t167;
t41 = pkin(11) * t242 + t48;
t169 = -pkin(4) * t209 + t199;
t78 = -pkin(5) * t294 - pkin(11) * t278 + t169;
t17 = -t246 * t41 + t251 * t78;
t18 = t246 * t78 + t251 * t41;
t386 = -t17 * t246 + t18 * t251;
t240 = pkin(3) * t317;
t201 = qJD(3) * t240 + t250 * t296;
t129 = pkin(4) * t165 + t201;
t77 = qJD(5) * t278 + t164 * t247 + t252 * t165;
t23 = pkin(5) * t77 - pkin(11) * t76 + t129;
t8 = t47 * qJD(5) + t252 * t36 + (-t176 * t313 - t183 * t312 - t248 * t275 - t253 * t276 + t260 - t360) * t247;
t2 = qJD(6) * t17 + t23 * t246 + t251 * t8;
t333 = qJD(6) * t18;
t3 = t23 * t251 - t246 * t8 - t333;
t385 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t55 + Ifges(7,6) * t56;
t384 = t242 * Ifges(6,6) / 0.2e1 + t353 / 0.2e1;
t113 = pkin(5) * t278 - pkin(11) * t294;
t379 = Ifges(7,5) * t373 + Ifges(7,6) * t374 + Ifges(7,3) * t371;
t283 = t17 * t251 + t18 * t246;
t339 = t242 * Ifges(6,5);
t343 = t278 * Ifges(6,1);
t98 = t145 + t339 + t343;
t383 = t169 * mrSges(6,2) + t98 / 0.2e1 + t339 / 0.2e1 + t145 / 0.2e1 - t283 * mrSges(7,3) + t399;
t381 = t55 / 0.2e1;
t380 = t56 / 0.2e1;
t378 = t77 / 0.2e1;
t344 = t294 * Ifges(6,2);
t377 = t344 / 0.2e1 + t384;
t375 = t9 * t281;
t369 = t209 / 0.2e1;
t368 = -t210 / 0.2e1;
t367 = t210 / 0.2e1;
t363 = pkin(4) * t210;
t326 = t244 * t250;
t204 = t245 * t254 - t249 * t326;
t205 = t245 * t249 + t254 * t326;
t149 = t204 * t253 - t205 * t248;
t150 = t204 * t248 + t205 * t253;
t280 = t252 * t149 - t150 * t247;
t361 = t280 * t9;
t359 = t2 * t251;
t358 = t3 * t246;
t357 = t47 * mrSges(6,3);
t356 = mrSges(5,3) * t209;
t355 = mrSges(5,3) * t210;
t354 = Ifges(4,4) * t249;
t350 = qJD(2) * pkin(2);
t345 = t278 * t48;
t340 = t210 * Ifges(5,4);
t335 = Ifges(4,5) * qJD(3);
t334 = Ifges(4,6) * qJD(3);
t331 = t294 * t246;
t330 = t294 * t251;
t325 = t244 * t255;
t203 = pkin(3) * t323 + t247 * t237;
t316 = qJD(2) * t250;
t309 = qJD(6) * t246;
t308 = qJD(6) * t251;
t307 = qJD(2) * qJD(3);
t300 = t244 * t316;
t299 = t255 * t318;
t298 = t335 / 0.2e1;
t297 = -t334 / 0.2e1;
t291 = mrSges(7,1) * t251 - mrSges(7,2) * t246;
t288 = Ifges(7,1) * t246 + t351;
t286 = Ifges(7,2) * t251 + t352;
t284 = Ifges(7,5) * t246 + Ifges(7,6) * t251;
t103 = t149 * t247 + t150 * t252;
t269 = qJD(3) * t245 + t299;
t157 = -t224 * t314 + t269 * t319;
t158 = -qJD(3) * t327 - t269 * t320;
t279 = t157 * t254 - t158 * t249;
t202 = -pkin(3) * t324 + t237 * t252;
t86 = -t103 * t246 - t251 * t325;
t273 = -t103 * t251 + t246 * t325;
t225 = -t302 - t350;
t271 = t189 * mrSges(4,3) + t317 * t398 + t397 - t335 / 0.2e1 - t225 * mrSges(4,2);
t270 = t190 * mrSges(4,3) + t334 / 0.2e1 + (t254 * Ifges(4,2) + t354) * qJD(2) / 0.2e1 - t225 * mrSges(4,1);
t88 = t113 + t363;
t265 = -qJD(6) * t283 - t358;
t264 = t265 * mrSges(7,3);
t12 = t55 * Ifges(7,4) + t56 * Ifges(7,2) + t77 * Ifges(7,6);
t13 = t55 * Ifges(7,1) + t56 * Ifges(7,4) + t77 * Ifges(7,5);
t262 = -t8 * mrSges(6,2) + mrSges(7,3) * t359 + t246 * t13 / 0.2e1 + t12 * t365 + t288 * t381 + t286 * t380 + t284 * t378 - Ifges(6,6) * t77 + Ifges(6,5) * t76 + (-mrSges(6,1) - t291) * t9 + t399 * qJD(6);
t261 = -t169 * mrSges(6,1) - t17 * mrSges(7,1) + t18 * mrSges(7,2) + t377 + 0.2e1 * t379 + t384;
t21 = mrSges(7,1) * t77 - mrSges(7,3) * t55;
t22 = -mrSges(7,2) * t77 + mrSges(7,3) * t56;
t89 = -mrSges(7,2) * t148 + mrSges(7,3) * t134;
t90 = mrSges(7,1) * t148 - mrSges(7,3) * t135;
t259 = -t246 * t21 + t251 * t22 + m(7) * (-t17 * t308 - t18 * t309 - t358 + t359) - t89 * t309 - t90 * t308;
t143 = t209 * Ifges(5,2) + t243 * Ifges(5,6) + t340;
t198 = Ifges(5,4) * t209;
t144 = t210 * Ifges(5,1) + t243 * Ifges(5,5) + t198;
t257 = Ifges(5,5) * t164 - Ifges(5,6) * t165 - (-Ifges(5,2) * t210 + t144 + t198) * t209 / 0.2e1 - (-Ifges(6,2) * t278 + t145 + t98) * t294 / 0.2e1 - t294 * t272 + t262 - t61 * mrSges(5,2) + t62 * mrSges(5,1) - t169 * (mrSges(6,1) * t278 + mrSges(6,2) * t294) + (Ifges(7,5) * t278 + t289 * t294) * t373 + (Ifges(7,6) * t278 + t287 * t294) * t374 + (Ifges(7,3) * t278 + t285 * t294) * t371 - t242 * (Ifges(6,5) * t294 - Ifges(6,6) * t278) / 0.2e1 + t294 * t357 - t199 * (mrSges(5,1) * t210 + mrSges(5,2) * t209) - t18 * (-mrSges(7,2) * t278 - mrSges(7,3) * t331) - t17 * (mrSges(7,1) * t278 - mrSges(7,3) * t330) + t278 * t377 + t278 * t379 + t143 * t367 + (Ifges(5,1) * t209 - t340) * t368 - t278 * (Ifges(6,1) * t294 - t353) / 0.2e1 - t243 * (Ifges(5,5) * t209 - Ifges(5,6) * t210) / 0.2e1 - t67 * t330 / 0.2e1 + t66 * t331 / 0.2e1 + t117 * t356;
t256 = qJD(2) ^ 2;
t227 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t315;
t226 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t317;
t214 = (mrSges(4,1) * t249 + mrSges(4,2) * t254) * t307;
t197 = pkin(11) + t203;
t196 = -pkin(5) - t202;
t192 = mrSges(5,1) * t243 - t355;
t191 = -mrSges(5,2) * t243 + t356;
t188 = t240 + t363;
t181 = -qJD(3) * t205 - t249 * t299;
t180 = qJD(3) * t204 + t254 * t299;
t138 = -mrSges(6,2) * t242 + mrSges(6,3) * t294;
t114 = mrSges(5,1) * t165 + mrSges(5,2) * t164;
t112 = -mrSges(6,1) * t294 + mrSges(6,2) * t278;
t83 = t240 + t88;
t81 = -qJD(4) * t150 - t180 * t248 + t181 * t253;
t80 = qJD(4) * t149 + t180 * t253 + t181 * t248;
t73 = Ifges(7,3) * t77;
t54 = t106 * t252 - t332;
t53 = t106 * t247 + t322;
t33 = mrSges(6,1) * t77 + mrSges(6,2) * t76;
t31 = t113 * t246 + t251 * t47;
t30 = t113 * t251 - t246 * t47;
t27 = t246 * t83 + t251 * t58;
t26 = -t246 * t58 + t251 * t83;
t25 = t246 * t88 + t251 * t54;
t24 = -t246 * t54 + t251 * t88;
t20 = qJD(5) * t103 + t247 * t80 - t252 * t81;
t19 = qJD(5) * t280 + t247 * t81 + t252 * t80;
t15 = qJD(6) * t273 - t19 * t246 + t251 * t300;
t14 = qJD(6) * t86 + t19 * t251 + t246 * t300;
t1 = [-t280 * t16 + t19 * t138 + t14 * t89 + t15 * t90 + t180 * t227 + t181 * t226 + t80 * t191 + t81 * t192 + t86 * t21 - t273 * t22 - t337 * t20 + (-t103 * t77 - t280 * t76) * mrSges(6,3) + (-t149 * t164 - t150 * t165) * mrSges(5,3) + (-t204 * t254 - t205 * t249) * mrSges(4,3) * t307 + ((-mrSges(3,2) * t256 - t114 - t214 - t33) * t255 + (-mrSges(3,1) * t256 + (t112 + t167 + qJD(2) * (-mrSges(4,1) * t254 + mrSges(4,2) * t249)) * qJD(2)) * t250) * t244 + m(7) * (t14 * t18 + t15 * t17 - t2 * t273 + t20 * t40 + t3 * t86 - t361) + m(6) * (-t361 + t103 * t8 + t19 * t48 - t20 * t47 + (-t129 * t255 + t169 * t316) * t244) + m(5) * (t117 * t81 + t118 * t80 + t149 * t62 + t150 * t61 + (t199 * t316 - t201 * t255) * t244) + m(4) * (t157 * t205 + t158 * t204 + t180 * t190 + t181 * t189 + (t225 - t302) * t300); t279 * mrSges(4,3) + m(4) * ((-t189 * t254 - t190 * t249) * qJD(3) + t279) * pkin(8) + t177 * t144 / 0.2e1 - t178 * t143 / 0.2e1 + t156 * t112 + t387 * t191 + (t117 * t388 + t118 * t387 + t186 * t62 + t187 * t61 + t199 * t241 + t201 * t238) * m(5) + t388 * t192 + (-t216 * t165 - t178 * t369) * Ifges(5,2) + (t216 * t164 - t217 * t165 + t177 * t369 - t178 * t367) * Ifges(5,4) + (-t117 * t177 - t118 * t178 - t164 * t186 - t165 * t187 + t216 * t61 - t217 * t62) * mrSges(5,3) - t281 * t16 + (-t100 * t77 - t281 * t76 - t47 * t84 - t48 * t85) * mrSges(6,3) + t393 * t337 + (t343 / 0.2e1 + t383) * t84 - (t73 / 0.2e1 - Ifges(6,4) * t76 + t129 * mrSges(6,1) - t8 * mrSges(6,3) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t77 + t385) * t277 + t199 * (mrSges(5,1) * t178 + mrSges(5,2) * t177) + t243 * (Ifges(5,5) * t177 - Ifges(5,6) * t178) / 0.2e1 + t45 * t21 + t46 * t22 + (-0.3e1 / 0.2e1 * t249 ^ 2 + 0.3e1 / 0.2e1 * t254 ^ 2) * Ifges(4,4) * t307 + t193 * t33 + t394 * t89 + (t17 * t395 + t18 * t394 + t2 * t46 + t3 * t45 - t393 * t40 - t375) * m(7) + t395 * t90 - pkin(2) * t214 + t201 * (-mrSges(5,1) * t216 + mrSges(5,2) * t217) + t238 * t114 + ((-pkin(8) * t226 - t271 + t298) * t254 + (pkin(3) * t167 - pkin(8) * t227 + t297 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t315 - t270) * t249) * qJD(3) + (-t344 / 0.2e1 - t261) * t85 + ((t226 * t249 - t227 * t254) * t255 + (-m(6) * t169 - t112 + t389) * t250 + ((t189 * t249 - t190 * t254) * t255 + (-t350 - t225) * t250) * m(4)) * t321 + (t217 * t164 + t177 * t367) * Ifges(5,1) + (t129 * mrSges(6,2) + Ifges(6,1) * t76 - Ifges(6,4) * t77 + t285 * t378 + t289 * t381 + t287 * t380 + t12 * t366 + t13 * t365 + (mrSges(6,3) + t290) * t9 + (-t2 * t246 - t251 * t3) * mrSges(7,3) + (t284 * t371 + t286 * t374 + t288 * t373 + t40 * t291 + t67 * t366 - t251 * t66 / 0.2e1 - t386 * mrSges(7,3)) * qJD(6)) * t172 + (t100 * t8 + t129 * t193 + t156 * t169 + t393 * t47 - t402 * t48 - t375) * m(6) - t402 * t138; -t188 * t112 - t157 * mrSges(4,2) + t158 * mrSges(4,1) + t391 * t138 + ((t191 * t253 - t192 * t248) * qJD(4) + (-t164 * t253 - t165 * t248) * mrSges(5,3)) * pkin(3) + (t162 * t89 + t197 * t22 + (-t17 * mrSges(7,3) - t197 * t90) * qJD(6)) * t251 + t257 - t27 * t89 - t26 * t90 - t120 * t191 - t119 * t192 + t196 * t16 + t190 * t226 - t189 * t227 + t118 * t355 + (-t162 * t90 + (-qJD(6) * t89 - t21) * t197 + (-t3 - t333) * mrSges(7,3)) * t246 + (-t202 * t76 - t203 * t77 + t345) * mrSges(6,3) + ((t298 + t397 + t271) * t254 + (t297 + (t354 / 0.2e1 + (t398 + Ifges(4,2) / 0.2e1) * t254) * qJD(2) + t389 * pkin(3) + t270) * t249) * qJD(2) - t337 * t390 + (-t169 * t188 - t202 * t9 + t203 * t8 - t390 * t47 + t391 * t48) * m(6) + ((t248 * t61 + t253 * t62 + (-t117 * t248 + t118 * t253) * qJD(4)) * pkin(3) - t117 * t119 - t118 * t120) * m(5) + ((t265 + t359) * t197 - t17 * t26 - t18 * t27 + t196 * t9 + t390 * t40 + t386 * t162) * m(7); (t192 + t355) * t118 + t264 - m(6) * (-t47 * t53 + t48 * t54) - m(7) * (t17 * t24 + t18 * t25 + t40 * t53) + t337 * t53 + t257 + (-t210 * t112 + (-t247 * t77 - t252 * t76) * mrSges(6,3) + (-t337 * t247 + (-t246 * t90 + t251 * t89 + t138) * t252 + m(7) * (t247 * t40 + t252 * t386)) * qJD(5) + (0.2e1 * t169 * t368 + t247 * t8 - t252 * t9 + (-t247 * t47 + t252 * t48) * qJD(5)) * m(6)) * pkin(4) + t259 * (pkin(4) * t247 + pkin(11)) - t54 * t138 - t25 * t89 - t24 * t90 - t117 * t191 + mrSges(6,3) * t345 + t400 * (-pkin(4) * t252 - pkin(5)); t264 + t337 * t48 + (t357 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t278 - t383) * t294 + t262 + t259 * pkin(11) + (t48 * mrSges(6,3) + t261) * t278 - t47 * t138 - m(7) * (t17 * t30 + t18 * t31 + t40 * t48) - t31 * t89 - t30 * t90 - t400 * pkin(5); t73 - t40 * (mrSges(7,1) * t135 + mrSges(7,2) * t134) + (Ifges(7,1) * t134 - t348) * t373 + t66 * t372 + (Ifges(7,5) * t134 - Ifges(7,6) * t135) * t371 - t17 * t89 + t18 * t90 + (t134 * t17 + t135 * t18) * mrSges(7,3) + (-Ifges(7,2) * t135 + t133 + t67) * t374 + t385;];
tauc  = t1(:);
