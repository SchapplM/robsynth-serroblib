% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:19
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:19:24
% EndTime: 2018-11-23 15:19:35
% DurationCPUTime: 11.15s
% Computational Cost: add. (7795->662), mult. (20946->918), div. (0->0), fcn. (15924->12), ass. (0->325)
t217 = sin(pkin(7));
t223 = sin(qJ(3));
t313 = qJD(3) * t223;
t293 = t217 * t313;
t207 = pkin(3) * t293;
t227 = cos(qJ(3));
t260 = pkin(10) * t223 - qJ(4) * t227;
t311 = qJD(4) * t223;
t236 = (qJD(3) * t260 - t311) * t217;
t117 = t207 + t236;
t335 = t217 * t223;
t211 = pkin(9) * t335;
t219 = cos(pkin(7));
t300 = -pkin(2) * t227 - pkin(3);
t120 = pkin(4) * t335 + t211 + (-pkin(10) + t300) * t219;
t229 = -pkin(3) - pkin(10);
t336 = qJ(4) * t223;
t248 = t227 * t229 - t336;
t132 = (-pkin(2) + t248) * t217;
t331 = t219 * t223;
t214 = pkin(2) * t331;
t334 = t217 * t227;
t375 = pkin(4) + pkin(9);
t138 = (t334 * t375 + t214) * qJD(3);
t224 = sin(qJ(2));
t325 = t224 * t227;
t228 = cos(qJ(2));
t326 = t223 * t228;
t246 = t219 * t325 + t326;
t218 = sin(pkin(6));
t319 = qJD(1) * t218;
t140 = t246 * t319;
t222 = sin(qJ(5));
t226 = cos(qJ(5));
t299 = t224 * t319;
t285 = t217 * t299;
t309 = qJD(5) * t226;
t310 = qJD(5) * t222;
t392 = t120 * t309 - t132 * t310 + (t117 - t285) * t226 + (t138 - t140) * t222;
t312 = qJD(3) * t227;
t291 = t219 * t312;
t208 = pkin(2) * t291;
t215 = t219 * qJD(4);
t119 = -t293 * t375 + t208 + t215;
t284 = t219 * t299;
t298 = t228 * t319;
t141 = -t223 * t284 + t227 * t298;
t404 = t119 - t141;
t292 = t217 * t312;
t403 = -pkin(11) * t292 - t392;
t173 = t219 * t222 + t226 * t334;
t290 = t222 * t313;
t126 = -qJD(5) * t173 + t217 * t290;
t305 = t222 * t334;
t127 = -qJD(5) * t305 + t219 * t309 - t226 * t293;
t402 = pkin(5) * t127 - pkin(11) * t126 + t404;
t317 = qJD(2) * t217;
t296 = t223 * t317;
t200 = qJD(5) + t296;
t210 = qJD(2) * t219 + qJD(3);
t177 = pkin(9) * t317 + t299;
t187 = qJD(2) * pkin(2) + t298;
t220 = cos(pkin(6));
t304 = t220 * t334;
t189 = qJD(1) * t304;
t330 = t219 * t227;
t321 = t187 * t330 + t189;
t86 = -t223 * (pkin(4) * t317 + t177) + t321;
t65 = t210 * t229 + qJD(4) - t86;
t318 = qJD(1) * t220;
t209 = t219 * t318;
t93 = t209 + (qJD(2) * t248 - t187) * t217;
t34 = -t222 * t93 + t226 * t65;
t24 = -pkin(5) * t200 - t34;
t221 = sin(qJ(6));
t225 = cos(qJ(6));
t35 = t222 * t65 + t226 * t93;
t25 = pkin(11) * t200 + t35;
t316 = qJD(2) * t227;
t295 = t217 * t316;
t152 = -t210 * t222 - t226 * t295;
t280 = t222 * t295;
t153 = t210 * t226 - t280;
t203 = t210 * qJ(4);
t247 = t227 * t177 + t187 * t331;
t297 = t223 * t318;
t99 = t217 * t297 + t247;
t87 = pkin(4) * t295 + t99;
t74 = t203 + t87;
t36 = -pkin(5) * t152 - pkin(11) * t153 + t74;
t11 = -t221 * t25 + t225 * t36;
t12 = t221 * t36 + t225 * t25;
t259 = t11 * t225 + t12 * t221;
t268 = mrSges(7,1) * t221 + mrSges(7,2) * t225;
t363 = t225 / 0.2e1;
t364 = -t221 / 0.2e1;
t105 = -t153 * t221 + t200 * t225;
t149 = qJD(6) - t152;
t106 = t153 * t225 + t200 * t221;
t357 = Ifges(7,4) * t106;
t38 = Ifges(7,2) * t105 + Ifges(7,6) * t149 + t357;
t104 = Ifges(7,4) * t105;
t39 = Ifges(7,1) * t106 + Ifges(7,5) * t149 + t104;
t401 = -t259 * mrSges(7,3) + t24 * t268 + t363 * t39 + t364 * t38;
t397 = -t210 / 0.2e1;
t400 = -t295 / 0.2e1;
t393 = mrSges(5,1) + mrSges(4,3);
t399 = mrSges(5,2) - mrSges(4,1);
t263 = Ifges(7,5) * t225 - Ifges(7,6) * t221;
t355 = Ifges(7,4) * t225;
t265 = -Ifges(7,2) * t221 + t355;
t356 = Ifges(7,4) * t221;
t267 = Ifges(7,1) * t225 - t356;
t367 = t149 / 0.2e1;
t371 = t106 / 0.2e1;
t373 = t105 / 0.2e1;
t398 = t263 * t367 + t265 * t373 + t267 * t371 + t401;
t396 = -t317 / 0.2e1;
t389 = t222 * t120 + t226 * t132;
t69 = pkin(11) * t335 + t389;
t176 = pkin(9) * t334 + t214;
t156 = -t219 * qJ(4) - t176;
t131 = pkin(4) * t334 - t156;
t174 = t219 * t226 - t305;
t75 = pkin(5) * t173 - pkin(11) * t174 + t131;
t26 = -t221 * t69 + t225 * t75;
t395 = qJD(6) * t26 + t402 * t221 - t225 * t403;
t27 = t221 * t75 + t225 * t69;
t394 = -qJD(6) * t27 + t221 * t403 + t402 * t225;
t333 = t218 * t224;
t294 = qJD(2) * t333;
t278 = qJD(1) * t294;
t315 = qJD(3) * t177;
t391 = t223 * (t219 * t278 + t315);
t98 = t177 * t223 - t321;
t390 = -qJD(4) - t98;
t168 = t176 * qJD(3);
t388 = t140 - t168;
t147 = Ifges(6,4) * t152;
t354 = Ifges(6,5) * t200;
t359 = Ifges(6,1) * t153;
t79 = t147 + t354 + t359;
t387 = t34 * mrSges(6,3) - t79 / 0.2e1 - t354 / 0.2e1 - t74 * mrSges(6,2) - t147 / 0.2e1;
t111 = -t210 * t310 + (-t227 * t309 + t290) * t317;
t314 = qJD(3) * t217;
t288 = qJD(2) * t314;
t277 = t223 * t288;
t112 = -qJD(5) * t280 + t210 * t309 - t226 * t277;
t332 = t218 * t228;
t302 = t227 * t332;
t279 = qJD(2) * t302;
t301 = qJD(1) * t279 + qJD(3) * t189 + t187 * t291;
t286 = -t210 * qJD(4) - t301;
t42 = (-t315 + (-pkin(4) * t314 - t284) * qJD(2)) * t223 - t286;
t21 = pkin(5) * t112 - pkin(11) * t111 + t42;
t276 = t227 * t288;
t243 = t246 * qJD(2);
t239 = t218 * t243;
t53 = qJD(1) * t239 + ((pkin(4) * t316 + t297) * t217 + t247) * qJD(3);
t320 = pkin(3) * t277 + t217 * t278;
t90 = qJD(2) * t236 + t320;
t9 = t222 * t53 + t226 * t90 + t65 * t309 - t310 * t93;
t7 = pkin(11) * t276 + t9;
t1 = qJD(6) * t11 + t21 * t221 + t225 * t7;
t2 = -qJD(6) * t12 + t21 * t225 - t221 * t7;
t271 = t1 * t225 - t2 * t221;
t10 = -qJD(5) * t35 - t222 * t90 + t226 * t53;
t386 = t10 * mrSges(6,1) - t9 * mrSges(6,2) + Ifges(6,5) * t111 - Ifges(6,6) * t112;
t385 = Ifges(4,4) * t400 + Ifges(4,5) * t397;
t29 = -qJD(5) * t389 - t117 * t222 + t138 * t226;
t383 = Ifges(4,5) / 0.2e1;
t382 = Ifges(5,5) / 0.2e1;
t381 = Ifges(6,2) / 0.2e1;
t49 = -qJD(6) * t106 - t111 * t221 + t225 * t276;
t43 = Ifges(7,6) * t49;
t48 = qJD(6) * t105 + t111 * t225 + t221 * t276;
t44 = Ifges(7,5) * t48;
t13 = Ifges(7,3) * t112 + t43 + t44;
t380 = t13 / 0.2e1;
t15 = Ifges(7,1) * t48 + Ifges(7,4) * t49 + Ifges(7,5) * t112;
t379 = t15 / 0.2e1;
t378 = -t38 / 0.2e1;
t377 = t48 / 0.2e1;
t376 = t49 / 0.2e1;
t374 = -t105 / 0.2e1;
t372 = -t106 / 0.2e1;
t370 = t112 / 0.2e1;
t368 = -t149 / 0.2e1;
t366 = -t173 / 0.2e1;
t365 = t174 / 0.2e1;
t16 = -mrSges(7,1) * t49 + mrSges(7,2) * t48;
t94 = mrSges(6,1) * t276 - mrSges(6,3) * t111;
t360 = t94 - t16;
t358 = Ifges(6,4) * t153;
t353 = Ifges(7,5) * t106;
t352 = Ifges(6,2) * t152;
t351 = Ifges(6,6) * t200;
t350 = Ifges(7,6) * t105;
t349 = Ifges(7,3) * t149;
t348 = t111 * Ifges(6,1);
t347 = t111 * Ifges(6,4);
t346 = t112 * Ifges(6,4);
t287 = t219 * t302;
t303 = t223 * t333;
t121 = -t287 + t303 - t304;
t281 = t220 * t293;
t61 = t247 * qJD(3) + (t239 + t281) * qJD(1);
t345 = t121 * t61;
t340 = t210 * Ifges(5,5);
t339 = t222 * t24;
t205 = pkin(3) * t296;
t137 = t260 * t317 + t205;
t52 = t226 * t137 + t222 * t87;
t115 = mrSges(6,1) * t200 - mrSges(6,3) * t153;
t57 = -mrSges(7,1) * t105 + mrSges(7,2) * t106;
t338 = -t57 + t115;
t162 = -mrSges(5,1) * t295 - mrSges(5,3) * t210;
t96 = -mrSges(6,1) * t152 + mrSges(6,2) * t153;
t337 = t96 - t162;
t329 = t221 * t223;
t328 = t222 * t229;
t327 = t223 * t225;
t324 = -t210 * t399 - t296 * t393;
t161 = -mrSges(4,2) * t210 + mrSges(4,3) * t295;
t323 = t161 - t162;
t164 = (mrSges(5,2) * t227 - mrSges(5,3) * t223) * t317;
t270 = -mrSges(4,1) * t227 + mrSges(4,2) * t223;
t322 = t270 * t317 + t164;
t308 = t210 / 0.2e1 - qJD(3);
t307 = 0.3e1 / 0.2e1 * Ifges(5,6) + 0.3e1 / 0.2e1 * Ifges(4,4);
t306 = -t96 - t323;
t289 = t229 * t309;
t283 = t217 * t294;
t282 = t226 * t296;
t275 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t272 = pkin(5) * t226 + pkin(11) * t222;
t269 = mrSges(7,1) * t225 - mrSges(7,2) * t221;
t266 = Ifges(7,1) * t221 + t355;
t264 = Ifges(7,2) * t225 + t356;
t262 = Ifges(7,5) * t221 + Ifges(7,6) * t225;
t261 = -pkin(3) * t227 - t336;
t258 = t11 * t221 - t12 * t225;
t32 = mrSges(7,1) * t112 - mrSges(7,3) * t48;
t33 = -mrSges(7,2) * t112 + mrSges(7,3) * t49;
t257 = -t221 * t32 + t225 * t33;
t72 = -mrSges(7,2) * t149 + mrSges(7,3) * t105;
t73 = mrSges(7,1) * t149 - mrSges(7,3) * t106;
t256 = -t221 * t72 - t225 * t73;
t255 = t222 * t34 - t226 * t35;
t245 = t219 * t326 + t325;
t122 = t218 * t245 + t220 * t335;
t172 = -t217 * t332 + t219 * t220;
t89 = t121 * t222 + t172 * t226;
t45 = t122 * t225 - t221 * t89;
t46 = t122 * t221 + t225 * t89;
t51 = -t137 * t222 + t226 * t87;
t70 = t120 * t226 - t132 * t222;
t250 = t226 * t121 - t172 * t222;
t167 = -pkin(9) * t293 + t208;
t128 = -t174 * t221 + t217 * t327;
t129 = t174 * t225 + t217 * t329;
t191 = pkin(5) * t222 - pkin(11) * t226 + qJ(4);
t151 = t191 * t221 + t225 * t328;
t150 = t191 * t225 - t221 * t328;
t242 = (-qJ(4) * t312 - t311) * t217;
t54 = t286 + t391;
t60 = t301 - t391;
t241 = -t60 * mrSges(4,2) - t54 * mrSges(5,3) + t399 * t61;
t103 = t209 + (qJD(2) * t261 - t187) * t217;
t146 = -t187 * t217 + t209;
t237 = t146 * mrSges(4,1) + Ifges(4,6) * t397 + (Ifges(4,4) * t223 + t227 * Ifges(4,2)) * t396 + t340 / 0.2e1 + (-Ifges(5,6) * t223 - t227 * Ifges(5,3)) * t317 / 0.2e1 - t103 * mrSges(5,2) - t99 * mrSges(4,3);
t235 = -t359 / 0.2e1 + t387;
t234 = t146 * mrSges(4,2) + t34 * mrSges(6,1) + t98 * mrSges(4,3) + Ifges(4,1) * t296 / 0.2e1 + Ifges(5,4) * t397 + (-t223 * Ifges(5,2) - Ifges(5,6) * t227) * t396 + t200 * Ifges(6,3) + t153 * Ifges(6,5) + t152 * Ifges(6,6) - t103 * mrSges(5,3) - t35 * mrSges(6,2) - t385;
t37 = t349 + t350 + t353;
t78 = t351 + t352 + t358;
t233 = -t350 / 0.2e1 - t353 / 0.2e1 - t349 / 0.2e1 + t351 / 0.2e1 + t12 * mrSges(7,2) - t11 * mrSges(7,1) - t74 * mrSges(6,1) + t35 * mrSges(6,3) - t37 / 0.2e1 + t78 / 0.2e1 + t358 / 0.2e1;
t232 = (-t352 / 0.2e1 - t233) * t226;
t198 = Ifges(4,5) * t276;
t197 = Ifges(5,5) * t277;
t196 = Ifges(6,3) * t276;
t178 = qJD(5) * t272 + qJD(4);
t175 = pkin(2) * t330 - t211;
t166 = -qJ(4) * t295 + t205;
t158 = t219 * t300 + t211;
t157 = (-pkin(2) + t261) * t217;
t155 = (mrSges(4,1) * t223 + mrSges(4,2) * t227) * t288;
t154 = (-mrSges(5,2) * t223 - mrSges(5,3) * t227) * t288;
t148 = -t167 - t215;
t145 = (t221 * t227 + t222 * t327) * t317;
t144 = (-t222 * t329 + t225 * t227) * t317;
t143 = t210 * t221 - t225 * t282;
t142 = t210 * t225 + t221 * t282;
t139 = t207 + t242;
t114 = -mrSges(6,2) * t200 + mrSges(6,3) * t152;
t101 = -t226 * t140 + t222 * t285;
t100 = qJD(2) * t242 + t320;
t97 = pkin(5) * t153 - pkin(11) * t152;
t95 = -mrSges(6,2) * t276 - mrSges(6,3) * t112;
t92 = -qJD(6) * t151 + t178 * t225 - t221 * t289;
t91 = qJD(6) * t150 + t178 * t221 + t225 * t289;
t84 = -t203 - t99;
t83 = -pkin(3) * t210 - t390;
t81 = t281 + (qJD(3) * t245 + t243) * t218;
t80 = -qJD(3) * t287 + t210 * t303 - t220 * t292 - t279;
t68 = -pkin(5) * t335 - t70;
t64 = -qJD(6) * t129 - t126 * t221 + t225 * t292;
t63 = qJD(6) * t128 + t126 * t225 + t221 * t292;
t62 = mrSges(6,1) * t112 + mrSges(6,2) * t111;
t59 = (-t177 + (-pkin(4) - t272) * t317) * t223 + t321;
t56 = Ifges(6,5) * t276 - t346 + t348;
t55 = -t112 * Ifges(6,2) + Ifges(6,6) * t276 + t347;
t41 = pkin(11) * t295 + t52;
t40 = -pkin(5) * t295 - t51;
t31 = qJD(5) * t89 + t222 * t283 - t226 * t81;
t30 = qJD(5) * t250 + t222 * t81 + t226 * t283;
t23 = -pkin(5) * t292 - t29;
t20 = t221 * t59 + t225 * t41;
t19 = -t221 * t41 + t225 * t59;
t18 = t221 * t97 + t225 * t34;
t17 = -t221 * t34 + t225 * t97;
t14 = Ifges(7,4) * t48 + Ifges(7,2) * t49 + Ifges(7,6) * t112;
t8 = -pkin(5) * t276 - t10;
t6 = -qJD(6) * t46 - t221 * t30 - t225 * t80;
t5 = qJD(6) * t45 - t221 * t80 + t225 * t30;
t3 = [t30 * t114 + t122 * t62 + t45 * t32 + t46 * t33 + t5 * t72 + t6 * t73 + t89 * t95 + t360 * t250 - t324 * t81 - t338 * t31 + (-mrSges(3,1) * t224 - mrSges(3,2) * t228) * qJD(2) ^ 2 * t218 + (t154 + t155) * t172 + t306 * t80 + (t322 * t333 + t393 * qJD(3) * (t121 * t227 - t122 * t223)) * t317 + m(4) * (t345 + t122 * t60 - t80 * t99 + t81 * t98 + (qJD(1) * t172 + t146) * t283) + m(5) * (t100 * t172 + t103 * t283 - t122 * t54 + t80 * t84 + t81 * t83 + t345) + m(6) * (t10 * t250 + t122 * t42 + t30 * t35 - t31 * t34 - t74 * t80 + t89 * t9) + m(7) * (t1 * t46 + t11 * t6 + t12 * t5 + t2 * t45 + t24 * t31 - t250 * t8); (t197 / 0.2e1 + t198 / 0.2e1 + t241) * t219 + t388 * t324 + (t100 * t157 + t156 * t54 + t158 * t61 + (t141 + t148) * t84 - t388 * t83 + (t139 - t285) * t103) * m(5) - m(4) * (t140 * t98 + t141 * t99 + t146 * t285) + t111 * (Ifges(6,1) * t174 - Ifges(6,4) * t173) / 0.2e1 + t153 * (Ifges(6,1) * t126 - Ifges(6,4) * t127) / 0.2e1 + t394 * t73 + (t1 * t27 + t2 * t26 + t68 * t8 + (-t101 + t23) * t24 + t395 * t12 + t394 * t11) * m(7) + t395 * t72 + (-pkin(2) * t155 - t322 * t299 + (-mrSges(5,1) * t54 + mrSges(5,2) * t100 + mrSges(4,3) * t60) * t227 + (t196 / 0.2e1 - t100 * mrSges(5,3) + t393 * t61 + t386) * t223) * t217 + t26 * t32 + t27 * t33 + t389 * t95 + t306 * t141 + (-t10 * t174 - t126 * t34 - t127 * t35 - t173 * t9) * mrSges(6,3) + t392 * t114 + ((-m(4) * pkin(2) + t270) * t285 + ((Ifges(6,5) * t365 + Ifges(6,6) * t366 + t158 * mrSges(5,1) - t175 * mrSges(4,3) + t307 * t334 + (-Ifges(5,4) + t383) * t219) * t227 + (t156 * mrSges(5,1) - t176 * mrSges(4,3) - t307 * t335 + (t382 - Ifges(4,6)) * t219 + (0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) + Ifges(6,3) / 0.2e1) * t334) * t223) * qJD(3)) * t317 + ((t84 * mrSges(5,1) + (-Ifges(4,6) / 0.2e1 + t382) * t210 + t237) * t223 + (t83 * mrSges(5,1) + (t383 - Ifges(5,4) / 0.2e1) * t210 + t234) * t227) * t314 + t338 * t101 + t200 * (Ifges(6,5) * t126 - Ifges(6,6) * t127) / 0.2e1 + (Ifges(7,4) * t129 + Ifges(7,2) * t128 + Ifges(7,6) * t173) * t376 + (Ifges(7,1) * t129 + Ifges(7,4) * t128 + Ifges(7,5) * t173) * t377 + t129 * t379 + t173 * t380 + t56 * t365 + t55 * t366 + (Ifges(7,5) * t63 + Ifges(7,6) * t64 + Ifges(7,3) * t127) * t367 + (Ifges(7,5) * t129 + Ifges(7,6) * t128 + Ifges(7,3) * t173) * t370 + (Ifges(7,1) * t63 + Ifges(7,4) * t64 + Ifges(7,5) * t127) * t371 + (Ifges(7,4) * t63 + Ifges(7,2) * t64 + Ifges(7,6) * t127) * t373 - t112 * (Ifges(6,4) * t174 - Ifges(6,2) * t173) / 0.2e1 + t152 * (Ifges(6,4) * t126 - Ifges(6,2) * t127) / 0.2e1 + m(4) * (t167 * t99 + t168 * t98 - t175 * t61 + t176 * t60) + (t10 * t70 + t131 * t42 + t389 * t9 + t404 * t74 + t392 * t35 + (t101 + t29) * t34) * m(6) + t23 * t57 + t63 * t39 / 0.2e1 + t64 * t38 / 0.2e1 + t24 * (-mrSges(7,1) * t64 + mrSges(7,2) * t63) + t68 * t16 + t70 * t94 + t29 * t115 + t119 * t96 + t126 * t79 / 0.2e1 + t127 * t37 / 0.2e1 + t11 * (mrSges(7,1) * t127 - mrSges(7,3) * t63) + t12 * (-mrSges(7,2) * t127 + mrSges(7,3) * t64) - t127 * t78 / 0.2e1 + t74 * (mrSges(6,1) * t127 + mrSges(6,2) * t126) + t128 * t14 / 0.2e1 + t8 * (-mrSges(7,1) * t128 + mrSges(7,2) * t129) + t131 * t62 + t157 * t154 + t148 * t162 + t139 * t164 + t167 * t161 + t1 * (-mrSges(7,2) * t173 + mrSges(7,3) * t128) + t2 * (mrSges(7,1) * t173 - mrSges(7,3) * t129) + t42 * (mrSges(6,1) * t173 + mrSges(6,2) * t174); (-pkin(3) * t61 - qJ(4) * t54 - t103 * t166 + t390 * t84 - t83 * t99) * m(5) + m(6) * (qJ(4) * t42 + qJD(4) * t74) + (t11 * t145 - t12 * t144) * mrSges(7,3) + m(7) * (t1 * t151 + t11 * t92 + t12 * t91 + t150 * t2) + (-t19 + t92) * t73 + t323 * t98 + t324 * t99 + ((t308 * Ifges(5,4) + qJD(3) * (Ifges(6,5) * t226 - Ifges(6,6) * t222) / 0.2e1 + (-pkin(3) * qJD(3) - t83) * mrSges(5,1) - t234 + Ifges(5,6) * t400 + t385) * t227 + (-t340 / 0.2e1 + (Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1) * t296 + t308 * Ifges(4,6) + (-qJ(4) * qJD(3) - t84) * mrSges(5,1) + (Ifges(4,2) / 0.2e1 - Ifges(5,2) / 0.2e1 + Ifges(5,3) / 0.2e1 - Ifges(4,1) / 0.2e1) * t295 + t235 * t222 + t232 - t237) * t223) * t317 + (-t9 * mrSges(6,3) + t380 - t55 / 0.2e1 + t42 * mrSges(6,1) - t347 / 0.2e1 + t44 / 0.2e1 + t43 / 0.2e1 + (m(6) * t9 + t95) * t229 + (t381 + Ifges(7,3) / 0.2e1) * t112 + t275) * t222 + (t15 * t363 + t14 * t364 - t10 * mrSges(6,3) + t42 * mrSges(6,2) + t348 / 0.2e1 - t346 / 0.2e1 + t8 * t268 + t267 * t377 + t265 * t376 + t263 * t370 + t56 / 0.2e1 + (-t1 * t221 - t2 * t225) * mrSges(7,3) + (m(6) * t10 - m(7) * t8 + t360) * t229 + (mrSges(7,3) * t258 + t225 * t378 + t24 * t269 + t262 * t368 + t264 * t374 + t266 * t372 + t364 * t39) * qJD(6)) * t226 + t197 + t198 + t337 * qJD(4) + t241 + t144 * t378 + (Ifges(7,5) * t145 + Ifges(7,6) * t144) * t368 + (Ifges(7,1) * t145 + Ifges(7,4) * t144) * t372 + (Ifges(7,4) * t145 + Ifges(7,2) * t144) * t374 + (-t20 + t91) * t72 - t40 * t57 + qJ(4) * t62 - t86 * t96 - t52 * t114 - t51 * t115 - t145 * t39 / 0.2e1 - t24 * (-mrSges(7,1) * t144 + mrSges(7,2) * t145) + t150 * t32 + t151 * t33 - t166 * t164 + (t232 + (t263 * t368 + t265 * t374 + t267 * t372 + t235 - t401) * t222 + (-m(6) * t255 + m(7) * t339 + t226 * t114 - t222 * t338) * t229) * qJD(5) - m(6) * (t34 * t51 + t35 * t52 + t74 * t86) - m(7) * (t11 * t19 + t12 * t20 + t24 * t40); -t142 * t73 - t143 * t72 - t337 * t210 + (mrSges(5,1) * t312 + t164 * t223) * t317 + (t114 * t296 + (-t221 * t73 + t225 * t72 + t114) * qJD(5) + t360) * t226 + (qJD(6) * t256 - t200 * t338 + t257 + t95) * t222 + ((-qJD(5) * t258 - t8) * t226 + (qJD(5) * t24 - qJD(6) * t259 + t271) * t222 - t11 * t142 - t12 * t143 + t296 * t339) * m(7) + (t10 * t226 - t200 * t255 - t210 * t74 + t222 * t9) * m(6) + (t103 * t296 + t210 * t84 + t61) * m(5); -pkin(5) * t16 + (-pkin(5) * t8 - t11 * t17 - t12 * t18 - t24 * t35) * m(7) - t8 * t269 + ((t381 - Ifges(6,1) / 0.2e1) * t153 + t387 - t398) * t152 + t271 * mrSges(7,3) + t196 + t338 * t35 + t386 + t233 * t153 + t398 * qJD(6) + (t257 + (-m(7) * t259 + t256) * qJD(6) + m(7) * t271) * pkin(11) + t264 * t376 + t266 * t377 + t221 * t379 + t14 * t363 + t262 * t370 - t18 * t72 - t17 * t73 - t34 * t114; -t24 * (mrSges(7,1) * t106 + mrSges(7,2) * t105) + (Ifges(7,1) * t105 - t357) * t372 + t38 * t371 + (Ifges(7,5) * t105 - Ifges(7,6) * t106) * t368 - t11 * t72 + t12 * t73 + (t105 * t11 + t106 * t12) * mrSges(7,3) + t275 + t13 + (-Ifges(7,2) * t106 + t104 + t39) * t374;];
tauc  = t3(:);
