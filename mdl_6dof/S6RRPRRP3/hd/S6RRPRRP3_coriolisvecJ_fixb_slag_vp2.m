% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:47:21
% EndTime: 2019-03-09 11:48:01
% DurationCPUTime: 20.31s
% Computational Cost: add. (12849->638), mult. (32819->858), div. (0->0), fcn. (24218->8), ass. (0->285)
t254 = cos(qJ(2));
t338 = -qJ(3) - pkin(7);
t239 = t338 * t254;
t233 = qJD(1) * t239;
t319 = sin(pkin(10));
t219 = t319 * t233;
t251 = sin(qJ(2));
t238 = t338 * t251;
t232 = qJD(1) * t238;
t224 = qJD(2) * pkin(2) + t232;
t320 = cos(pkin(10));
t171 = t224 * t320 + t219;
t163 = -qJD(2) * pkin(3) - t171;
t277 = qJD(1) * t319;
t278 = qJD(1) * t320;
t216 = -t251 * t278 - t254 * t277;
t250 = sin(qJ(4));
t253 = cos(qJ(4));
t182 = qJD(2) * t253 + t216 * t250;
t118 = -t182 * pkin(4) + t163;
t183 = qJD(2) * t250 - t216 * t253;
t249 = sin(qJ(5));
t252 = cos(qJ(5));
t124 = t182 * t249 + t183 * t252;
t215 = -t251 * t277 + t254 * t278;
t209 = qJD(4) - t215;
t247 = -pkin(2) * t254 - pkin(1);
t304 = qJD(1) * t247;
t234 = qJD(3) + t304;
t143 = -t215 * pkin(3) + t216 * pkin(8) + t234;
t280 = t320 * t233;
t172 = t319 * t224 - t280;
t164 = qJD(2) * pkin(8) + t172;
t96 = t253 * t143 - t164 * t250;
t83 = -pkin(9) * t183 + t96;
t69 = pkin(4) * t209 + t83;
t97 = t143 * t250 + t164 * t253;
t84 = pkin(9) * t182 + t97;
t79 = t249 * t84;
t28 = t252 * t69 - t79;
t432 = qJ(6) * t124;
t18 = t28 - t432;
t200 = qJD(5) + t209;
t15 = pkin(5) * t200 + t18;
t81 = t252 * t84;
t29 = t249 * t69 + t81;
t275 = t252 * t182 - t183 * t249;
t398 = qJ(6) * t275;
t19 = t29 + t398;
t355 = -t200 / 0.2e1;
t368 = t124 / 0.2e1;
t369 = -t124 / 0.2e1;
t228 = t251 * t320 + t254 * t319;
t217 = t228 * qJD(2);
t201 = qJD(1) * t217;
t227 = t251 * t319 - t254 * t320;
t218 = t227 * qJD(2);
t202 = qJD(1) * t218;
t137 = qJD(4) * t182 - t202 * t253;
t138 = -qJD(4) * t183 + t202 * t250;
t57 = qJD(5) * t275 + t137 * t252 + t138 * t249;
t282 = qJD(2) * t338;
t213 = qJD(3) * t254 + t251 * t282;
t188 = t213 * qJD(1);
t214 = -t251 * qJD(3) + t254 * t282;
t256 = qJD(1) * t214;
t136 = t188 * t320 + t256 * t319;
t296 = qJD(1) * qJD(2);
t285 = t251 * t296;
t274 = pkin(2) * t285;
t140 = pkin(3) * t201 + pkin(8) * t202 + t274;
t45 = -qJD(4) * t97 - t136 * t250 + t253 * t140;
t26 = pkin(4) * t201 - pkin(9) * t137 + t45;
t299 = qJD(4) * t253;
t300 = qJD(4) * t250;
t44 = t253 * t136 + t250 * t140 + t143 * t299 - t164 * t300;
t32 = pkin(9) * t138 + t44;
t8 = -qJD(5) * t29 - t249 * t32 + t252 * t26;
t2 = pkin(5) * t201 - qJ(6) * t57 - qJD(6) * t124 + t8;
t58 = -qJD(5) * t124 - t137 * t249 + t138 * t252;
t297 = qJD(5) * t252;
t298 = qJD(5) * t249;
t7 = t249 * t26 + t252 * t32 + t69 * t297 - t298 * t84;
t3 = qJ(6) * t58 + qJD(6) * t275 + t7;
t384 = t8 * mrSges(6,1) + t2 * mrSges(7,1) - t7 * mrSges(6,2) - t3 * mrSges(7,2);
t418 = -Ifges(7,3) - Ifges(6,3);
t419 = Ifges(6,6) + Ifges(7,6);
t421 = -Ifges(7,5) - Ifges(6,5);
t391 = -t201 * t418 + t419 * t58 - t421 * t57;
t423 = Ifges(6,1) + Ifges(7,1);
t422 = Ifges(6,4) + Ifges(7,4);
t437 = t422 * t275;
t414 = t124 * t423 - t421 * t200 + t437;
t420 = Ifges(6,2) + Ifges(7,2);
t435 = t124 * t422;
t415 = t200 * t419 + t275 * t420 + t435;
t425 = -t275 / 0.2e1;
t74 = -pkin(5) * t275 + qJD(6) + t118;
t442 = t384 + t391 + (-t124 * t419 - t275 * t421) * t355 + (t124 * t19 + t15 * t275) * mrSges(7,3) + (t124 * t29 + t275 * t28) * mrSges(6,3) - t118 * (mrSges(6,1) * t124 + mrSges(6,2) * t275) - t74 * (mrSges(7,1) * t124 + mrSges(7,2) * t275) + t415 * t368 + (-t124 * t420 + t414 + t437) * t425 + (t423 * t275 - t435) * t369;
t303 = qJD(1) * t251;
t293 = pkin(2) * t303;
t154 = -pkin(3) * t216 - pkin(8) * t215 + t293;
t174 = t232 * t320 + t219;
t103 = t250 * t154 + t253 * t174;
t289 = t319 * pkin(2);
t244 = t289 + pkin(8);
t339 = pkin(9) + t244;
t281 = qJD(4) * t339;
t313 = t215 * t250;
t441 = -pkin(9) * t313 + t250 * t281 + t103;
t102 = t253 * t154 - t174 * t250;
t312 = t215 * t253;
t440 = pkin(4) * t216 + pkin(9) * t312 - t253 * t281 - t102;
t439 = t300 - t313;
t346 = -t250 / 0.2e1;
t351 = -t215 / 0.2e1;
t438 = Ifges(4,2) * t351;
t231 = t249 * t253 + t250 * t252;
t149 = t231 * t215;
t388 = qJD(4) + qJD(5);
t177 = t388 * t231;
t393 = t149 - t177;
t261 = t249 * t250 - t252 * t253;
t150 = t261 * t215;
t176 = t388 * t261;
t436 = t176 - t150;
t433 = -t420 * t58 / 0.2e1 - t422 * t57 / 0.2e1 - t419 * t201 / 0.2e1;
t225 = t339 * t250;
t226 = t339 * t253;
t168 = -t249 * t225 + t252 * t226;
t401 = -qJD(5) * t168 + t249 * t441 + t252 * t440;
t400 = -t225 * t297 - t226 * t298 + t249 * t440 - t252 * t441;
t431 = -t218 * t250 + t228 * t299;
t173 = t232 * t319 - t280;
t392 = pkin(4) * t439 - t173;
t328 = t183 * Ifges(5,4);
t112 = t182 * Ifges(5,2) + t209 * Ifges(5,6) + t328;
t272 = mrSges(5,1) * t250 + mrSges(5,2) * t253;
t429 = t112 * t346 + t163 * t272;
t426 = -Ifges(4,6) / 0.2e1;
t301 = qJD(2) * t251;
t424 = pkin(2) * t301;
t416 = -t201 * t421 + t422 * t58 + t423 * t57;
t413 = qJ(6) * t393 - qJD(6) * t261 + t400;
t412 = pkin(5) * t216 + qJ(6) * t436 - qJD(6) * t231 + t401;
t404 = t234 * mrSges(4,2);
t402 = -pkin(5) * t393 + t392;
t399 = Ifges(4,5) * qJD(2);
t160 = t261 * t228;
t105 = -mrSges(7,2) * t200 + mrSges(7,3) * t275;
t106 = -mrSges(6,2) * t200 + mrSges(6,3) * t275;
t397 = t105 + t106;
t107 = mrSges(7,1) * t200 - mrSges(7,3) * t124;
t108 = mrSges(6,1) * t200 - mrSges(6,3) * t124;
t396 = t108 + t107;
t395 = Ifges(5,5) * t137 + Ifges(5,6) * t138;
t170 = t227 * pkin(3) - t228 * pkin(8) + t247;
t180 = t238 * t319 - t239 * t320;
t175 = t253 * t180;
t117 = t250 * t170 + t175;
t335 = mrSges(4,3) * t216;
t305 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t182 + mrSges(5,2) * t183 - t335;
t390 = -t250 * t45 + t253 * t44;
t389 = t45 * mrSges(5,1) - t44 * mrSges(5,2);
t302 = qJD(1) * t254;
t317 = Ifges(3,6) * qJD(2);
t334 = Ifges(3,4) * t251;
t387 = t317 / 0.2e1 + (t254 * Ifges(3,2) + t334) * qJD(1) / 0.2e1 + pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t302);
t386 = qJD(2) * t426;
t323 = t216 * Ifges(4,4);
t385 = t323 / 0.2e1 + t438;
t267 = Ifges(5,5) * t253 - Ifges(5,6) * t250;
t332 = Ifges(5,4) * t253;
t269 = -Ifges(5,2) * t250 + t332;
t333 = Ifges(5,4) * t250;
t271 = Ifges(5,1) * t253 - t333;
t181 = Ifges(5,4) * t182;
t113 = t183 * Ifges(5,1) + t209 * Ifges(5,5) + t181;
t308 = t253 * t113;
t356 = t183 / 0.2e1;
t383 = t209 * t267 / 0.2e1 + t271 * t356 + t182 * t269 / 0.2e1 + t308 / 0.2e1 + t429;
t382 = t234 * mrSges(4,1) + t96 * mrSges(5,1) + t28 * mrSges(6,1) + t15 * mrSges(7,1) - t97 * mrSges(5,2) - t29 * mrSges(6,2) - t19 * mrSges(7,2) + t385 + t386;
t381 = -0.2e1 * pkin(1);
t379 = t57 / 0.2e1;
t378 = t58 / 0.2e1;
t371 = t275 / 0.2e1;
t367 = t137 / 0.2e1;
t366 = t138 / 0.2e1;
t358 = -t182 / 0.2e1;
t357 = -t183 / 0.2e1;
t354 = t200 / 0.2e1;
t353 = t201 / 0.2e1;
t352 = -t209 / 0.2e1;
t350 = t216 / 0.2e1;
t345 = t253 / 0.2e1;
t344 = pkin(7) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t303);
t342 = pkin(9) * t253;
t50 = -mrSges(7,2) * t201 + mrSges(7,3) * t58;
t51 = -mrSges(6,2) * t201 + mrSges(6,3) * t58;
t337 = t50 + t51;
t34 = t252 * t83 - t79;
t116 = t253 * t170 - t180 * t250;
t92 = pkin(4) * t227 - t228 * t342 + t116;
t310 = t228 * t250;
t99 = -pkin(9) * t310 + t117;
t47 = t249 * t92 + t252 * t99;
t336 = mrSges(4,3) * t215;
t205 = Ifges(4,4) * t215;
t329 = t182 * Ifges(5,6);
t327 = t183 * Ifges(5,5);
t326 = t209 * Ifges(5,3);
t324 = t216 * Ifges(4,1);
t318 = Ifges(3,5) * qJD(2);
t135 = t188 * t319 - t320 * t256;
t179 = -t320 * t238 - t239 * t319;
t316 = t135 * t179;
t291 = Ifges(5,3) * t201 + t395;
t290 = t320 * pkin(2);
t16 = -t58 * mrSges(7,1) + t57 * mrSges(7,2);
t284 = t254 * t296;
t33 = -t249 * t83 - t81;
t46 = -t249 * t99 + t252 * t92;
t279 = t201 * mrSges(4,1) - t202 * mrSges(4,2);
t153 = t213 * t320 + t214 * t319;
t155 = pkin(3) * t217 + pkin(8) * t218 + t424;
t276 = -t153 * t250 + t253 * t155;
t167 = -t252 * t225 - t226 * t249;
t245 = -t290 - pkin(3);
t273 = mrSges(5,1) * t253 - mrSges(5,2) * t250;
t270 = Ifges(5,1) * t250 + t332;
t268 = Ifges(5,2) * t253 + t333;
t266 = Ifges(5,5) * t250 + Ifges(5,6) * t253;
t265 = -t250 * t44 - t253 * t45;
t264 = -t250 * t97 - t253 * t96;
t263 = t250 * t96 - t253 * t97;
t152 = t213 * t319 - t320 * t214;
t141 = -mrSges(5,2) * t209 + mrSges(5,3) * t182;
t142 = mrSges(5,1) * t209 - mrSges(5,3) * t183;
t262 = t141 * t253 - t142 * t250;
t151 = pkin(4) * t310 + t179;
t40 = t218 * t342 + pkin(4) * t217 + (-t175 + (pkin(9) * t228 - t170) * t250) * qJD(4) + t276;
t59 = t253 * t153 + t250 * t155 + t170 * t299 - t180 * t300;
t43 = -pkin(9) * t431 + t59;
t9 = t249 * t40 + t252 * t43 + t92 * t297 - t298 * t99;
t235 = -t253 * pkin(4) + t245;
t104 = pkin(4) * t431 + t152;
t93 = -t138 * pkin(4) + t135;
t10 = -qJD(5) * t47 - t249 * t43 + t252 * t40;
t248 = Ifges(3,4) * t302;
t246 = pkin(4) * t252 + pkin(5);
t223 = Ifges(3,1) * t303 + t248 + t318;
t186 = -qJD(2) * mrSges(4,2) + t336;
t185 = pkin(5) * t261 + t235;
t166 = -mrSges(4,1) * t215 - mrSges(4,2) * t216;
t159 = t231 * t228;
t158 = t205 - t324 + t399;
t134 = -qJ(6) * t261 + t168;
t133 = -qJ(6) * t231 + t167;
t115 = -mrSges(5,2) * t201 + mrSges(5,3) * t138;
t114 = mrSges(5,1) * t201 - mrSges(5,3) * t137;
t111 = t326 + t327 + t329;
t101 = t159 * pkin(5) + t151;
t100 = pkin(4) * t183 + pkin(5) * t124;
t91 = -mrSges(5,1) * t138 + mrSges(5,2) * t137;
t86 = -mrSges(6,1) * t275 + mrSges(6,2) * t124;
t85 = -mrSges(7,1) * t275 + mrSges(7,2) * t124;
t76 = t137 * Ifges(5,1) + t138 * Ifges(5,4) + t201 * Ifges(5,5);
t75 = t137 * Ifges(5,4) + t138 * Ifges(5,2) + t201 * Ifges(5,6);
t71 = t160 * t388 + t231 * t218;
t70 = -t177 * t228 + t218 * t261;
t64 = t124 * Ifges(6,5) + Ifges(6,6) * t275 + t200 * Ifges(6,3);
t63 = t124 * Ifges(7,5) + Ifges(7,6) * t275 + t200 * Ifges(7,3);
t60 = -qJD(4) * t117 + t276;
t49 = mrSges(6,1) * t201 - mrSges(6,3) * t57;
t48 = mrSges(7,1) * t201 - mrSges(7,3) * t57;
t41 = -t71 * pkin(5) + t104;
t35 = -qJ(6) * t159 + t47;
t31 = pkin(5) * t227 + qJ(6) * t160 + t46;
t27 = -t58 * pkin(5) + t93;
t21 = t34 - t432;
t20 = t33 - t398;
t17 = -mrSges(6,1) * t58 + mrSges(6,2) * t57;
t5 = qJ(6) * t71 - qJD(6) * t159 + t9;
t4 = pkin(5) * t217 - qJ(6) * t70 + qJD(6) * t160 + t10;
t1 = [-(t404 + t158 / 0.2e1 + t205 / 0.2e1 - t324 / 0.2e1 + t264 * mrSges(5,3) + t383) * t218 + (t291 + t391 + t395) * t227 / 0.2e1 + (t382 + t64 / 0.2e1 + t63 / 0.2e1 + t327 / 0.2e1 + t329 / 0.2e1 + (Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1) * t275 + t111 / 0.2e1 + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t200 + (Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1) * t124 + t326 / 0.2e1 + t385) * t217 + (-t317 / 0.2e1 + (mrSges(3,1) * t381 - 0.3e1 / 0.2e1 * t334 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t254) * qJD(1) + (qJD(1) * mrSges(4,2) * t228 + m(4) * (t234 + t304) + t166) * pkin(2) - t387) * t301 + (-t159 * t7 + t160 * t8 - t28 * t70 + t29 * t71) * mrSges(6,3) + t93 * (mrSges(6,1) * t159 - mrSges(6,2) * t160) + t27 * (mrSges(7,1) * t159 - mrSges(7,2) * t160) + (t171 * t218 - t172 * t217 - t180 * t201) * mrSges(4,3) + (-t15 * t70 - t159 * t3 + t160 * t2 + t19 * t71) * mrSges(7,3) + (mrSges(4,1) * qJD(1) * t424 - mrSges(4,3) * t136 + Ifges(4,4) * t202 - t418 * t353 + t419 * t378 - t421 * t379 + t384 + t389) * t227 + (Ifges(5,3) / 0.2e1 + Ifges(4,2)) * t201 * t227 + (-mrSges(4,3) * t202 + t91) * t179 + (t76 * t345 + t75 * t346 - t202 * Ifges(4,1) - Ifges(4,4) * t201 + t271 * t367 + t269 * t366 + (mrSges(4,3) + t272) * t135 + t265 * mrSges(5,3) + (-t253 * t112 / 0.2e1 + t163 * t273 + t268 * t358 + t270 * t357 + t266 * t352 + t113 * t346 + t263 * mrSges(5,3)) * qJD(4)) * t228 + m(7) * (t101 * t27 + t15 * t4 + t19 * t5 + t2 * t31 + t3 * t35 + t41 * t74) + m(6) * (t10 * t28 + t104 * t118 + t151 * t93 + t29 * t9 + t46 * t8 + t47 * t7) + t153 * t186 + t151 * t17 + t59 * t141 + t60 * t142 + t116 * t114 + t117 * t115 + t118 * (-mrSges(6,1) * t71 + mrSges(6,2) * t70) + t5 * t105 + t9 * t106 + t4 * t107 + t10 * t108 + t101 * t16 + t104 * t86 + t41 * t85 + t74 * (-mrSges(7,1) * t71 + mrSges(7,2) * t70) + t31 * t48 + t46 * t49 + t35 * t50 + t47 * t51 + m(5) * (t116 * t45 + t117 * t44 + t152 * t163 + t59 * t97 + t60 * t96 + t316) + m(4) * (t136 * t180 - t152 * t171 + t153 * t172 + t316) + t305 * t152 + t247 * t279 + t414 * t70 / 0.2e1 + t415 * t71 / 0.2e1 - t416 * t160 / 0.2e1 + (-t159 * t419 + t160 * t421 + t228 * t267) * t353 + (t419 * t71 - t421 * t70) * t354 + (-t159 * t420 - t160 * t422) * t378 + (t420 * t71 + t422 * t70) * t371 + (-t159 * t422 - t160 * t423) * t379 + (t422 * t71 + t423 * t70) * t368 + (-Ifges(4,5) * t218 / 0.2e1 + t217 * t426 + (t223 / 0.2e1 - t344 + t318 / 0.2e1 + (mrSges(3,2) * t381 + 0.3e1 / 0.2e1 * Ifges(3,4) * t254) * qJD(1)) * t254) * qJD(2) + t159 * t433; t400 * t106 + (t118 * t392 + t167 * t8 + t168 * t7 + t235 * t93 + t28 * t401 + t29 * t400) * m(6) + t401 * t108 + t402 * t85 + (-m(5) * t163 - t305) * t173 + (m(5) * t390 - t250 * t114 + t253 * t115 - t141 * t300 - t142 * t299) * t244 + t392 * t86 + t387 * t303 + (m(5) * t244 * t264 + t383) * qJD(4) + (-mrSges(7,1) * t393 - mrSges(7,2) * t436) * t74 + (-mrSges(6,1) * t393 - mrSges(6,2) * t436) * t118 + (t15 * t436 + t19 * t393 - t2 * t231 - t261 * t3) * mrSges(7,3) + (-t231 * t8 - t261 * t7 + t28 * t436 + t29 * t393) * mrSges(6,3) + (Ifges(4,1) * t350 - t404 - t399 / 0.2e1 + t267 * t352 + t271 * t357 + t269 * t358 - t429) * t215 + (-t439 * t97 + (-t299 + t312) * t96 + t390) * mrSges(5,3) + t171 * t336 + ((-t135 * t320 + t136 * t319) * pkin(2) + t171 * t173 - t172 * t174 - t234 * t293) * m(4) + (t323 + t111 + t64 + t63) * t350 + (-Ifges(5,5) * t357 - Ifges(5,6) * t358 - Ifges(5,3) * t352 + t355 * t418 + t369 * t421 - t419 * t425 + t382 + t386 + t438) * t216 + t27 * (mrSges(7,1) * t261 + mrSges(7,2) * t231) + (-t231 * t421 - t261 * t419 + t266) * t353 + (t231 * t422 - t261 * t420) * t378 + (t231 * t423 - t261 * t422) * t379 + (-t149 * t420 - t150 * t422) * t425 + (-t176 / 0.2e1 + t150 / 0.2e1) * t414 + (t176 * t421 - t177 * t419) * t354 + (-t176 * t422 - t177 * t420) * t371 + (-t176 * t423 - t177 * t422) * t368 + (-t177 / 0.2e1 + t149 / 0.2e1) * t415 + (-t201 * t289 + t202 * t290) * mrSges(4,3) - (-Ifges(3,2) * t303 + t223 + t248) * t302 / 0.2e1 + (m(5) * t245 - mrSges(4,1) - t273) * t135 + (pkin(1) * (mrSges(3,1) * t251 + mrSges(3,2) * t254) - t251 * (Ifges(3,1) * t254 - t334) / 0.2e1) * qJD(1) ^ 2 - m(5) * (t102 * t96 + t103 * t97) + t93 * (mrSges(6,1) * t261 + mrSges(6,2) * t231) + (t205 + t308 + t158) * t351 + (-mrSges(3,1) * t284 + mrSges(3,2) * t285) * pkin(7) - Ifges(4,6) * t201 - Ifges(4,5) * t202 + t185 * t16 - t174 * t186 + t167 * t49 + t168 * t51 - t103 * t141 - t102 * t142 - t136 * mrSges(4,2) + t133 * t48 + t134 * t50 - t172 * t335 - (Ifges(3,5) * t254 - Ifges(3,6) * t251) * t296 / 0.2e1 - t166 * t293 - Ifges(3,6) * t285 + t302 * t344 + t75 * t345 + t412 * t107 + t413 * t105 + (t133 * t2 + t134 * t3 + t15 * t412 + t185 * t27 + t19 * t413 + t402 * t74) * m(7) + t416 * t231 / 0.2e1 + t235 * t17 + t245 * t91 + (-t149 * t419 + t150 * t421) * t355 + Ifges(3,5) * t284 + (-t149 * t422 - t150 * t423) * t369 + t250 * t76 / 0.2e1 + t268 * t366 + t270 * t367 + t261 * t433; t253 * t114 + t250 * t115 + t337 * t231 - (t48 + t49) * t261 + t262 * qJD(4) + (-t186 - t262) * t215 + (t85 + t86 + t305) * t216 + t279 - t397 * t436 + t396 * t393 + (t15 * t393 - t19 * t436 - t2 * t261 + t216 * t74 + t231 * t3) * m(7) + (t118 * t216 + t231 * t7 - t261 * t8 + t28 * t393 - t29 * t436) * m(6) + (t163 * t216 - t209 * t263 - t265) * m(5) + (-t171 * t216 - t172 * t215 + t274) * m(4); t389 - m(6) * (t28 * t33 + t29 * t34) + (-Ifges(5,2) * t183 + t113 + t181) * t358 + (-t183 * t86 + t252 * t49 + t337 * t249 + (-t249 * t396 + t252 * t397) * qJD(5) + (-t118 * t183 + t249 * t7 + t252 * t8 - t28 * t298 + t29 * t297) * m(6)) * pkin(4) + (-t100 * t74 - t15 * t20 - t19 * t21 + t2 * t246 + (-t15 * t298 + t19 * t297 + t249 * t3) * pkin(4)) * m(7) + t291 - t163 * (mrSges(5,1) * t183 + mrSges(5,2) * t182) + t97 * t142 - t96 * t141 - t21 * t105 - t34 * t106 - t20 * t107 - t33 * t108 - t100 * t85 + (t182 * t96 + t183 * t97) * mrSges(5,3) + t246 * t48 + (Ifges(5,5) * t182 - Ifges(5,6) * t183) * t352 + t112 * t356 + (Ifges(5,1) * t182 - t328) * t357 + t442; (-(-t15 + t18) * t19 + (-t124 * t74 + t2) * pkin(5)) * m(7) + (-t124 * t85 + t48) * pkin(5) - t18 * t105 - t28 * t106 + t19 * t107 + t29 * t108 + t442; -t275 * t105 + t124 * t107 + 0.2e1 * (t27 / 0.2e1 + t19 * t425 + t15 * t368) * m(7) + t16;];
tauc  = t1(:);
