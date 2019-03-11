% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:37:02
% EndTime: 2019-03-09 08:37:26
% DurationCPUTime: 12.32s
% Computational Cost: add. (4980->611), mult. (12633->798), div. (0->0), fcn. (8087->6), ass. (0->268)
t217 = sin(qJ(2));
t271 = qJD(1) * qJD(2);
t252 = t217 * t271;
t373 = -t252 / 0.2e1;
t213 = sin(pkin(9));
t214 = cos(pkin(9));
t216 = sin(qJ(5));
t218 = cos(qJ(5));
t168 = t213 * t216 + t214 * t218;
t219 = cos(qJ(2));
t229 = t168 * t219;
t227 = qJD(2) * t229;
t283 = qJD(1) * t217;
t260 = t213 * t283;
t272 = t214 * qJD(2);
t165 = t260 - t272;
t258 = t214 * t283;
t273 = t213 * qJD(2);
t166 = t258 + t273;
t233 = t218 * t165 - t166 * t216;
t56 = qJD(1) * t227 + qJD(5) * t233;
t343 = t56 / 0.2e1;
t105 = t165 * t216 + t166 * t218;
t251 = t219 * t271;
t249 = t214 * t251;
t250 = t213 * t251;
t57 = qJD(5) * t105 + t216 * t249 - t218 * t250;
t341 = t57 / 0.2e1;
t371 = Ifges(4,1) + Ifges(5,1);
t355 = Ifges(6,1) + Ifges(7,1);
t354 = -Ifges(6,4) + Ifges(7,5);
t353 = Ifges(7,4) + Ifges(6,5);
t282 = qJD(1) * t219;
t257 = t214 * t282;
t259 = t213 * t282;
t136 = t216 * t257 - t218 * t259;
t274 = qJD(5) * t218;
t275 = qJD(5) * t216;
t153 = t213 * t274 - t214 * t275;
t288 = t136 - t153;
t137 = qJD(1) * t229;
t152 = t168 * qJD(5);
t287 = t137 + t152;
t372 = t354 * t341 + t355 * t343 + t353 * t373;
t366 = Ifges(5,4) + Ifges(4,5);
t206 = pkin(7) * t282;
t285 = qJ(4) * t257 - t206;
t339 = pkin(3) + pkin(4);
t109 = -t259 * t339 + t285;
t370 = qJD(4) * t213 - t109;
t369 = -qJD(2) / 0.2e1;
t368 = qJD(2) / 0.2e1;
t367 = -Ifges(4,4) + Ifges(5,5);
t365 = Ifges(7,2) + Ifges(6,3);
t352 = -Ifges(6,6) + Ifges(7,6);
t100 = Ifges(6,4) * t233;
t201 = qJD(5) + t282;
t303 = Ifges(7,5) * t233;
t351 = t105 * t355 + t201 * t353 + t100 - t303;
t296 = t213 * t218;
t169 = -t214 * t216 + t296;
t363 = -pkin(5) * t288 + qJ(6) * t287 - qJD(6) * t169 + t370;
t305 = Ifges(5,5) * t213;
t310 = Ifges(4,4) * t213;
t362 = t214 * t371 + t305 - t310;
t335 = -t233 / 0.2e1;
t336 = t233 / 0.2e1;
t334 = -t105 / 0.2e1;
t359 = t166 / 0.2e1;
t329 = -t201 / 0.2e1;
t358 = t271 / 0.2e1;
t255 = Ifges(3,6) * t369;
t356 = Ifges(3,5) * t368;
t177 = -pkin(2) * t219 - t217 * qJ(3) - pkin(1);
t160 = t177 * qJD(1);
t188 = qJD(2) * qJ(3) + t206;
t107 = t214 * t160 - t213 * t188;
t88 = pkin(3) * t282 + qJD(4) - t107;
t55 = pkin(4) * t282 - t166 * pkin(8) + t88;
t108 = t213 * t160 + t214 * t188;
t90 = -qJ(4) * t282 + t108;
t61 = t165 * pkin(8) + t90;
t20 = -t216 * t61 + t218 * t55;
t350 = qJD(6) - t20;
t276 = qJD(4) * t219;
t281 = qJD(2) * t217;
t349 = qJ(4) * t281 - t276;
t348 = t213 * t339 + pkin(7);
t21 = t216 * t55 + t218 * t61;
t291 = t214 * t219;
t270 = pkin(8) * t291;
t241 = pkin(2) * t217 - qJ(3) * t219;
t150 = qJD(2) * t241 - t217 * qJD(3);
t138 = t150 * qJD(1);
t205 = pkin(7) * t283;
t174 = (qJD(3) - t205) * qJD(2);
t85 = t214 * t138 - t213 * t174;
t48 = (-t217 * t339 - t270) * t271 - t85;
t86 = t213 * t138 + t214 * t174;
t263 = qJ(4) * t252 + t86;
t49 = (pkin(8) * t273 - qJD(4)) * t282 + t263;
t4 = -qJD(5) * t21 - t216 * t49 + t218 * t48;
t131 = pkin(7) * t291 + t213 * t177;
t124 = -qJ(4) * t219 + t131;
t324 = pkin(8) * t217;
t106 = t213 * t324 + t124;
t295 = t213 * t219;
t199 = pkin(7) * t295;
t212 = t219 * pkin(3);
t91 = pkin(4) * t219 + t199 + t212 + (-t177 - t324) * t214;
t315 = t218 * t106 + t216 * t91;
t262 = -pkin(7) * t213 - pkin(3);
t225 = -t270 + (-pkin(4) + t262) * t217;
t294 = t214 * t150;
t65 = qJD(2) * t225 - t294;
t141 = t213 * t150;
t292 = t214 * t217;
t230 = -pkin(7) * t292 + pkin(8) * t295;
t66 = qJD(2) * t230 + t141 + t349;
t9 = -qJD(5) * t315 - t216 * t66 + t218 * t65;
t176 = -qJD(2) * pkin(2) + qJD(3) + t205;
t204 = Ifges(3,4) * t282;
t304 = Ifges(5,5) * t214;
t242 = Ifges(5,3) * t213 + t304;
t309 = Ifges(4,4) * t214;
t243 = -Ifges(4,2) * t213 + t309;
t314 = mrSges(5,3) * t214;
t325 = t214 / 0.2e1;
t326 = t213 / 0.2e1;
t327 = -t213 / 0.2e1;
t82 = t165 * pkin(3) - t166 * qJ(4) + t176;
t347 = -(t213 * t90 - t214 * t88) * mrSges(5,2) - (t107 * t214 + t108 * t213) * mrSges(4,3) + (m(4) * t176 - qJD(2) * mrSges(3,1) + mrSges(4,1) * t165 + mrSges(4,2) * t166 + mrSges(3,3) * t283) * pkin(7) + t176 * (mrSges(4,1) * t213 + mrSges(4,2) * t214) + t82 * (mrSges(5,1) * t213 - t314) + Ifges(3,1) * t283 / 0.2e1 + t204 / 0.2e1 + t356 + (Ifges(5,5) * t166 - Ifges(5,6) * t282) * t326 + (Ifges(4,4) * t166 - Ifges(4,6) * t282) * t327 + t362 * t359 + (t166 * t371 - t282 * t366) * t325 + (Ifges(5,3) * t326 - Ifges(4,2) * t327 + t367 * t325 - t243 / 0.2e1 + t242 / 0.2e1) * t165;
t15 = -pkin(5) * t201 + t350;
t16 = qJ(6) * t201 + t21;
t264 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t267 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t311 = Ifges(3,4) * t217;
t346 = t264 * t233 + t267 * t105 + 0.2e1 * (Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * t165 + (Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t166 - (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t201 + t107 * mrSges(4,1) + t15 * mrSges(7,1) + t21 * mrSges(6,2) + t90 * mrSges(5,3) + t255 - (t219 * Ifges(3,2) + t311) * qJD(1) / 0.2e1 + Ifges(6,6) * t335 + Ifges(7,6) * t336 - t108 * mrSges(4,2) - t16 * mrSges(7,3) - t20 * mrSges(6,1) - t88 * mrSges(5,1) + t366 * t359 - (Ifges(4,3) + Ifges(5,2)) * t282 / 0.2e1 + t353 * t334 + t365 * t329;
t345 = Ifges(7,5) * t343 + Ifges(7,6) * t373 + Ifges(7,3) * t341;
t344 = -t56 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t341 + Ifges(6,6) * t252 / 0.2e1;
t342 = -t57 / 0.2e1;
t338 = pkin(1) * mrSges(3,1);
t337 = pkin(1) * mrSges(3,2);
t333 = t105 / 0.2e1;
t328 = t201 / 0.2e1;
t320 = -pkin(8) + qJ(3);
t40 = -mrSges(7,2) * t57 - mrSges(7,3) * t252;
t43 = mrSges(6,2) * t252 - mrSges(6,3) * t57;
t319 = t40 + t43;
t41 = -mrSges(6,1) * t252 - mrSges(6,3) * t56;
t42 = mrSges(7,1) * t252 + t56 * mrSges(7,2);
t318 = t42 - t41;
t313 = mrSges(6,3) * t233;
t72 = -mrSges(6,2) * t201 + t313;
t75 = mrSges(7,2) * t233 + mrSges(7,3) * t201;
t317 = t72 + t75;
t312 = mrSges(6,3) * t105;
t73 = mrSges(6,1) * t201 - t312;
t74 = -mrSges(7,1) * t201 + mrSges(7,2) * t105;
t316 = t73 - t74;
t172 = t241 * qJD(1);
t293 = t214 * t172;
t76 = qJD(1) * t225 - t293;
t154 = t213 * t172;
t202 = qJ(4) * t283;
t89 = qJD(1) * t230 + t154 + t202;
t25 = t216 * t76 + t218 * t89;
t308 = Ifges(5,4) * t214;
t307 = Ifges(6,4) * t105;
t306 = Ifges(4,5) * t214;
t302 = Ifges(4,6) * t213;
t301 = Ifges(5,6) * t213;
t112 = mrSges(5,1) * t165 - mrSges(5,3) * t166;
t39 = -mrSges(6,1) * t233 + mrSges(6,2) * t105;
t300 = t39 - t112;
t297 = qJD(2) * mrSges(3,2);
t132 = -t165 * mrSges(5,2) - mrSges(5,3) * t282;
t133 = mrSges(4,2) * t282 - t165 * mrSges(4,3);
t290 = t132 + t133;
t134 = -mrSges(4,1) * t282 - t166 * mrSges(4,3);
t135 = mrSges(5,1) * t282 + t166 * mrSges(5,2);
t289 = t134 - t135;
t286 = -qJ(4) * t249 - t166 * qJD(4);
t140 = mrSges(4,1) * t250 + mrSges(4,2) * t249;
t256 = t219 * t272;
t284 = -qJ(4) * t256 - qJD(4) * t292;
t280 = qJD(2) * t219;
t279 = qJD(3) * t213;
t278 = qJD(3) * t214;
t175 = -t214 * pkin(3) - t213 * qJ(4) - pkin(2);
t269 = pkin(7) * t281;
t261 = pkin(3) * t213 + pkin(7);
t130 = t214 * t177 - t199;
t155 = t214 * pkin(4) - t175;
t246 = t262 * t217;
t24 = -t216 * t89 + t218 * t76;
t35 = -t216 * t106 + t218 * t91;
t181 = t320 * t213;
t182 = t320 * t214;
t232 = t218 * t181 - t182 * t216;
t122 = t181 * t216 + t182 * t218;
t127 = -pkin(7) * t258 + t154;
t116 = -t214 * t269 + t141;
t231 = t261 * t280;
t3 = t216 * t48 + t218 * t49 + t55 * t274 - t275 * t61;
t8 = -t106 * t275 + t216 * t65 + t218 * t66 + t91 * t274;
t145 = t168 * t217;
t196 = qJ(4) * t292;
t123 = -t217 * t348 + t196;
t59 = -pkin(4) * t165 - t82;
t228 = -t216 * t316 + t218 * t317;
t1 = -qJ(6) * t252 + qJD(6) * t201 + t3;
t2 = pkin(5) * t252 - t4;
t226 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t87 = -t280 * t348 - t284;
t60 = -qJD(1) * t276 + t263;
t77 = qJD(1) * t231 + t286;
t224 = (Ifges(5,6) * t217 + t219 * t242) * t358 - (Ifges(4,6) * t217 + t219 * t243) * t271 / 0.2e1 + t77 * mrSges(5,1) - t60 * mrSges(5,2) - t86 * mrSges(4,3);
t68 = -pkin(3) * t252 - t85;
t223 = t68 * mrSges(5,2) - t85 * mrSges(4,3) - t77 * mrSges(5,3) + (t217 * t366 + t219 * t362) * t358;
t62 = t251 * t348 + t286;
t222 = t308 / 0.2e1 + t301 / 0.2e1 + t306 / 0.2e1 - t302 / 0.2e1;
t189 = mrSges(3,3) * t282 - t297;
t185 = mrSges(5,2) * t249;
t183 = mrSges(5,1) * t250;
t149 = (-mrSges(5,2) * t295 + mrSges(5,3) * t217) * t271;
t148 = -mrSges(5,1) * t252 + t185;
t147 = (mrSges(4,1) * t217 - mrSges(4,3) * t291) * t271;
t146 = (-mrSges(4,2) * t217 - mrSges(4,3) * t295) * t271;
t144 = t216 * t292 - t217 * t296;
t142 = t217 * t261 - t196;
t139 = -mrSges(5,3) * t249 + t183;
t129 = pkin(3) * t259 - t285;
t126 = pkin(7) * t260 + t293;
t125 = -t130 + t212;
t115 = t213 * t269 + t294;
t114 = qJD(1) * t246 - t293;
t111 = t127 + t202;
t110 = t231 + t284;
t101 = qJD(2) * t246 - t294;
t99 = Ifges(7,5) * t105;
t84 = t116 + t349;
t81 = qJD(5) * t169 * t217 + t227;
t80 = -t218 * t219 * t273 + qJD(5) * t145 + t216 * t256;
t71 = pkin(5) * t168 - qJ(6) * t169 + t155;
t70 = qJD(5) * t122 + t216 * t278 - t218 * t279;
t69 = qJD(3) * t168 + qJD(5) * t232;
t54 = Ifges(7,4) * t56;
t53 = Ifges(6,5) * t56;
t52 = Ifges(6,6) * t57;
t51 = Ifges(7,6) * t57;
t45 = pkin(5) * t144 - qJ(6) * t145 + t123;
t38 = -mrSges(7,1) * t233 - mrSges(7,3) * t105;
t37 = pkin(5) * t105 - qJ(6) * t233;
t31 = Ifges(6,2) * t233 + t201 * Ifges(6,6) + t307;
t28 = t201 * Ifges(7,6) - Ifges(7,3) * t233 + t99;
t27 = -pkin(5) * t219 - t35;
t26 = qJ(6) * t219 + t315;
t23 = pkin(5) * t283 - t24;
t22 = -qJ(6) * t283 + t25;
t19 = t57 * mrSges(6,1) + t56 * mrSges(6,2);
t18 = t57 * mrSges(7,1) - t56 * mrSges(7,3);
t17 = -pkin(5) * t233 - qJ(6) * t105 + t59;
t10 = t80 * pkin(5) - t81 * qJ(6) - t145 * qJD(6) + t87;
t7 = pkin(5) * t281 - t9;
t6 = -qJ(6) * t281 + qJD(6) * t219 + t8;
t5 = -t57 * pkin(5) + t56 * qJ(6) + t105 * qJD(6) + t62;
t11 = [(t144 * t354 + t145 * t355) * t343 + m(5) * (t101 * t88 + t110 * t82 + t124 * t60 + t125 * t68 + t142 * t77 + t84 * t90) + m(7) * (t1 * t26 + t10 * t17 + t15 * t7 + t16 * t6 + t2 * t27 - t45 * t5) + (-t60 * mrSges(5,3) + t86 * mrSges(4,2) + t68 * mrSges(5,1) - t85 * mrSges(4,1) + t54 / 0.2e1 + t53 / 0.2e1 - t52 / 0.2e1 + t51 / 0.2e1 + t264 * t57 - t267 * t56 + ((-0.2e1 * t337 + (-0.3e1 / 0.2e1 * t308 - 0.3e1 / 0.2e1 * t301 - 0.3e1 / 0.2e1 * t306 + 0.3e1 / 0.2e1 * t302 + 0.3e1 / 0.2e1 * Ifges(3,4)) * t219 + (0.3e1 / 0.2e1 * Ifges(3,1) + m(4) * pkin(7) ^ 2 - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + (pkin(7) * mrSges(4,2) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t214) * t214 + (pkin(7) * mrSges(4,1) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t213 + t367 * t214) * t213 - t365) * t217) * qJD(1) + t356 + t347) * qJD(2) + t226) * t219 + (pkin(7) * t140 + t223 * t214 + t224 * t213 + (-pkin(7) * t189 + ((-0.3e1 / 0.2e1 * Ifges(3,4) + t222) * t217 - 0.2e1 * t338 + t267 * t145 - t264 * t144) * qJD(1) + t255 + t346) * qJD(2)) * t217 + (t352 * t328 + t354 * t333 - t21 * mrSges(6,3) - t16 * mrSges(7,2) + t28 / 0.2e1 - t31 / 0.2e1 + t17 * mrSges(7,1) + t59 * mrSges(6,1) + Ifges(7,3) * t335 - Ifges(6,2) * t336) * t80 + (-t1 * t144 + t145 * t2) * mrSges(7,2) + (t351 / 0.2e1 + t353 * t328 + t355 * t333 + mrSges(6,2) * t59 + mrSges(7,2) * t15 - mrSges(6,3) * t20 - mrSges(7,3) * t17 + Ifges(6,4) * t336 + Ifges(7,5) * t335) * t81 + (-t144 * t3 - t145 * t4) * mrSges(6,3) + m(6) * (-t123 * t62 + t20 * t9 + t21 * t8 + t3 * t315 + t35 * t4 + t59 * t87) + t315 * t43 + t145 * t372 + t10 * t38 + t26 * t40 + t35 * t41 + t27 * t42 + t45 * t18 + t8 * t72 + t9 * t73 + t7 * t74 + t6 * t75 + t87 * t39 + (Ifges(7,5) * t145 + Ifges(7,3) * t144) * t341 + (Ifges(6,4) * t145 - Ifges(6,2) * t144) * t342 + t144 * t344 + t144 * t345 + t110 * t112 + t123 * t19 + t84 * t132 + t116 * t133 + t115 * t134 + t101 * t135 + t142 * t139 - t5 * (mrSges(7,1) * t144 - mrSges(7,3) * t145) - t62 * (mrSges(6,1) * t144 + mrSges(6,2) * t145) + t131 * t146 + t130 * t147 + t125 * t148 + t124 * t149 + m(4) * (t107 * t115 + t108 * t116 + t85 * t130 + t86 * t131); (t136 * t352 + t137 * t353) * t329 + (-t152 * t353 + t153 * t352) * t328 + (t136 * t354 + t137 * t355) * t334 + (t168 * t354 + t169 * t355) * t343 + (-t152 * t355 + t153 * t354) * t333 - m(4) * (t107 * t126 + t108 * t127) + (t232 * t4 + t122 * t3 - t155 * t62 + t370 * t59 + (-t25 + t69) * t21 + (-t24 - t70) * t20) * m(6) + (t31 - t28) * (t136 / 0.2e1 - t153 / 0.2e1) - t318 * t232 + (((t189 + t297) * pkin(7) + (t338 + t311 / 0.2e1) * qJD(1) + t255 + (t352 * t168 + t353 * t169) * t369 + ((Ifges(4,6) - Ifges(5,6)) * t214 + t366 * t213) * t368 - t346) * t217 + ((t219 * t222 + t337) * qJD(1) + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t283 - t204 / 0.2e1 + ((Ifges(4,2) * t214 + t310) * t327 + (-Ifges(5,3) * t214 + t305) * t326 + Ifges(3,5) / 0.2e1 + (-m(4) * pkin(2) - mrSges(4,1) * t214 + mrSges(4,2) * t213 - mrSges(3,1)) * pkin(7) + (t213 * t371 - t304 + t309) * t325) * qJD(2) - t347) * t219) * qJD(1) + (-Ifges(6,4) * t152 + Ifges(7,5) * t137 - Ifges(6,2) * t153 + Ifges(7,3) * t136) * t336 + (-t111 * t90 - t114 * t88 - t129 * t82 + t175 * t77 + (qJ(3) * t60 + qJD(3) * t90) * t214 + (qJ(3) * t68 + qJD(3) * t88 - qJD(4) * t82) * t213) * m(5) + (Ifges(6,4) * t137 - Ifges(7,5) * t152 - Ifges(6,2) * t136 + Ifges(7,3) * t153) * t335 + (t1 * t122 - t232 * t2 - t5 * t71 + t363 * t17 + (-t22 + t69) * t16 + (-t23 + t70) * t15) * m(7) + t363 * t38 + t351 * (-t137 / 0.2e1 - t152 / 0.2e1) + (t300 * qJD(4) - t289 * qJD(3) + (-t147 + t148) * qJ(3) + t223) * t213 + (-t168 * t3 - t169 * t4 + t20 * t287 + t21 * t288) * mrSges(6,3) + (-t1 * t168 - t15 * t287 + t16 * t288 + t169 * t2) * mrSges(7,2) + (-mrSges(7,1) * t288 + mrSges(7,3) * t287) * t17 + (-mrSges(6,1) * t288 - mrSges(6,2) * t287) * t59 + (t290 * qJD(3) + (t146 + t149) * qJ(3) - t224) * t214 + m(4) * (-t107 * t279 + t108 * t278 + (-t213 * t85 + t214 * t86) * qJ(3)) + t169 * t372 - t316 * t70 + t317 * t69 + t319 * t122 + t71 * t18 - t25 * t72 - t24 * t73 - t23 * t74 - t22 * t75 + (Ifges(7,5) * t169 + Ifges(7,3) * t168) * t341 + (Ifges(6,4) * t169 - Ifges(6,2) * t168) * t342 + t168 * t344 + t168 * t345 - t109 * t39 - t129 * t112 - t111 * t132 - t127 * t133 - t126 * t134 - t114 * t135 - pkin(2) * t140 + t155 * t19 - t5 * (mrSges(7,1) * t168 - mrSges(7,3) * t169) - t62 * (mrSges(6,1) * t168 + mrSges(6,2) * t169) + t175 * t139; t183 + (-mrSges(6,1) - mrSges(7,1)) * t57 + (-mrSges(6,2) + mrSges(7,3)) * t56 + t289 * t166 + t290 * t165 + t317 * t233 - t316 * t105 + (m(4) * pkin(7) - t314) * t251 - m(4) * (-t107 * t166 - t108 * t165) + t140 + (t105 * t15 + t16 * t233 + t5) * m(7) + (-t105 * t20 + t21 * t233 + t62) * m(6) + (t165 * t90 - t166 * t88 + t77) * m(5); t185 - t318 * t218 + t319 * t216 + (-t38 - t300) * t166 + t228 * qJD(5) + (-mrSges(5,1) * t281 + (t132 + t228) * t219) * qJD(1) + (t1 * t216 - t17 * t166 - t2 * t218 + t201 * (t15 * t216 + t16 * t218)) * m(7) + (-t59 * t166 + t216 * t3 + t218 * t4 - t201 * (t20 * t216 - t21 * t218)) * m(6) + (t82 * t166 + t282 * t90 + t68) * m(5); (t105 * t16 - t15 * t233) * mrSges(7,2) + t54 + t53 - t52 + t51 + t226 + (t312 + t316) * t21 + (t313 - t317) * t20 - t365 * t252 - t37 * t38 + qJ(6) * t40 - pkin(5) * t42 + qJD(6) * t75 + t31 * t333 + (Ifges(7,3) * t105 + t303) * t336 - t17 * (mrSges(7,1) * t105 - mrSges(7,3) * t233) - t59 * (mrSges(6,1) * t105 + mrSges(6,2) * t233) + (t105 * t352 + t233 * t353) * t329 + (-pkin(5) * t2 + qJ(6) * t1 - t15 * t21 + t16 * t350 - t17 * t37) * m(7) + (-Ifges(6,2) * t105 + t100 + t351) * t335 + (t233 * t355 + t28 - t307 + t99) * t334; t105 * t38 - t201 * t75 + 0.2e1 * (t2 / 0.2e1 + t17 * t333 + t16 * t329) * m(7) + t42;];
tauc  = t11(:);
