% Calculate vector of inverse dynamics joint torques for
% S5RRRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:48
% EndTime: 2019-12-05 18:38:08
% DurationCPUTime: 10.92s
% Computational Cost: add. (8996->528), mult. (21577->712), div. (0->0), fcn. (15753->16), ass. (0->253)
t266 = cos(qJ(2));
t268 = -pkin(7) - pkin(6);
t226 = t268 * t266;
t214 = qJD(1) * t226;
t265 = cos(qJ(3));
t199 = t265 * t214;
t262 = sin(qJ(2));
t225 = t268 * t262;
t213 = qJD(1) * t225;
t261 = sin(qJ(3));
t161 = -t213 * t261 + t199;
t207 = -t261 * t262 + t265 * t266;
t194 = t207 * qJD(1);
t314 = qJ(4) * t194;
t131 = t161 - t314;
t196 = t261 * t214;
t162 = t265 * t213 + t196;
t208 = t261 * t266 + t262 * t265;
t195 = t208 * qJD(1);
t189 = t195 * qJ(4);
t132 = -t189 + t162;
t258 = sin(pkin(9));
t259 = cos(pkin(9));
t308 = t259 * t261;
t321 = pkin(2) * qJD(3);
t371 = -t259 * t131 + t132 * t258 + (-t258 * t265 - t308) * t321;
t310 = t258 * t261;
t370 = -t258 * t131 - t259 * t132 + (t259 * t265 - t310) * t321;
t287 = t259 * t194 - t195 * t258;
t334 = pkin(8) * t287;
t378 = t334 + t371;
t151 = t194 * t258 + t195 * t259;
t143 = pkin(8) * t151;
t377 = t143 + t370;
t256 = qJD(2) + qJD(3);
t376 = t256 / 0.2e1;
t375 = t287 / 0.2e1;
t239 = pkin(2) * t265 + pkin(3);
t185 = -pkin(2) * t310 + t259 * t239;
t181 = pkin(4) + t185;
t187 = pkin(2) * t308 + t239 * t258;
t260 = sin(qJ(5));
t264 = cos(qJ(5));
t141 = t181 * t260 + t187 * t264;
t374 = -qJD(5) * t141 - t377 * t260 + t378 * t264;
t140 = t181 * t264 - t187 * t260;
t373 = qJD(5) * t140 + t378 * t260 + t377 * t264;
t372 = Ifges(5,4) * t151;
t222 = -mrSges(3,1) * t266 + mrSges(3,2) * t262;
t257 = qJ(2) + qJ(3);
t249 = pkin(9) + t257;
t234 = sin(t249);
t235 = cos(t249);
t250 = sin(t257);
t251 = cos(t257);
t238 = qJ(5) + t249;
t230 = sin(t238);
t231 = cos(t238);
t288 = t231 * mrSges(6,1) - t230 * mrSges(6,2);
t273 = -t251 * mrSges(4,1) - t235 * mrSges(5,1) + t250 * mrSges(4,2) + t234 * mrSges(5,2) - t288;
t369 = t222 + t273;
t300 = qJD(1) * qJD(2);
t217 = qJDD(1) * t266 - t262 * t300;
t263 = sin(qJ(1));
t267 = cos(qJ(1));
t355 = g(1) * t267 + g(2) * t263;
t368 = -t151 * t260 + t264 * t287;
t96 = t151 * t264 + t260 * t287;
t367 = t372 / 0.2e1 + Ifges(5,2) * t375 + Ifges(5,6) * t376;
t254 = qJDD(2) + qJDD(3);
t343 = t254 / 0.2e1;
t340 = t262 / 0.2e1;
t253 = t266 * pkin(2);
t240 = t253 + pkin(1);
t365 = Ifges(5,4) * t287;
t364 = Ifges(4,5) * t208;
t363 = Ifges(4,6) * t207;
t236 = pkin(3) * t259 + pkin(4);
t335 = pkin(3) * t258;
t188 = t236 * t260 + t264 * t335;
t202 = qJD(2) * pkin(2) + t213;
t155 = t265 * t202 + t196;
t126 = t155 - t189;
t156 = t202 * t261 - t199;
t127 = t156 + t314;
t309 = t259 * t127;
t65 = -t126 * t258 - t309;
t51 = t65 - t334;
t116 = t258 * t127;
t66 = t259 * t126 - t116;
t52 = -t143 + t66;
t362 = -t188 * qJD(5) + t260 * t52 - t264 * t51;
t186 = t236 * t264 - t260 * t335;
t361 = t186 * qJD(5) - t260 * t51 - t264 * t52;
t313 = qJDD(1) * pkin(1);
t170 = t261 * t225 - t265 * t226;
t205 = t217 * pkin(6);
t218 = qJDD(1) * t262 + t266 * t300;
t206 = t218 * pkin(6);
t357 = t205 * t266 + t206 * t262;
t354 = 0.2e1 * t343;
t229 = pkin(4) * t235;
t237 = pkin(3) * t251;
t307 = t237 + t253;
t297 = t229 + t307;
t353 = mrSges(2,1) + m(6) * (pkin(1) + t297) + m(5) * (pkin(1) + t307) + m(4) * t240 + m(3) * pkin(1) - t369;
t255 = -qJ(4) + t268;
t352 = mrSges(2,2) + m(6) * (-pkin(8) + t255) - mrSges(6,3) + m(5) * t255 - mrSges(5,3) + m(4) * t268 - mrSges(4,3) - m(3) * pkin(6) - mrSges(3,3);
t351 = m(4) * pkin(2);
t350 = m(5) * pkin(3);
t349 = -t368 / 0.2e1;
t348 = -t96 / 0.2e1;
t347 = t96 / 0.2e1;
t345 = t195 / 0.2e1;
t248 = qJD(5) + t256;
t344 = -t248 / 0.2e1;
t338 = pkin(2) * t262;
t337 = pkin(3) * t195;
t336 = pkin(3) * t250;
t113 = pkin(3) * t256 + t126;
t61 = t259 * t113 - t116;
t48 = pkin(4) * t256 - t143 + t61;
t62 = t258 * t113 + t309;
t50 = t62 + t334;
t12 = -t260 * t50 + t264 * t48;
t331 = t12 * mrSges(6,3);
t13 = t260 * t48 + t264 * t50;
t330 = t13 * mrSges(6,3);
t329 = t61 * mrSges(5,3);
t328 = t62 * mrSges(5,3);
t327 = t96 * Ifges(6,4);
t278 = t207 * qJD(3);
t133 = qJD(1) * t278 + t217 * t261 + t218 * t265;
t166 = qJDD(2) * pkin(2) - pkin(7) * t218 - t206;
t168 = pkin(7) * t217 + t205;
t82 = -qJD(3) * t156 + t265 * t166 - t168 * t261;
t45 = pkin(3) * t254 - qJ(4) * t133 - qJD(4) * t195 + t82;
t279 = t208 * qJD(3);
t134 = -qJD(1) * t279 + t217 * t265 - t218 * t261;
t301 = qJD(3) * t265;
t302 = qJD(3) * t261;
t81 = t261 * t166 + t265 * t168 + t202 * t301 + t214 * t302;
t47 = qJ(4) * t134 + qJD(4) * t194 + t81;
t11 = t258 * t45 + t259 * t47;
t296 = qJD(2) * t268;
t215 = t262 * t296;
t216 = t266 * t296;
t108 = t265 * t215 + t261 * t216 + t225 * t301 + t226 * t302;
t164 = -qJD(2) * t208 - t279;
t74 = qJ(4) * t164 + qJD(4) * t207 + t108;
t109 = -qJD(3) * t170 - t215 * t261 + t265 * t216;
t163 = qJD(2) * t207 + t278;
t75 = -qJ(4) * t163 - qJD(4) * t208 + t109;
t36 = t258 * t75 + t259 * t74;
t326 = mrSges(3,2) * t266;
t325 = mrSges(4,3) * t195;
t324 = Ifges(3,4) * t262;
t323 = Ifges(3,4) * t266;
t322 = Ifges(4,4) * t195;
t320 = t155 * mrSges(4,3);
t319 = t156 * mrSges(4,3);
t315 = t266 * Ifges(3,2);
t169 = t265 * t225 + t226 * t261;
t146 = -qJ(4) * t208 + t169;
t147 = qJ(4) * t207 + t170;
t86 = t258 * t146 + t259 * t147;
t306 = qJD(1) * t262;
t305 = qJD(1) * t266;
t304 = qJD(2) * t262;
t303 = qJD(2) * t266;
t244 = pkin(2) * t304;
t78 = -t133 * t258 + t134 * t259;
t79 = t133 * t259 + t134 * t258;
t292 = -t78 * mrSges(5,1) + t79 * mrSges(5,2);
t21 = qJD(5) * t368 + t260 * t78 + t264 * t79;
t22 = -qJD(5) * t96 - t260 * t79 + t264 * t78;
t291 = -t22 * mrSges(6,1) + t21 * mrSges(6,2);
t152 = -pkin(3) * t164 + t244;
t10 = -t258 * t47 + t259 * t45;
t35 = -t258 * t74 + t259 * t75;
t85 = t259 * t146 - t147 * t258;
t176 = -pkin(3) * t207 - t240;
t110 = pkin(4) * t151 + t337;
t203 = -pkin(4) * t234 - t336;
t286 = -g(1) * t263 + g(2) * t267;
t224 = t240 * qJD(1);
t285 = mrSges(6,1) * t230 + mrSges(6,2) * t231;
t284 = t315 + t324;
t283 = Ifges(3,5) * t266 - Ifges(3,6) * t262;
t190 = -pkin(2) * t217 - t313;
t158 = t207 * t258 + t208 * t259;
t57 = -pkin(8) * t158 + t85;
t157 = t207 * t259 - t208 * t258;
t58 = pkin(8) * t157 + t86;
t31 = -t260 * t58 + t264 * t57;
t32 = t260 * t57 + t264 * t58;
t99 = t157 * t264 - t158 * t260;
t100 = t157 * t260 + t158 * t264;
t6 = pkin(4) * t254 - pkin(8) * t79 + t10;
t7 = pkin(8) * t78 + t11;
t2 = qJD(5) * t12 + t260 * t6 + t264 * t7;
t247 = qJDD(5) + t254;
t3 = -qJD(5) * t13 - t260 * t7 + t264 * t6;
t282 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t21 + Ifges(6,6) * t22 + Ifges(6,3) * t247;
t281 = pkin(1) * (mrSges(3,1) * t262 + t326);
t280 = t262 * (Ifges(3,1) * t266 - t324);
t101 = -pkin(3) * t134 + qJDD(4) + t190;
t167 = -pkin(3) * t194 + qJD(4) - t224;
t272 = mrSges(5,1) * t234 + mrSges(4,2) * t251 + mrSges(5,2) * t235 + t285;
t104 = -pkin(4) * t287 + t167;
t144 = t194 * Ifges(4,2) + t256 * Ifges(4,6) + t322;
t191 = Ifges(4,4) * t194;
t145 = t195 * Ifges(4,1) + t256 * Ifges(4,5) + t191;
t38 = Ifges(6,2) * t368 + t248 * Ifges(6,6) + t327;
t92 = Ifges(6,4) * t368;
t39 = Ifges(6,1) * t96 + Ifges(6,5) * t248 + t92;
t91 = t151 * Ifges(5,1) + t256 * Ifges(5,5) + t365;
t269 = t282 + t224 * (mrSges(4,1) * t195 + mrSges(4,2) * t194) + t287 * t329 - t167 * (mrSges(5,1) * t151 + mrSges(5,2) * t287) - t151 * (Ifges(5,1) * t287 - t372) / 0.2e1 - t195 * (Ifges(4,1) * t194 - t322) / 0.2e1 + t151 * t328 + t144 * t345 + t194 * t320 + t151 * t367 + t10 * mrSges(5,1) - t11 * mrSges(5,2) + Ifges(5,6) * t78 + Ifges(5,5) * t79 - t81 * mrSges(4,2) + t82 * mrSges(4,1) + Ifges(4,5) * t133 + Ifges(4,6) * t134 + (-t104 * mrSges(6,2) + Ifges(6,5) * t344 + Ifges(6,1) * t348 + Ifges(6,4) * t349 + t331 - t39 / 0.2e1) * t368 - (t104 * mrSges(6,1) - t330 + Ifges(6,6) * t344 + Ifges(6,4) * t348 + Ifges(6,2) * t349 - t38 / 0.2e1) * t96 - (-Ifges(5,2) * t151 + t365 + t91) * t287 / 0.2e1 - (-Ifges(4,2) * t195 + t145 + t191) * t194 / 0.2e1 - (Ifges(4,5) * t194 + Ifges(5,5) * t287 - Ifges(4,6) * t195 - Ifges(5,6) * t151) * t256 / 0.2e1 + (Ifges(4,3) + Ifges(5,3)) * t254;
t243 = pkin(2) * t306;
t242 = Ifges(3,4) * t305;
t221 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t305;
t220 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t306;
t193 = Ifges(3,1) * t306 + Ifges(3,5) * qJD(2) + t242;
t192 = Ifges(3,6) * qJD(2) + qJD(1) * t284;
t173 = mrSges(4,1) * t256 - t325;
t172 = -mrSges(4,2) * t256 + mrSges(4,3) * t194;
t171 = t243 + t337;
t154 = -mrSges(4,1) * t194 + mrSges(4,2) * t195;
t136 = mrSges(5,1) * t256 - mrSges(5,3) * t151;
t135 = -mrSges(5,2) * t256 + mrSges(5,3) * t287;
t120 = -pkin(4) * t157 + t176;
t115 = -mrSges(4,2) * t254 + mrSges(4,3) * t134;
t114 = mrSges(4,1) * t254 - mrSges(4,3) * t133;
t105 = t110 + t243;
t103 = t163 * t259 + t164 * t258;
t102 = -t163 * t258 + t164 * t259;
t98 = -mrSges(5,1) * t287 + mrSges(5,2) * t151;
t84 = mrSges(6,1) * t248 - mrSges(6,3) * t96;
t83 = -mrSges(6,2) * t248 + mrSges(6,3) * t368;
t70 = -pkin(4) * t102 + t152;
t64 = mrSges(5,1) * t254 - mrSges(5,3) * t79;
t63 = -mrSges(5,2) * t254 + mrSges(5,3) * t78;
t44 = -mrSges(6,1) * t368 + mrSges(6,2) * t96;
t42 = -pkin(4) * t78 + t101;
t34 = -qJD(5) * t100 + t102 * t264 - t103 * t260;
t33 = qJD(5) * t99 + t102 * t260 + t103 * t264;
t26 = pkin(8) * t102 + t36;
t25 = -pkin(8) * t103 + t35;
t15 = -mrSges(6,2) * t247 + mrSges(6,3) * t22;
t14 = mrSges(6,1) * t247 - mrSges(6,3) * t21;
t5 = -qJD(5) * t32 + t25 * t264 - t26 * t260;
t4 = qJD(5) * t31 + t25 * t260 + t26 * t264;
t1 = [(-mrSges(4,2) * t240 + Ifges(4,1) * t208 + Ifges(4,4) * t207) * t133 - t163 * t320 + t266 * (Ifges(3,4) * t218 + Ifges(3,2) * t217) / 0.2e1 + (t363 / 0.2e1 + t364 / 0.2e1) * t254 + (t363 + t364) * t343 + (-t101 * mrSges(5,1) + t11 * mrSges(5,3) + Ifges(5,4) * t79 + Ifges(5,2) * t78 + Ifges(5,6) * t354) * t157 + (t101 * mrSges(5,2) - t10 * mrSges(5,3) + Ifges(5,1) * t79 + Ifges(5,4) * t78 + Ifges(5,5) * t354) * t158 + (m(3) * t357 - t220 * t303 - t221 * t304 - t262 * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t218) + t266 * (-qJDD(2) * mrSges(3,2) + mrSges(3,3) * t217)) * pkin(6) + t357 * mrSges(3,3) + (t263 * t353 + t267 * t352) * g(1) + (t263 * t352 - t267 * t353) * g(2) + (t266 * (-Ifges(3,2) * t262 + t323) + t280) * t300 / 0.2e1 + (mrSges(4,1) * t240 + Ifges(4,4) * t208 + Ifges(4,2) * t207) * t134 + (-mrSges(6,1) * t42 + mrSges(6,3) * t2 + Ifges(6,4) * t21 + Ifges(6,2) * t22 + Ifges(6,6) * t247) * t99 - t281 * t300 + (mrSges(6,2) * t42 - mrSges(6,3) * t3 + Ifges(6,1) * t21 + Ifges(6,4) * t22 + Ifges(6,5) * t247) * t100 - t33 * t331 - t222 * t313 + t218 * (t262 * Ifges(3,1) + t323) / 0.2e1 - t103 * t329 + t193 * t303 / 0.2e1 - t192 * t304 / 0.2e1 + (Ifges(4,1) * t163 + Ifges(4,4) * t164) * t345 + (Ifges(6,1) * t33 + Ifges(6,4) * t34) * t347 + t102 * t328 + t34 * t330 + t164 * t319 + t154 * t244 + (t207 * t81 - t208 * t82) * mrSges(4,3) + qJD(2) ^ 2 * t283 / 0.2e1 + t217 * t284 / 0.2e1 + (Ifges(3,1) * t218 + Ifges(3,4) * t217) * t340 + t120 * t291 + t176 * t292 + t102 * t367 + (m(3) * t313 + mrSges(3,1) * t217 - mrSges(3,2) * t218) * pkin(1) + (Ifges(5,4) * t103 + Ifges(5,2) * t102) * t375 + (Ifges(4,5) * t163 + Ifges(5,5) * t103 + Ifges(4,6) * t164 + Ifges(5,6) * t102) * t376 + (0.2e1 * Ifges(3,5) * t340 + Ifges(3,6) * t266) * qJDD(2) + t31 * t14 + t32 * t15 + t34 * t38 / 0.2e1 + t33 * t39 / 0.2e1 + t70 * t44 + t4 * t83 + t5 * t84 + t85 * t64 + t86 * t63 + t103 * t91 / 0.2e1 + t104 * (-mrSges(6,1) * t34 + mrSges(6,2) * t33) + t36 * t135 + t35 * t136 + m(6) * (t104 * t70 + t12 * t5 + t120 * t42 + t13 * t4 + t2 * t32 + t3 * t31) + m(5) * (t10 * t85 + t101 * t176 + t11 * t86 + t152 * t167 + t35 * t61 + t36 * t62) + t151 * (Ifges(5,1) * t103 + Ifges(5,4) * t102) / 0.2e1 + t152 * t98 + t163 * t145 / 0.2e1 + t164 * t144 / 0.2e1 + t167 * (-mrSges(5,1) * t102 + mrSges(5,2) * t103) + t169 * t114 + t170 * t115 + t108 * t172 + t109 * t173 + t194 * (Ifges(4,4) * t163 + Ifges(4,2) * t164) / 0.2e1 + t190 * (-mrSges(4,1) * t207 + mrSges(4,2) * t208) + Ifges(2,3) * qJDD(1) - t224 * (-mrSges(4,1) * t164 + mrSges(4,2) * t163) + t248 * (Ifges(6,5) * t33 + Ifges(6,6) * t34) / 0.2e1 + t368 * (Ifges(6,4) * t33 + Ifges(6,2) * t34) / 0.2e1 + m(4) * (t108 * t156 + t109 * t155 + t169 * t82 + t170 * t81 - t190 * t240 - t224 * t244); t370 * t135 + (-t307 * g(3) + t10 * t185 + t11 * t187 - t167 * t171 + t370 * t62 + t371 * t61) * m(5) + t371 * t136 - m(4) * (t155 * t161 + t156 * t162) + t369 * g(3) + t373 * t83 + (-t297 * g(3) - t104 * t105 + t374 * t12 + t373 * t13 + t140 * t3 + t141 * t2) * m(6) + t374 * t84 + t355 * (-m(5) * (-t336 - t338) - m(6) * (t203 - t338) + mrSges(4,1) * t250 + t326 + (mrSges(3,1) + t351) * t262 + t272) + t269 + (t261 * t81 + t265 * t82 + (-t155 * t261 + t156 * t265) * qJD(3)) * t351 + t195 * t319 + (-g(3) * m(4) * t266 + t265 * t114 + t261 * t115 + (t172 * t265 - t173 * t261) * qJD(3)) * pkin(2) + Ifges(3,3) * qJDD(2) - t105 * t44 + t140 * t14 + t141 * t15 - t171 * t98 - t162 * t172 - t161 * t173 + t185 * t64 + t187 * t63 - t205 * mrSges(3,2) - t206 * mrSges(3,1) + Ifges(3,6) * t217 + Ifges(3,5) * t218 + (t192 * t340 - qJD(2) * t283 / 0.2e1 + (-t280 / 0.2e1 + t281 + t315 * t340) * qJD(1) + (t220 * t266 + t221 * t262) * pkin(6) + (m(4) * t224 - t154) * t338 - (t193 + t242) * t266 / 0.2e1) * qJD(1); (t10 * t259 + t11 * t258) * t350 - m(5) * (t167 * t337 + t61 * t65 + t62 * t66) + t362 * t84 + t273 * g(3) + (t173 + t325) * t156 + t361 * t83 + t269 - t110 * t44 - t66 * t135 - t65 * t136 - t155 * t172 + t186 * t14 + t188 * t15 + (-m(5) * g(3) * t251 - t195 * t98 + t258 * t63 + t259 * t64) * pkin(3) + t355 * (-m(6) * t203 + (mrSges(4,1) + t350) * t250 + t272) + (-t104 * t110 + t186 * t3 + t188 * t2 + (-t229 - t237) * g(3) + t361 * t13 + t362 * t12) * m(6); -t287 * t135 + t151 * t136 - t368 * t83 + t96 * t84 + t291 + t292 + (t12 * t96 - t13 * t368 + t286 + t42) * m(6) + (t151 * t61 - t287 * t62 + t101 + t286) * m(5); -t104 * (mrSges(6,1) * t96 + mrSges(6,2) * t368) + (Ifges(6,1) * t368 - t327) * t348 + t38 * t347 + (Ifges(6,5) * t368 - Ifges(6,6) * t96) * t344 - t12 * t83 + t13 * t84 - g(3) * t288 + (t12 * t368 + t13 * t96) * mrSges(6,3) + t282 + (-Ifges(6,2) * t96 + t39 + t92) * t349 + t355 * t285;];
tau = t1;
