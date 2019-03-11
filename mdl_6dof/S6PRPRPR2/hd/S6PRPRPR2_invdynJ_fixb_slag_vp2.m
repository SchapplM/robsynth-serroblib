% Calculate vector of inverse dynamics joint torques for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:29:35
% EndTime: 2019-03-08 19:30:03
% DurationCPUTime: 15.70s
% Computational Cost: add. (5911->601), mult. (13749->846), div. (0->0), fcn. (10896->16), ass. (0->279)
t217 = sin(pkin(6));
t215 = sin(pkin(11));
t219 = cos(pkin(11));
t225 = sin(qJ(2));
t228 = cos(qJ(2));
t255 = t215 * t228 + t219 * t225;
t150 = t255 * t217;
t132 = qJD(1) * t150;
t224 = sin(qJ(4));
t227 = cos(qJ(4));
t259 = pkin(4) * t224 - qJ(5) * t227;
t300 = qJD(5) * t224;
t409 = -qJD(4) * t259 + t132 + t300;
t268 = mrSges(5,1) * t227 - mrSges(5,2) * t224;
t218 = cos(pkin(12));
t207 = pkin(5) * t218 + pkin(4);
t213 = pkin(12) + qJ(6);
t210 = sin(t213);
t211 = cos(t213);
t214 = sin(pkin(12));
t267 = -mrSges(6,1) * t218 + mrSges(6,2) * t214;
t381 = m(6) * pkin(4) + m(7) * t207 + mrSges(7,1) * t211 - mrSges(7,2) * t210 - t267;
t339 = pkin(9) + qJ(5);
t399 = -m(6) * qJ(5) - m(7) * t339 - mrSges(6,3) - mrSges(7,3);
t408 = -t224 * t399 + t227 * t381 + t268;
t315 = t217 * t225;
t288 = qJD(1) * t315;
t188 = t215 * t288;
t304 = qJD(1) * t228;
t287 = t217 * t304;
t135 = t219 * t287 - t188;
t206 = pkin(2) * t215 + pkin(8);
t302 = qJD(4) * t224;
t285 = t206 * t302;
t319 = t214 * t227;
t388 = t135 * t319 + t214 * t285 - t409 * t218;
t311 = t218 * t227;
t407 = t135 * t311 + t409 * t214;
t405 = m(5) + m(6);
t253 = pkin(5) * t224 - pkin(9) * t311;
t404 = qJD(4) * t253 + t388;
t312 = t218 * t224;
t403 = -(-pkin(9) * t319 - t206 * t312) * qJD(4) + t407;
t223 = sin(qJ(6));
t226 = cos(qJ(6));
t177 = t214 * t226 + t218 * t223;
t246 = t177 * t227;
t138 = qJD(2) * t246;
t163 = t177 * qJD(6);
t401 = t138 - t163;
t256 = t214 * t223 - t218 * t226;
t245 = t256 * t227;
t139 = qJD(2) * t245;
t162 = t256 * qJD(6);
t400 = t139 - t162;
t299 = t227 * qJD(2);
t203 = qJD(6) - t299;
t303 = qJD(2) * t224;
t172 = qJD(4) * t218 - t214 * t303;
t173 = qJD(4) * t214 + t218 * t303;
t273 = t226 * t172 - t173 * t223;
t113 = t172 * t223 + t173 * t226;
t334 = Ifges(7,4) * t113;
t47 = Ifges(7,2) * t273 + Ifges(7,6) * t203 + t334;
t358 = t47 / 0.2e1;
t107 = Ifges(7,4) * t273;
t48 = Ifges(7,1) * t113 + Ifges(7,5) * t203 + t107;
t357 = t48 / 0.2e1;
t266 = t214 * mrSges(6,1) + t218 * mrSges(6,2);
t343 = pkin(5) * t214;
t398 = -m(7) * (pkin(8) + t343) - t210 * mrSges(7,1) - t211 * mrSges(7,2) - mrSges(5,3) - t266 - t405 * pkin(8);
t216 = sin(pkin(10));
t220 = cos(pkin(10));
t221 = cos(pkin(6));
t308 = t221 * t228;
t397 = -t216 * t225 + t220 * t308;
t298 = qJD(2) * qJD(4);
t185 = qJDD(2) * t224 + t227 * t298;
t140 = qJDD(4) * t218 - t185 * t214;
t141 = qJDD(4) * t214 + t185 * t218;
t41 = qJD(6) * t273 + t140 * t223 + t141 * t226;
t360 = t41 / 0.2e1;
t42 = -qJD(6) * t113 + t140 * t226 - t141 * t223;
t359 = t42 / 0.2e1;
t351 = t140 / 0.2e1;
t350 = t141 / 0.2e1;
t184 = -t227 * qJDD(2) + t224 * t298;
t174 = qJDD(6) + t184;
t349 = t174 / 0.2e1;
t396 = -t184 / 0.2e1;
t348 = t184 / 0.2e1;
t395 = t185 / 0.2e1;
t394 = qJD(4) / 0.2e1;
t252 = -pkin(4) * t227 - qJ(5) * t224 - pkin(3);
t344 = pkin(2) * t219;
t170 = t252 - t344;
t117 = t214 * t170 + t206 * t311;
t320 = t214 * t224;
t101 = -pkin(9) * t320 + t117;
t154 = t218 * t170;
t89 = -pkin(9) * t312 + t154 + (-t206 * t214 - pkin(5)) * t227;
t34 = -t101 * t223 + t226 * t89;
t393 = qJD(6) * t34 + t223 * t404 - t226 * t403;
t35 = t101 * t226 + t223 * t89;
t392 = -qJD(6) * t35 + t223 * t403 + t226 * t404;
t189 = t339 * t214;
t190 = t339 * t218;
t123 = -t189 * t226 - t190 * t223;
t181 = t259 * qJD(2);
t186 = qJD(2) * pkin(2) + t287;
t127 = t215 * t186 + t219 * t288;
t122 = qJD(2) * pkin(8) + t127;
t115 = t224 * t122;
t201 = qJD(1) * t221 + qJD(3);
t86 = t201 * t227 - t115;
t52 = t218 * t181 - t214 * t86;
t49 = qJD(2) * t253 + t52;
t286 = t214 * t299;
t53 = t214 * t181 + t218 * t86;
t50 = -pkin(9) * t286 + t53;
t391 = -qJD(5) * t256 + qJD(6) * t123 - t223 * t49 - t226 * t50;
t124 = -t189 * t223 + t190 * t226;
t390 = -qJD(5) * t177 - qJD(6) * t124 + t223 * t50 - t226 * t49;
t389 = t173 * Ifges(6,5) + t113 * Ifges(7,5) + t172 * Ifges(6,6) + Ifges(7,6) * t273 - Ifges(6,3) * t299 + t203 * Ifges(7,3);
t387 = -t218 * t285 - t407;
t77 = -t140 * mrSges(6,1) + t141 * mrSges(6,2);
t386 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t185 + t77;
t292 = mrSges(5,3) * t303;
t385 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t172 + mrSges(6,2) * t173 + t292;
t335 = Ifges(6,4) * t218;
t262 = -Ifges(6,2) * t214 + t335;
t336 = Ifges(6,4) * t214;
t264 = Ifges(6,1) * t218 - t336;
t382 = t172 * (Ifges(6,6) * t224 + t227 * t262) + t173 * (Ifges(6,5) * t224 + t227 * t264);
t307 = t228 * t219;
t149 = t215 * t315 - t217 * t307;
t380 = (mrSges(3,1) * t228 - mrSges(3,2) * t225) * t217 - t149 * mrSges(4,1) - t150 * mrSges(4,2);
t136 = mrSges(6,2) * t299 + mrSges(6,3) * t172;
t137 = -mrSges(6,1) * t299 - mrSges(6,3) * t173;
t378 = t136 * t218 - t137 * t214;
t102 = -mrSges(6,2) * t184 + mrSges(6,3) * t140;
t103 = mrSges(6,1) * t184 - mrSges(6,3) * t141;
t377 = t102 * t218 - t103 * t214;
t200 = qJDD(1) * t221 + qJDD(3);
t301 = qJD(4) * t227;
t272 = qJD(2) * t288;
t313 = t217 * t228;
t155 = qJDD(1) * t313 - t272;
t329 = qJDD(2) * pkin(2);
t144 = t155 + t329;
t281 = qJD(2) * t304;
t156 = (qJDD(1) * t225 + t281) * t217;
t79 = t215 * t144 + t219 * t156;
t72 = qJDD(2) * pkin(8) + t79;
t294 = t224 * t200 + t201 * t301 + t227 * t72;
t28 = -t122 * t302 + t294;
t29 = -t122 * t301 + t200 * t227 - t201 * t302 - t224 * t72;
t376 = -t224 * t29 + t227 * t28;
t22 = qJDD(4) * qJ(5) + (qJD(5) - t115) * qJD(4) + t294;
t78 = t144 * t219 - t215 * t156;
t71 = -qJDD(2) * pkin(3) - t78;
t45 = pkin(4) * t184 - qJ(5) * t185 - qJD(2) * t300 + t71;
t10 = t214 * t45 + t218 * t22;
t9 = -t214 * t22 + t218 * t45;
t270 = t10 * t218 - t214 * t9;
t375 = -m(7) - t405;
t374 = mrSges(5,1) + t381;
t373 = mrSges(5,2) + t399;
t209 = Ifges(5,4) * t299;
t372 = t218 * (t173 * Ifges(6,1) + t172 * Ifges(6,4) - Ifges(6,5) * t299) + Ifges(5,1) * t303 + Ifges(5,5) * qJD(4) + t209;
t176 = t215 * t225 - t307;
t309 = t221 * t225;
t305 = -t215 * t308 - t219 * t309;
t100 = -t220 * t176 + t216 * t305;
t371 = t216 * t176 + t220 * t305;
t370 = -mrSges(4,1) - t408;
t126 = t186 * t219 - t188;
t121 = -qJD(2) * pkin(3) - t126;
t322 = t201 * t224;
t87 = t122 * t227 + t322;
t76 = qJD(4) * qJ(5) + t87;
t94 = qJD(2) * t252 - t126;
t32 = -t214 * t76 + t218 * t94;
t33 = t214 * t94 + t218 * t76;
t74 = -qJD(4) * pkin(4) + qJD(5) - t86;
t368 = -t227 * t74 * t266 - t33 * (-mrSges(6,2) * t224 - mrSges(6,3) * t319) - t32 * (mrSges(6,1) * t224 - mrSges(6,3) * t311) - t121 * (mrSges(5,1) * t224 + mrSges(5,2) * t227);
t367 = -m(6) * t74 - t385;
t5 = pkin(5) * t184 - pkin(9) * t141 + t9;
t23 = -pkin(5) * t299 - pkin(9) * t173 + t32;
t27 = pkin(9) * t172 + t33;
t6 = -t223 * t27 + t226 * t23;
t8 = pkin(9) * t140 + t10;
t1 = qJD(6) * t6 + t223 * t5 + t226 * t8;
t7 = t223 * t23 + t226 * t27;
t2 = -qJD(6) * t7 - t223 * t8 + t226 * t5;
t366 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t365 = -mrSges(4,2) - t398;
t338 = Ifges(5,4) * t224;
t263 = t227 * Ifges(5,2) + t338;
t364 = t6 * mrSges(7,1) - Ifges(5,6) * qJD(4) / 0.2e1 - qJD(2) * t263 / 0.2e1 - t7 * mrSges(7,2);
t363 = qJD(2) ^ 2;
t362 = Ifges(7,4) * t360 + Ifges(7,2) * t359 + Ifges(7,6) * t349;
t361 = Ifges(7,1) * t360 + Ifges(7,4) * t359 + Ifges(7,5) * t349;
t356 = Ifges(6,1) * t350 + Ifges(6,4) * t351 + Ifges(6,5) * t348;
t355 = -t273 / 0.2e1;
t354 = t273 / 0.2e1;
t353 = -t113 / 0.2e1;
t352 = t113 / 0.2e1;
t347 = -t203 / 0.2e1;
t346 = t203 / 0.2e1;
t337 = Ifges(5,4) * t227;
t26 = -qJDD(4) * pkin(4) + qJDD(5) - t29;
t332 = t224 * t26;
t316 = t217 * t224;
t314 = t217 * t227;
t202 = pkin(2) * t313;
t15 = -t42 * mrSges(7,1) + t41 * mrSges(7,2);
t297 = t15 + t386;
t296 = Ifges(7,5) * t41 + Ifges(7,6) * t42 + Ifges(7,3) * t174;
t51 = -mrSges(7,1) * t273 + mrSges(7,2) * t113;
t295 = t51 + t385;
t293 = m(4) - t375;
t291 = mrSges(5,3) * t299;
t276 = t206 + t343;
t275 = -t298 / 0.2e1;
t274 = t298 / 0.2e1;
t271 = t397 * pkin(2);
t261 = Ifges(5,5) * t227 - Ifges(5,6) * t224;
t260 = Ifges(6,5) * t218 - Ifges(6,6) * t214;
t258 = -t214 * t32 + t218 * t33;
t120 = t150 * t227 + t221 * t224;
t59 = -t120 * t214 + t149 * t218;
t60 = t120 * t218 + t149 * t214;
t20 = -t223 * t60 + t226 * t59;
t21 = t223 * t59 + t226 * t60;
t257 = -t224 * t86 + t227 * t87;
t119 = t150 * t224 - t221 * t227;
t251 = -t216 * t308 - t220 * t225;
t249 = t224 * (Ifges(5,1) * t227 - t338);
t63 = t220 * t314 - t224 * t371;
t65 = t100 * t224 - t216 * t314;
t247 = -g(1) * t65 - g(2) * t63 - g(3) * t119;
t244 = t176 * t221;
t237 = t251 * pkin(2);
t232 = t227 * (Ifges(6,3) * t224 + t227 * t260);
t208 = -pkin(3) - t344;
t193 = -qJD(4) * mrSges(5,2) + t291;
t180 = t268 * qJD(2);
t159 = t276 * t224;
t157 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t184;
t152 = t256 * t224;
t151 = t177 * t224;
t148 = t276 * t301;
t134 = t176 * t217 * qJD(2);
t133 = qJD(2) * t150;
t125 = mrSges(5,1) * t184 + mrSges(5,2) * t185;
t116 = -t206 * t319 + t154;
t105 = t173 * Ifges(6,4) + t172 * Ifges(6,2) - Ifges(6,6) * t299;
t99 = t216 * t244 - t220 * t255;
t96 = -t216 * t255 - t220 * t244;
t93 = -qJD(4) * t246 + t162 * t224;
t92 = -qJD(4) * t245 - t163 * t224;
t83 = mrSges(7,1) * t203 - mrSges(7,3) * t113;
t82 = -mrSges(7,2) * t203 + mrSges(7,3) * t273;
t69 = t322 + (qJD(2) * t343 + t122) * t227;
t66 = t100 * t227 + t216 * t316;
t64 = -t220 * t316 - t227 * t371;
t58 = -qJD(4) * t119 - t134 * t227;
t57 = qJD(4) * t120 - t134 * t224;
t55 = t141 * Ifges(6,4) + t140 * Ifges(6,2) + t184 * Ifges(6,6);
t54 = -pkin(5) * t172 + t74;
t40 = t133 * t214 + t218 * t58;
t39 = t133 * t218 - t214 * t58;
t31 = -mrSges(7,2) * t174 + mrSges(7,3) * t42;
t30 = mrSges(7,1) * t174 - mrSges(7,3) * t41;
t18 = -pkin(5) * t140 + t26;
t4 = -qJD(6) * t21 - t223 * t40 + t226 * t39;
t3 = qJD(6) * t20 + t223 * t39 + t226 * t40;
t11 = [m(2) * qJDD(1) + t60 * t102 + t59 * t103 + t120 * t157 + t149 * t125 - t133 * t180 + t40 * t136 + t39 * t137 + t58 * t193 + t20 * t30 + t21 * t31 + t3 * t82 + t4 * t83 + (-mrSges(3,1) * t225 - mrSges(3,2) * t228) * t363 * t217 + (-mrSges(4,1) * t133 + mrSges(4,2) * t134) * qJD(2) + t295 * t57 + t297 * t119 + t380 * qJDD(2) + (-m(2) - m(3) - t293) * g(3) + m(7) * (t1 * t21 + t119 * t18 + t2 * t20 + t3 * t7 + t4 * t6 + t54 * t57) + m(6) * (t10 * t60 + t119 * t26 + t32 * t39 + t33 * t40 + t57 * t74 + t59 * t9) + m(5) * (-t119 * t29 + t120 * t28 + t121 * t133 + t149 * t71 - t57 * t86 + t58 * t87) + m(4) * (-t126 * t133 - t127 * t134 - t149 * t78 + t150 * t79 + t200 * t221) + m(3) * (qJDD(1) * t221 ^ 2 + (t155 * t228 + t156 * t225) * t217); (Ifges(7,4) * t92 + Ifges(7,2) * t93) * t354 + (qJD(2) * t132 + t219 * t329 + t78) * mrSges(4,1) - (Ifges(6,5) * t141 + Ifges(6,6) * t140 + Ifges(6,3) * t184 + t296) * t227 / 0.2e1 + (t272 + t155) * mrSges(3,1) + (t217 * t281 - t156) * mrSges(3,2) + (qJD(2) * t135 - t215 * t329 - t79) * mrSges(4,2) + (-m(4) * t237 - t251 * mrSges(3,1) - (t216 * t309 - t220 * t228) * mrSges(3,2) + t375 * (t99 * pkin(3) + t237) + t370 * t99 - t365 * t100) * g(1) + ((t215 * t79 + t219 * t78) * pkin(2) + t126 * t132 - t127 * t135) * m(4) + (Ifges(7,1) * t92 + Ifges(7,4) * t93) * t352 + (-t121 * t132 - t135 * t257 + t208 * t71) * m(5) + (t10 * t117 + t116 * t9 + t388 * t32 + t387 * t33) * m(6) + (-t9 * mrSges(6,1) + Ifges(5,4) * t395 + Ifges(5,2) * t396 + t10 * mrSges(6,2) + (-Ifges(5,2) * t224 + t337) * t274 - Ifges(6,3) * t348 - Ifges(7,3) * t349 - Ifges(6,5) * t350 - Ifges(6,6) * t351 - Ifges(7,6) * t359 - Ifges(7,5) * t360 - t135 * t193 - t366) * t227 + (t261 * t394 - t368) * qJD(4) - t55 * t320 / 0.2e1 + t382 * t394 + t337 * t395 + t263 * t396 + (-t397 * mrSges(3,1) - (-t216 * t228 - t220 * t309) * mrSges(3,2) - m(4) * t271 + t375 * (t96 * pkin(3) + t271) + t370 * t96 + t365 * t371) * g(2) + (-t10 * t320 - t312 * t9) * mrSges(6,3) + (-Ifges(7,4) * t152 - Ifges(7,2) * t151) * t359 + (-t1 * t151 + t152 * t2 - t6 * t92 + t7 * t93) * mrSges(7,3) + (-Ifges(7,1) * t152 - Ifges(7,4) * t151) * t360 + t18 * (mrSges(7,1) * t151 - mrSges(7,2) * t152) + (-Ifges(7,5) * t152 - Ifges(7,6) * t151) * t349 + (-m(4) * t202 + t375 * (-t149 * pkin(3) + t202) + t398 * t150 + t408 * t149 - t380) * g(3) + t376 * mrSges(5,3) - t71 * t268 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (Ifges(7,5) * t352 + Ifges(7,6) * t354 + Ifges(7,3) * t346 + t364 - t87 * mrSges(5,3) + t389 / 0.2e1) * t302 + qJDD(4) * (Ifges(5,5) * t224 + Ifges(5,6) * t227) + t208 * t125 + (t372 / 0.2e1 - t214 * t105 / 0.2e1 - t86 * mrSges(5,3)) * t301 + t132 * t180 + t159 * t15 + t148 * t51 + t116 * t103 + t117 * t102 + t54 * (-mrSges(7,1) * t93 + mrSges(7,2) * t92) + (t385 * t301 + t386 * t224 + m(6) * (t301 * t74 + t332) + t157 * t227 + (t376 + (-t87 * t224 - t86 * t227) * qJD(4)) * m(5)) * t206 + (Ifges(7,5) * t92 + Ifges(7,6) * t93) * t346 + t387 * t136 + t388 * t137 + t392 * t83 + t393 * t82 + (t1 * t35 + t148 * t54 + t159 * t18 + t2 * t34 + t392 * t6 + t393 * t7) * m(7) + (Ifges(5,1) * t185 + Ifges(5,4) * t396 + t260 * t348 + t262 * t351 + t264 * t350) * t224 + t249 * t274 + t232 * t275 + t266 * t332 + t312 * t356 + t92 * t357 + t93 * t358 - t152 * t361 - t151 * t362 + t34 * t30 + t35 * t31 - t193 * t285 + (-m(7) * t54 + t367 - t51) * t135 * t224; -t151 * t30 - t152 * t31 + t92 * t82 + t93 * t83 - t297 * t227 + (t157 + t377) * t224 + ((t193 + t378) * t227 + t295 * t224) * qJD(4) + m(4) * t200 + m(5) * (qJD(4) * t257 + t224 * t28 + t227 * t29) + m(7) * (-t1 * t152 - t151 * t2 - t18 * t227 + t302 * t54 + t6 * t93 + t7 * t92) + m(6) * (-t227 * t26 + t270 * t224 + (t224 * t74 + t227 * t258) * qJD(4)) + (-t221 * g(3) + (-g(1) * t216 + g(2) * t220) * t217) * t293; t400 * t357 + t401 * t358 - (-Ifges(5,2) * t303 + t209 + t372) * t299 / 0.2e1 + (t119 * t374 + t120 * t373) * g(3) + (t373 * t66 + t374 * t65) * g(1) + (t373 * t64 + t374 * t63) * g(2) + (-t193 + t291) * t86 + (Ifges(7,5) * t353 + Ifges(7,6) * t355 + Ifges(7,3) * t347 - t364) * t303 + (t232 / 0.2e1 - t249 / 0.2e1) * t363 + (t292 + t367) * t87 + (-Ifges(7,4) * t162 - Ifges(7,2) * t163) * t354 + (-Ifges(7,5) * t162 - Ifges(7,6) * t163) * t346 + (-Ifges(7,1) * t162 - Ifges(7,4) * t163) * t352 + (-mrSges(7,1) * t401 + mrSges(7,2) * t400) * t54 + (-t1 * t256 - t177 * t2 - t400 * t6 + t401 * t7) * mrSges(7,3) + (t368 - t382 / 0.2e1) * qJD(2) + (-pkin(4) * t26 + qJ(5) * t270 + qJD(5) * t258 - t32 * t52 - t33 * t53) * m(6) + t26 * t267 + Ifges(5,3) * qJDD(4) + (-Ifges(7,5) * t139 - Ifges(7,6) * t138) * t347 + (-Ifges(7,4) * t139 - Ifges(7,2) * t138) * t355 + (-Ifges(7,1) * t139 - Ifges(7,4) * t138) * t353 + t18 * (mrSges(7,1) * t256 + mrSges(7,2) * t177) + (Ifges(7,5) * t177 - Ifges(7,6) * t256) * t349 + (Ifges(7,4) * t177 - Ifges(7,2) * t256) * t359 + (Ifges(7,1) * t177 - Ifges(7,4) * t256) * t360 - t256 * t362 + t218 * t55 / 0.2e1 - t207 * t15 - Ifges(5,6) * t184 + Ifges(5,5) * t185 - t53 * t136 - t52 * t137 + t123 * t30 + t124 * t31 - pkin(4) * t77 - t69 * t51 + t270 * mrSges(6,3) + t377 * qJ(5) + t378 * qJD(5) - t389 * t303 / 0.2e1 + t390 * t83 + t391 * t82 + (t1 * t124 + t123 * t2 - t18 * t207 + t390 * t6 + t391 * t7 - t54 * t69) * m(7) + t261 * t275 + (Ifges(6,5) * t214 + Ifges(6,6) * t218) * t348 + (Ifges(6,1) * t214 + t335) * t350 + (Ifges(6,2) * t218 + t336) * t351 + t214 * t356 + t177 * t361 - t28 * mrSges(5,2) + t29 * mrSges(5,1) + t105 * t286 / 0.2e1; t113 * t83 - t273 * t82 - t172 * t136 + t173 * t137 + t15 + t77 + (t113 * t6 - t273 * t7 + t18 + t247) * m(7) + (-t172 * t33 + t173 * t32 + t247 + t26) * m(6); -t54 * (mrSges(7,1) * t113 + mrSges(7,2) * t273) + (Ifges(7,1) * t273 - t334) * t353 + t47 * t352 + (Ifges(7,5) * t273 - Ifges(7,6) * t113) * t347 - t6 * t82 + t7 * t83 - g(1) * ((-t210 * t66 - t211 * t99) * mrSges(7,1) + (t210 * t99 - t211 * t66) * mrSges(7,2)) - g(2) * ((-t210 * t64 - t211 * t96) * mrSges(7,1) + (t210 * t96 - t211 * t64) * mrSges(7,2)) - g(3) * ((-t120 * t210 + t149 * t211) * mrSges(7,1) + (-t120 * t211 - t149 * t210) * mrSges(7,2)) + (t113 * t7 + t273 * t6) * mrSges(7,3) + t296 + (-Ifges(7,2) * t113 + t107 + t48) * t355 + t366;];
tau  = t11;
