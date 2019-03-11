% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:38:00
% EndTime: 2019-03-09 02:38:28
% DurationCPUTime: 17.47s
% Computational Cost: add. (8070->657), mult. (17459->873), div. (0->0), fcn. (12375->18), ass. (0->292)
t216 = sin(pkin(11));
t219 = cos(pkin(11));
t253 = -t219 * mrSges(6,1) + t216 * mrSges(6,2);
t387 = -m(6) * pkin(4) - mrSges(5,1) + t253;
t338 = m(6) + m(7);
t278 = m(5) + t338;
t386 = mrSges(5,2) - mrSges(7,3);
t224 = sin(qJ(3));
t227 = cos(qJ(3));
t190 = -mrSges(4,1) * t227 + mrSges(4,2) * t224;
t214 = qJ(3) + pkin(10);
t203 = sin(t214);
t206 = cos(t214);
t385 = t386 * t203 + t206 * t387 + t190;
t384 = m(6) * qJ(5) + mrSges(6,3);
t383 = Ifges(5,5) * qJD(3);
t382 = Ifges(5,6) * qJD(3);
t218 = sin(pkin(9));
t194 = pkin(1) * t218 + pkin(7);
t184 = t194 * qJDD(1);
t381 = qJD(2) * qJD(3) + t184;
t215 = qJ(1) + pkin(9);
t204 = sin(t215);
t207 = cos(t215);
t380 = g(1) * t207 + g(2) * t204;
t379 = -m(4) * pkin(2) - mrSges(3,1) + t385;
t217 = sin(pkin(10));
t301 = cos(pkin(10));
t176 = t217 * t227 + t224 * t301;
t163 = t176 * qJD(1);
t144 = qJD(3) * t216 + t163 * t219;
t223 = sin(qJ(6));
t226 = cos(qJ(6));
t260 = t219 * qJD(3) - t163 * t216;
t377 = -t144 * t223 + t226 * t260;
t85 = t144 * t226 + t223 * t260;
t277 = qJD(1) * qJD(3);
t179 = qJDD(1) * t227 - t224 * t277;
t180 = qJDD(1) * t224 + t227 * t277;
t131 = t217 * t179 + t180 * t301;
t110 = qJDD(3) * t219 - t131 * t216;
t111 = qJDD(3) * t216 + t131 * t219;
t27 = qJD(6) * t377 + t110 * t223 + t111 * t226;
t347 = t27 / 0.2e1;
t28 = -qJD(6) * t85 + t110 * t226 - t111 * t223;
t346 = t28 / 0.2e1;
t376 = -m(4) - m(3);
t375 = -m(5) - m(7);
t263 = t301 * t227;
t280 = t224 * qJD(1);
t161 = -qJD(1) * t263 + t217 * t280;
t210 = t227 * qJD(2);
t186 = t194 * qJD(1);
t258 = qJ(4) * qJD(1) + t186;
t147 = -t224 * t258 + t210;
t136 = qJD(3) * pkin(3) + t147;
t283 = qJD(2) * t224;
t148 = t227 * t258 + t283;
t264 = t301 * t148;
t81 = t217 * t136 + t264;
t72 = qJD(3) * qJ(5) + t81;
t220 = cos(pkin(9));
t197 = -pkin(1) * t220 - pkin(2);
t211 = t227 * pkin(3);
t183 = t197 - t211;
t160 = qJD(1) * t183 + qJD(4);
t99 = pkin(4) * t161 - qJ(5) * t163 + t160;
t40 = -t216 * t72 + t219 * t99;
t21 = pkin(5) * t161 - pkin(8) * t144 + t40;
t41 = t216 * t99 + t219 * t72;
t26 = pkin(8) * t260 + t41;
t8 = t21 * t226 - t223 * t26;
t374 = t8 * mrSges(7,1);
t9 = t21 * t223 + t226 * t26;
t373 = t9 * mrSges(7,2);
t337 = t110 / 0.2e1;
t336 = t111 / 0.2e1;
t130 = -t301 * t179 + t180 * t217;
t129 = qJDD(6) + t130;
t335 = t129 / 0.2e1;
t334 = t130 / 0.2e1;
t372 = t40 * mrSges(6,1);
t371 = t41 * mrSges(6,2);
t157 = qJD(6) + t161;
t366 = t260 * Ifges(6,6);
t367 = t144 * Ifges(6,5);
t370 = t85 * Ifges(7,5) + Ifges(7,6) * t377 + t161 * Ifges(6,3) + t157 * Ifges(7,3) + t366 + t367;
t322 = pkin(3) * t217;
t192 = qJ(5) + t322;
t316 = pkin(8) + t192;
t169 = t316 * t216;
t170 = t316 * t219;
t120 = -t169 * t226 - t170 * t223;
t242 = t216 * t223 - t219 * t226;
t293 = t161 * t219;
t272 = pkin(3) * t280;
t109 = pkin(4) * t163 + qJ(5) * t161 + t272;
t133 = t217 * t148;
t89 = t147 * t301 - t133;
t51 = t219 * t109 - t216 * t89;
t29 = pkin(5) * t163 + pkin(8) * t293 + t51;
t294 = t161 * t216;
t52 = t216 * t109 + t219 * t89;
t39 = pkin(8) * t294 + t52;
t369 = -qJD(5) * t242 + qJD(6) * t120 - t223 * t29 - t226 * t39;
t121 = -t169 * t223 + t170 * t226;
t177 = t216 * t226 + t219 * t223;
t368 = -qJD(5) * t177 - qJD(6) * t121 + t223 * t39 - t226 * t29;
t252 = mrSges(6,1) * t216 + mrSges(6,2) * t219;
t80 = t136 * t301 - t133;
t71 = -qJD(3) * pkin(4) + qJD(5) - t80;
t365 = t71 * t252;
t302 = qJDD(3) / 0.2e1;
t125 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t131;
t57 = -t110 * mrSges(6,1) + t111 * mrSges(6,2);
t364 = t57 - t125;
t314 = mrSges(5,3) * t163;
t363 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t260 + mrSges(6,2) * t144 + t314;
t362 = t203 * t380;
t104 = t177 * t161;
t166 = t177 * qJD(6);
t359 = t104 + t166;
t105 = t242 * t161;
t165 = t242 * qJD(6);
t358 = t105 + t165;
t282 = qJD(3) * t224;
t118 = t224 * qJDD(2) - t186 * t282 + t227 * t381;
t154 = t186 * t227 + t283;
t209 = t227 * qJDD(2);
t119 = -qJD(3) * t154 - t184 * t224 + t209;
t357 = t118 * t227 - t119 * t224;
t100 = -mrSges(6,2) * t161 + mrSges(6,3) * t260;
t101 = mrSges(6,1) * t161 - mrSges(6,3) * t144;
t356 = t100 * t219 - t101 * t216;
t63 = -mrSges(6,2) * t130 + mrSges(6,3) * t110;
t64 = mrSges(6,1) * t130 - mrSges(6,3) * t111;
t355 = -t216 * t64 + t219 * t63;
t276 = qJD(1) * qJD(4);
t281 = qJD(3) * t227;
t75 = -t186 * t281 + qJDD(3) * pkin(3) - qJ(4) * t180 + t209 + (-t276 - t381) * t224;
t82 = qJ(4) * t179 + t227 * t276 + t118;
t38 = t217 * t75 + t301 * t82;
t35 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t38;
t185 = t197 * qJDD(1);
t146 = -pkin(3) * t179 + qJDD(4) + t185;
t50 = pkin(4) * t130 - qJ(5) * t131 - qJD(5) * t163 + t146;
t14 = -t216 * t35 + t219 * t50;
t15 = t216 * t50 + t219 * t35;
t245 = -t14 * t216 + t15 * t219;
t354 = 0.2e1 * t302;
t11 = pkin(8) * t110 + t15;
t7 = pkin(5) * t130 - pkin(8) * t111 + t14;
t1 = qJD(6) * t8 + t11 * t226 + t223 * t7;
t2 = -qJD(6) * t9 - t11 * t223 + t226 * t7;
t351 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t221 = -qJ(4) - pkin(7);
t350 = -m(7) * pkin(5) * t216 - m(4) * pkin(7) + m(6) * t221 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - t252;
t349 = Ifges(7,4) * t347 + Ifges(7,2) * t346 + Ifges(7,6) * t335;
t348 = Ifges(7,1) * t347 + Ifges(7,4) * t346 + Ifges(7,5) * t335;
t324 = Ifges(7,4) * t85;
t31 = Ifges(7,2) * t377 + Ifges(7,6) * t157 + t324;
t345 = t31 / 0.2e1;
t79 = Ifges(7,4) * t377;
t32 = Ifges(7,1) * t85 + Ifges(7,5) * t157 + t79;
t344 = t32 / 0.2e1;
t343 = Ifges(6,1) * t336 + Ifges(6,4) * t337 + Ifges(6,5) * t334;
t342 = -t377 / 0.2e1;
t341 = t377 / 0.2e1;
t340 = -t85 / 0.2e1;
t339 = t85 / 0.2e1;
t333 = -t157 / 0.2e1;
t332 = t157 / 0.2e1;
t331 = -t161 / 0.2e1;
t330 = t161 / 0.2e1;
t327 = t163 / 0.2e1;
t325 = t219 / 0.2e1;
t225 = sin(qJ(1));
t323 = pkin(1) * t225;
t317 = t219 * pkin(5);
t228 = cos(qJ(1));
t212 = t228 * pkin(1);
t162 = t176 * qJD(3);
t238 = -t217 * t224 + t263;
t164 = t238 * qJD(3);
t271 = pkin(3) * t282;
t87 = pkin(4) * t162 - qJ(5) * t164 - qJD(5) * t176 + t271;
t285 = qJ(4) + t194;
t259 = qJD(3) * t285;
t149 = qJD(4) * t227 - t224 * t259;
t150 = -qJD(4) * t224 - t227 * t259;
t98 = t149 * t301 + t217 * t150;
t47 = t216 * t87 + t219 * t98;
t315 = mrSges(5,3) * t161;
t313 = Ifges(4,4) * t224;
t312 = Ifges(4,4) * t227;
t311 = Ifges(5,4) * t163;
t310 = Ifges(6,4) * t216;
t309 = Ifges(6,4) * t219;
t306 = t203 * mrSges(6,3);
t299 = qJ(5) * t203;
t292 = t164 * t216;
t291 = t164 * t219;
t290 = t176 * t216;
t289 = t176 * t219;
t222 = -pkin(8) - qJ(5);
t288 = t203 * t222;
t287 = t204 * t206;
t286 = t206 * t207;
t113 = -pkin(4) * t238 - qJ(5) * t176 + t183;
t172 = t285 * t224;
t173 = t285 * t227;
t123 = -t217 * t172 + t173 * t301;
t59 = t216 * t113 + t219 * t123;
t199 = t211 + pkin(2);
t284 = t207 * t199 + t212;
t279 = t227 * qJD(1);
t274 = Ifges(7,5) * t27 + Ifges(7,6) * t28 + Ifges(7,3) * t129;
t42 = -mrSges(7,1) * t377 + mrSges(7,2) * t85;
t273 = t42 + t363;
t270 = mrSges(4,3) * t280;
t269 = mrSges(4,3) * t279;
t268 = -t216 * (t144 * Ifges(6,4) + Ifges(6,2) * t260 + t161 * Ifges(6,6)) / 0.2e1;
t267 = (t144 * Ifges(6,1) + Ifges(6,4) * t260 + t161 * Ifges(6,5)) * t325;
t266 = t301 * pkin(3);
t10 = -t28 * mrSges(7,1) + t27 * mrSges(7,2);
t46 = -t216 * t98 + t219 * t87;
t261 = t130 * mrSges(5,1) + t131 * mrSges(5,2);
t58 = t219 * t113 - t123 * t216;
t88 = t147 * t217 + t264;
t97 = t149 * t217 - t301 * t150;
t122 = t301 * t172 + t173 * t217;
t196 = -t266 - pkin(4);
t37 = -t217 * t82 + t301 * t75;
t255 = mrSges(4,1) * t224 + mrSges(4,2) * t227;
t251 = Ifges(6,1) * t219 - t310;
t250 = t227 * Ifges(4,2) + t313;
t249 = -Ifges(6,2) * t216 + t309;
t248 = Ifges(4,5) * t227 - Ifges(4,6) * t224;
t247 = Ifges(6,5) * t219 - Ifges(6,6) * t216;
t244 = t216 * t40 - t219 * t41;
t45 = -pkin(5) * t238 - pkin(8) * t289 + t58;
t53 = -pkin(8) * t290 + t59;
t16 = -t223 * t53 + t226 * t45;
t17 = t223 * t45 + t226 * t53;
t195 = pkin(4) + t317;
t243 = t195 * t206 - t288;
t151 = -qJD(3) * mrSges(5,2) - t315;
t241 = t151 + t356;
t240 = t197 * qJD(1) * t255;
t239 = t224 * (Ifges(4,1) * t227 - t313);
t213 = pkin(11) + qJ(6);
t202 = sin(t213);
t205 = cos(t213);
t234 = m(7) * t195 + mrSges(7,1) * t205 - mrSges(7,2) * t202;
t36 = -qJDD(3) * pkin(4) + qJDD(5) - t37;
t200 = Ifges(4,4) * t279;
t189 = -qJD(3) * mrSges(4,2) + t269;
t187 = qJD(3) * mrSges(4,1) - t270;
t182 = t196 - t317;
t168 = Ifges(4,1) * t280 + Ifges(4,5) * qJD(3) + t200;
t167 = Ifges(4,6) * qJD(3) + qJD(1) * t250;
t159 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t180;
t158 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t179;
t156 = Ifges(5,4) * t161;
t153 = -t186 * t224 + t210;
t140 = t202 * t204 + t205 * t286;
t139 = -t202 * t286 + t204 * t205;
t138 = t202 * t207 - t205 * t287;
t137 = t202 * t287 + t205 * t207;
t126 = mrSges(5,1) * t161 + mrSges(5,2) * t163;
t124 = -qJDD(3) * mrSges(5,2) - mrSges(5,3) * t130;
t117 = t163 * Ifges(5,1) - t156 + t383;
t116 = -t161 * Ifges(5,2) + t311 + t382;
t115 = t242 * t176;
t114 = t177 * t176;
t95 = pkin(5) * t290 + t122;
t65 = pkin(5) * t292 + t97;
t62 = -pkin(5) * t294 + t88;
t61 = mrSges(7,1) * t157 - mrSges(7,3) * t85;
t60 = -mrSges(7,2) * t157 + mrSges(7,3) * t377;
t56 = -t164 * t177 + t165 * t176;
t55 = -t164 * t242 - t166 * t176;
t54 = -pkin(5) * t260 + t71;
t43 = t111 * Ifges(6,4) + t110 * Ifges(6,2) + t130 * Ifges(6,6);
t34 = -pkin(8) * t292 + t47;
t22 = pkin(5) * t162 - pkin(8) * t291 + t46;
t20 = -t110 * pkin(5) + t36;
t19 = -mrSges(7,2) * t129 + mrSges(7,3) * t28;
t18 = mrSges(7,1) * t129 - mrSges(7,3) * t27;
t4 = -qJD(6) * t17 + t22 * t226 - t223 * t34;
t3 = qJD(6) * t16 + t22 * t223 + t226 * t34;
t5 = [(Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t220 - 0.2e1 * mrSges(3,2) * t218 + m(3) * (t218 ^ 2 + t220 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) - t43 * t290 / 0.2e1 + t180 * t312 / 0.2e1 + t224 * (Ifges(4,4) * t179 + Ifges(4,5) * qJDD(3)) / 0.2e1 + (Ifges(7,4) * t55 + Ifges(7,2) * t56) * t341 + t126 * t271 + (m(4) * t197 + t190) * t185 + t168 * t281 / 0.2e1 - t167 * t282 / 0.2e1 - (t146 * mrSges(5,1) + t14 * mrSges(6,1) - t15 * mrSges(6,2) - t38 * mrSges(5,3) - Ifges(5,4) * t131 + Ifges(6,5) * t336 + Ifges(7,5) * t347 + Ifges(5,2) * t130 - Ifges(5,6) * t354 + Ifges(6,6) * t337 + Ifges(7,6) * t346 + Ifges(6,3) * t334 + Ifges(7,3) * t335 + t351) * t238 - (Ifges(6,5) * t111 + Ifges(6,6) * t110 + Ifges(6,3) * t130 + t274) * t238 / 0.2e1 + (-m(5) * t80 + m(6) * t71 + t363) * t97 + (-m(5) * t37 + m(6) * t36 + t364) * t122 + t180 * t224 * Ifges(4,1) + m(6) * (t14 * t58 + t15 * t59 + t40 * t46 + t41 * t47) + m(5) * (t123 * t38 + t146 * t183 + t160 * t271 + t81 * t98) + (-t14 * t289 - t15 * t290 - t291 * t40 - t292 * t41) * mrSges(6,3) + m(7) * (t1 * t17 + t16 * t2 + t20 * t95 + t3 * t9 + t4 * t8 + t54 * t65) + (m(4) * ((-t153 * t227 - t154 * t224) * qJD(3) + t357) - t187 * t281 - t189 * t282 + t227 * t158 - t224 * t159) * t194 + (-t153 * t281 - t154 * t282 + t357) * mrSges(4,3) + (t146 * mrSges(5,2) - t37 * mrSges(5,3) + Ifges(5,1) * t131 - Ifges(5,4) * t130 + Ifges(5,5) * t354 + t247 * t334 + t249 * t337 + t251 * t336 + t36 * t252) * t176 + (Ifges(7,5) * t55 + Ifges(7,6) * t56) * t332 + (Ifges(4,5) * t224 + Ifges(4,6) * t227) * t302 + t179 * t250 / 0.2e1 + (Ifges(7,1) * t55 + Ifges(7,4) * t56) * t339 + (t240 + t248 * qJD(3) / 0.2e1) * qJD(3) + t183 * t261 + (t227 * (-Ifges(4,2) * t224 + t312) + t239) * t277 / 0.2e1 + (-m(6) * t284 - mrSges(2,1) * t228 - t140 * mrSges(7,1) + mrSges(2,2) * t225 - t139 * mrSges(7,2) + t375 * (-t204 * t221 + t284) + t376 * t212 + t350 * t204 + (-m(7) * t243 - t203 * t384 + t379) * t207) * g(2) + t98 * t151 - t115 * t348 - t114 * t349 + t123 * t124 + t47 * t100 + t46 * t101 + t95 * t10 + t65 * t42 + t3 * t60 + t4 * t61 + t59 * t63 + t58 * t64 + t54 * (-mrSges(7,1) * t56 + mrSges(7,2) * t55) + t16 * t18 + t17 * t19 + t197 * (-mrSges(4,1) * t179 + mrSges(4,2) * t180) + t289 * t343 + t55 * t344 + t56 * t345 + t227 * (Ifges(4,4) * t180 + Ifges(4,2) * t179 + Ifges(4,6) * qJDD(3)) / 0.2e1 + (-Ifges(7,4) * t115 - Ifges(7,2) * t114) * t346 + (-Ifges(7,5) * t115 - Ifges(7,6) * t114) * t335 + (-Ifges(7,1) * t115 - Ifges(7,4) * t114) * t347 + t20 * (mrSges(7,1) * t114 - mrSges(7,2) * t115) + (-t1 * t114 + t115 * t2 - t55 * t8 + t56 * t9) * mrSges(7,3) + (mrSges(2,1) * t225 - t138 * mrSges(7,1) + mrSges(2,2) * t228 - t137 * mrSges(7,2) + t375 * (-t207 * t221 - t323) + (m(6) - t376) * t323 + t350 * t207 + (m(5) * t199 - m(7) * (-t199 - t243) - m(6) * (-t199 - t299) + t306 - t379) * t204) * g(1) + (-t81 * mrSges(5,3) - t373 + t374 + t366 / 0.2e1 + t367 / 0.2e1 - t371 + t372 - Ifges(5,4) * t327 + Ifges(6,3) * t330 - Ifges(5,2) * t331 + Ifges(7,3) * t332 + Ifges(7,5) * t339 + Ifges(7,6) * t341 - t116 / 0.2e1 + t160 * mrSges(5,1) + t370 / 0.2e1 - t382 / 0.2e1) * t162 + (-t80 * mrSges(5,3) + t260 * t249 / 0.2e1 + t144 * t251 / 0.2e1 + t267 + t268 + t365 + Ifges(5,1) * t327 + t247 * t330 + Ifges(5,4) * t331 + t117 / 0.2e1 + t160 * mrSges(5,2) + t383 / 0.2e1) * t164; m(3) * qJDD(2) - t114 * t18 - t115 * t19 + t224 * t158 + t227 * t159 + t55 * t60 + t56 * t61 + (-t187 * t224 + t189 * t227) * qJD(3) + (t124 + t355) * t176 - (t10 + t364) * t238 + t241 * t164 + t273 * t162 + (-t278 + t376) * g(3) + m(4) * (t118 * t224 + t119 * t227 + (-t153 * t224 + t154 * t227) * qJD(3)) + m(7) * (-t1 * t115 - t114 * t2 + t162 * t54 - t20 * t238 + t55 * t9 + t56 * t8) + m(5) * (-t162 * t80 + t164 * t81 + t176 * t38 + t238 * t37) + m(6) * (t162 * t71 - t164 * t244 + t176 * t245 - t238 * t36); (-qJD(5) * t244 + t192 * t245 + t196 * t36 - t40 * t51 - t41 * t52 - t71 * t88) * m(6) - t163 * t374 - t163 * t372 + t368 * t61 + (t1 * t121 + t120 * t2 + t182 * t20 + t368 * t8 + t369 * t9 - t54 * t62) * m(7) + t369 * t60 - (-Ifges(5,1) * t161 - t311 + t370) * t163 / 0.2e1 + (-m(7) * (t211 - t288) - t234 * t206 - m(6) * (t211 + t299) - t306 - m(5) * t211 + t385) * g(3) + ((t217 * t38 + t301 * t37) * pkin(3) - t160 * t272 + t80 * t88 - t81 * t89) * m(5) + (-Ifges(7,5) * t165 - Ifges(7,6) * t166) * t332 + (-Ifges(7,1) * t165 - Ifges(7,4) * t166) * t339 + (-Ifges(7,4) * t165 - Ifges(7,2) * t166) * t341 + (-t240 - t239 * qJD(1) / 0.2e1) * qJD(1) + t167 * t280 / 0.2e1 - t248 * t277 / 0.2e1 + t20 * (mrSges(7,1) * t242 + mrSges(7,2) * t177) + (Ifges(7,5) * t177 - Ifges(7,6) * t242) * t335 + (Ifges(7,4) * t177 - Ifges(7,2) * t242) * t346 + (Ifges(7,1) * t177 - Ifges(7,4) * t242) * t347 + (-t1 * t242 - t177 * t2 + t358 * t8 - t359 * t9) * mrSges(7,3) - t242 * t349 - t363 * t88 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + t163 * t373 + t161 * t365 + t163 * t371 + (mrSges(7,1) * t359 - mrSges(7,2) * t358) * t54 + (-t293 * t40 - t294 * t41 + t245) * mrSges(6,3) + t355 * t192 + t356 * qJD(5) + t380 * (t255 + t278 * pkin(3) * t224 + (m(7) * t222 - t384 + t386) * t206 + (t234 - t387) * t203) + (t269 - t189) * t153 + t125 * t266 + t36 * t253 - t144 * (Ifges(6,5) * t163 - t161 * t251) / 0.2e1 + t161 * t267 + t161 * t268 - (-Ifges(4,2) * t280 + t168 + t200) * t279 / 0.2e1 - t260 * (Ifges(6,6) * t163 - t161 * t249) / 0.2e1 - t89 * t151 + (-Ifges(5,2) * t163 + t117 - t156) * t330 + Ifges(5,5) * t131 - t126 * t272 - Ifges(5,6) * t130 + (t270 + t187) * t154 - t118 * mrSges(4,2) + t119 * mrSges(4,1) + t120 * t18 + t121 * t19 - t52 * t100 - t51 * t101 - t104 * t31 / 0.2e1 - t105 * t32 / 0.2e1 - t62 * t42 + t37 * mrSges(5,1) - t38 * mrSges(5,2) + t81 * t314 - t80 * t315 + Ifges(4,6) * t179 + Ifges(4,5) * t180 + t182 * t10 + t196 * t57 + t124 * t322 + t43 * t325 + t116 * t327 + (Ifges(6,3) * t163 - t161 * t247) * t331 + (Ifges(7,5) * t105 + Ifges(7,6) * t104 + Ifges(7,3) * t163) * t333 + (Ifges(6,5) * t216 + Ifges(6,6) * t219) * t334 + (Ifges(6,1) * t216 + t309) * t336 + (Ifges(6,2) * t219 + t310) * t337 + (Ifges(7,1) * t105 + Ifges(7,4) * t104 + Ifges(7,5) * t163) * t340 + (Ifges(7,4) * t105 + Ifges(7,2) * t104 + Ifges(7,6) * t163) * t342 + t216 * t343 - t165 * t344 - t166 * t345 + t177 * t348 - t160 * (mrSges(5,1) * t163 - mrSges(5,2) * t161) - qJD(3) * (-Ifges(5,5) * t161 - Ifges(5,6) * t163) / 0.2e1; -t242 * t18 + t177 * t19 + t216 * t63 + t219 * t64 - t359 * t61 - t358 * t60 - t273 * t163 + t241 * t161 + t261 + (-g(1) * t204 + g(2) * t207) * t278 + (t1 * t177 - t163 * t54 - t2 * t242 - t358 * t9 - t359 * t8) * m(7) + (t14 * t219 + t15 * t216 - t161 * t244 - t163 * t71) * m(6) + (t161 * t81 + t163 * t80 + t146) * m(5); t338 * t206 * g(3) - t260 * t100 + t144 * t101 - t377 * t60 + t85 * t61 + t10 + t57 + (-t377 * t9 + t8 * t85 + t20 - t362) * m(7) + (t144 * t40 - t260 * t41 + t36 - t362) * m(6); -t54 * (mrSges(7,1) * t85 + mrSges(7,2) * t377) + (Ifges(7,1) * t377 - t324) * t340 + t31 * t339 + (Ifges(7,5) * t377 - Ifges(7,6) * t85) * t333 - t8 * t60 + t9 * t61 - g(1) * (mrSges(7,1) * t139 - mrSges(7,2) * t140) - g(2) * (-mrSges(7,1) * t137 + mrSges(7,2) * t138) - g(3) * (-mrSges(7,1) * t202 - mrSges(7,2) * t205) * t203 + (t377 * t8 + t85 * t9) * mrSges(7,3) + t274 + (-Ifges(7,2) * t85 + t32 + t79) * t342 + t351;];
tau  = t5;
