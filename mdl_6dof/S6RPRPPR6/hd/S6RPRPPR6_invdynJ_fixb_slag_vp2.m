% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:53:22
% EndTime: 2019-03-09 02:53:44
% DurationCPUTime: 15.37s
% Computational Cost: add. (7986->659), mult. (16118->898), div. (0->0), fcn. (11177->14), ass. (0->286)
t234 = cos(qJ(3));
t399 = t234 / 0.2e1;
t223 = qJ(3) + pkin(9);
t214 = sin(t223);
t216 = cos(t223);
t231 = sin(qJ(3));
t263 = mrSges(4,1) * t231 + mrSges(4,2) * t234;
t398 = mrSges(5,1) * t214 + t263 + (mrSges(5,2) - mrSges(7,3)) * t216;
t226 = sin(pkin(9));
t311 = cos(pkin(9));
t271 = t311 * t234;
t290 = qJD(1) * t231;
t165 = qJD(1) * t271 - t226 * t290;
t225 = sin(pkin(10));
t227 = cos(pkin(10));
t136 = qJD(3) * t225 + t165 * t227;
t230 = sin(qJ(6));
t233 = cos(qJ(6));
t266 = t227 * qJD(3) - t165 * t225;
t397 = -t136 * t230 + t233 * t266;
t77 = t136 * t233 + t230 * t266;
t264 = mrSges(4,1) * t234 - mrSges(4,2) * t231;
t337 = pkin(3) * t234;
t333 = t227 * pkin(5);
t204 = pkin(4) + t333;
t222 = pkin(10) + qJ(6);
t213 = sin(t222);
t215 = cos(t222);
t261 = -mrSges(6,1) * t227 + mrSges(6,2) * t225;
t242 = m(6) * pkin(4) - t261;
t376 = m(7) * t204 + mrSges(7,1) * t215 - mrSges(7,2) * t213 + t242;
t396 = m(5) * t337 - mrSges(5,2) * t214 + t264 + (mrSges(5,1) + t376) * t216;
t282 = qJD(1) * qJD(3);
t183 = qJDD(1) * t234 - t231 * t282;
t184 = -qJDD(1) * t231 - t234 * t282;
t131 = t183 * t311 + t226 * t184;
t110 = qJDD(3) * t227 - t131 * t225;
t111 = qJDD(3) * t225 + t131 * t227;
t24 = qJD(6) * t397 + t110 * t230 + t111 * t233;
t362 = t24 / 0.2e1;
t25 = -qJD(6) * t77 + t110 * t233 - t111 * t230;
t361 = t25 / 0.2e1;
t353 = -m(3) - m(4);
t395 = -m(7) - m(5);
t394 = m(7) + m(6);
t351 = t110 / 0.2e1;
t350 = t111 / 0.2e1;
t130 = t183 * t226 - t311 * t184;
t129 = qJDD(6) + t130;
t349 = t129 / 0.2e1;
t348 = t130 / 0.2e1;
t232 = sin(qJ(1));
t335 = g(1) * t232;
t176 = t226 * t234 + t231 * t311;
t161 = t176 * qJD(1);
t158 = qJD(6) + t161;
t393 = t136 * Ifges(6,5) + t77 * Ifges(7,5) + Ifges(6,6) * t266 + Ifges(7,6) * t397 + t161 * Ifges(6,3) + t158 * Ifges(7,3);
t338 = pkin(3) * t226;
t202 = qJ(5) + t338;
t332 = pkin(8) + t202;
t171 = t332 * t225;
t172 = t332 * t227;
t120 = -t171 * t233 - t172 * t230;
t250 = t225 * t230 - t227 * t233;
t307 = t161 * t227;
t236 = -pkin(1) - pkin(7);
t194 = qJD(1) * t236 + qJD(2);
t300 = t194 * t231;
t155 = -qJ(4) * t290 + t300;
t149 = t226 * t155;
t185 = t234 * t194;
t289 = qJD(1) * t234;
t156 = -qJ(4) * t289 + t185;
t107 = t156 * t311 - t149;
t277 = pkin(3) * t289;
t109 = pkin(4) * t165 + qJ(5) * t161 + t277;
t51 = -t107 * t225 + t227 * t109;
t32 = pkin(5) * t165 + pkin(8) * t307 + t51;
t308 = t161 * t225;
t52 = t227 * t107 + t225 * t109;
t37 = pkin(8) * t308 + t52;
t392 = -qJD(5) * t250 + qJD(6) * t120 - t230 * t32 - t233 * t37;
t121 = -t171 * t230 + t172 * t233;
t178 = t225 * t233 + t227 * t230;
t391 = -qJD(5) * t178 - qJD(6) * t121 + t230 * t37 - t233 * t32;
t288 = qJD(3) * t231;
t162 = -qJD(3) * t271 + t226 * t288;
t167 = t250 * qJD(6);
t387 = t250 * qJD(1) + t162 * t178 + t167 * t176;
t168 = t178 * qJD(6);
t386 = -t178 * qJD(1) + t162 * t250 - t168 * t176;
t123 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t131;
t58 = -t110 * mrSges(6,1) + t111 * mrSges(6,2);
t385 = t58 - t123;
t330 = mrSges(5,3) * t165;
t384 = qJD(3) * mrSges(5,1) + mrSges(6,1) * t266 - mrSges(6,2) * t136 - t330;
t383 = t214 * t232;
t235 = cos(qJ(1));
t334 = g(2) * t235;
t372 = t334 - t335;
t382 = t216 * t372;
t101 = t178 * t161;
t379 = t101 + t168;
t102 = t250 * t161;
t378 = t102 + t167;
t220 = t231 * pkin(3);
t228 = -qJ(4) - pkin(7);
t377 = t235 * t220 + t232 * t228;
t283 = qJD(1) * qJD(2);
t195 = qJDD(1) * qJ(2) + t283;
t193 = qJDD(1) * t236 + qJDD(2);
t138 = t234 * t193 - t194 * t288;
t287 = qJD(3) * t234;
t139 = t231 * t193 + t194 * t287;
t251 = t138 * t234 + t139 * t231;
t97 = -mrSges(6,2) * t161 + mrSges(6,3) * t266;
t98 = mrSges(6,1) * t161 - mrSges(6,3) * t136;
t374 = -t225 * t98 + t227 * t97;
t62 = -mrSges(6,2) * t130 + mrSges(6,3) * t110;
t63 = mrSges(6,1) * t130 - mrSges(6,3) * t111;
t373 = -t225 * t63 + t227 * t62;
t281 = qJD(1) * qJD(4);
t103 = qJ(4) * t184 - t231 * t281 + t139;
t91 = qJDD(3) * pkin(3) - qJ(4) * t183 - t234 * t281 + t138;
t50 = t311 * t103 + t226 * t91;
t43 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t50;
t141 = -pkin(3) * t184 + qJDD(4) + t195;
t44 = pkin(4) * t130 - qJ(5) * t131 - qJD(5) * t165 + t141;
t14 = -t225 * t43 + t227 * t44;
t15 = t225 * t44 + t227 * t43;
t253 = -t14 * t225 + t15 * t227;
t328 = Ifges(4,4) * t234;
t371 = (-Ifges(4,1) * t231 - t328) * t399 + qJ(2) * t264;
t274 = m(6) * qJ(5) + mrSges(6,3);
t190 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t290;
t191 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t289;
t370 = (t190 * t234 - t191 * t231) * qJD(3);
t163 = t176 * qJD(3);
t175 = t226 * t231 - t271;
t49 = -t226 * t103 + t311 * t91;
t152 = qJD(3) * pkin(3) + t156;
t94 = t152 * t311 - t149;
t272 = t311 * t155;
t95 = t226 * t152 + t272;
t369 = t162 * t95 + t163 * t94 + t175 * t49 - t176 * t50;
t368 = m(4) * t251 + t234 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t183) + t231 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t184);
t229 = -pkin(8) - qJ(5);
t295 = t216 * t229;
t367 = -t242 * t214 + t216 * t274 + mrSges(2,2) - mrSges(3,3) - m(7) * (t204 * t214 + t295) - t398;
t260 = -mrSges(6,1) * t225 - mrSges(6,2) * t227;
t366 = m(7) * pkin(5) * t225 + mrSges(2,1) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) - t260;
t365 = qJD(1) ^ 2;
t364 = Ifges(7,4) * t362 + Ifges(7,2) * t361 + Ifges(7,6) * t349;
t363 = Ifges(7,1) * t362 + Ifges(7,4) * t361 + Ifges(7,5) * t349;
t339 = Ifges(7,4) * t77;
t29 = Ifges(7,2) * t397 + Ifges(7,6) * t158 + t339;
t360 = t29 / 0.2e1;
t73 = Ifges(7,4) * t397;
t30 = Ifges(7,1) * t77 + Ifges(7,5) * t158 + t73;
t359 = t30 / 0.2e1;
t358 = Ifges(6,4) * t350 + Ifges(6,2) * t351 + Ifges(6,6) * t348;
t357 = -t397 / 0.2e1;
t356 = t397 / 0.2e1;
t355 = -t77 / 0.2e1;
t354 = t77 / 0.2e1;
t347 = -t158 / 0.2e1;
t346 = t158 / 0.2e1;
t345 = -t161 / 0.2e1;
t344 = t161 / 0.2e1;
t341 = t165 / 0.2e1;
t340 = t176 / 0.2e1;
t85 = qJD(3) * qJ(5) + t95;
t186 = pkin(3) * t290 + qJD(1) * qJ(2) + qJD(4);
t96 = pkin(4) * t161 - qJ(5) * t165 + t186;
t48 = t225 * t96 + t227 * t85;
t292 = qJ(4) - t236;
t153 = -qJD(4) * t234 + t288 * t292;
t188 = t292 * t234;
t154 = -qJD(3) * t188 - qJD(4) * t231;
t105 = t226 * t153 + t154 * t311;
t196 = pkin(3) * t287 + qJD(2);
t74 = -pkin(4) * t162 + qJ(5) * t163 + qJD(5) * t175 + t196;
t40 = t227 * t105 + t225 * t74;
t331 = mrSges(5,3) * t161;
t329 = Ifges(4,4) * t231;
t327 = Ifges(5,4) * t165;
t326 = Ifges(6,4) * t225;
t325 = Ifges(6,4) * t227;
t46 = -qJDD(3) * pkin(4) + qJDD(5) - t49;
t320 = t175 * t46;
t306 = t163 * t225;
t305 = t163 * t227;
t304 = t175 * t225;
t303 = t175 * t227;
t299 = t213 * t232;
t298 = t213 * t235;
t297 = t215 * t232;
t296 = t215 * t235;
t206 = qJ(2) + t220;
t125 = pkin(4) * t176 + qJ(5) * t175 + t206;
t187 = t292 * t231;
t133 = -t187 * t311 - t226 * t188;
t65 = t225 * t125 + t227 * t133;
t291 = t235 * pkin(1) + t232 * qJ(2);
t285 = qJDD(1) * mrSges(3,2);
t284 = -m(5) - t394;
t279 = Ifges(7,5) * t24 + Ifges(7,6) * t25 + Ifges(7,3) * t129;
t34 = -mrSges(7,1) * t397 + mrSges(7,2) * t77;
t278 = -t34 + t384;
t276 = t232 * t220 + t291;
t275 = t311 * pkin(3);
t7 = -t25 * mrSges(7,1) + t24 * mrSges(7,2);
t219 = t235 * qJ(2);
t273 = -pkin(1) * t232 + t219;
t47 = -t225 * t85 + t227 * t96;
t270 = -t282 / 0.2e1;
t269 = t130 * mrSges(5,1) + t131 * mrSges(5,2);
t39 = -t105 * t225 + t227 * t74;
t267 = (t195 + t283) * qJ(2);
t64 = t227 * t125 - t133 * t225;
t104 = -t311 * t153 + t154 * t226;
t106 = t156 * t226 + t272;
t132 = -t187 * t226 + t311 * t188;
t205 = -t275 - pkin(4);
t259 = t234 * Ifges(4,1) - t329;
t258 = -Ifges(6,1) * t227 + t326;
t257 = -t231 * Ifges(4,2) + t328;
t256 = Ifges(6,2) * t225 - t325;
t255 = -Ifges(4,5) * t231 - Ifges(4,6) * t234;
t254 = -Ifges(6,5) * t227 + Ifges(6,6) * t225;
t252 = t225 * t47 - t227 * t48;
t26 = pkin(5) * t161 - pkin(8) * t136 + t47;
t33 = pkin(8) * t266 + t48;
t9 = -t230 * t33 + t233 * t26;
t10 = t230 * t26 + t233 * t33;
t45 = pkin(5) * t176 + pkin(8) * t303 + t64;
t53 = pkin(8) * t304 + t65;
t16 = -t230 * t53 + t233 * t45;
t17 = t230 * t45 + t233 * t53;
t147 = -qJD(3) * mrSges(5,2) - t331;
t247 = -t147 - t374;
t84 = -qJD(3) * pkin(4) + qJD(5) - t94;
t246 = t84 * t260;
t244 = t231 * (-Ifges(4,2) * t234 - t329);
t113 = t178 * t175;
t211 = -pkin(1) * qJDD(1) + qJDD(2);
t200 = t232 * t337;
t189 = t205 - t333;
t180 = t263 * qJD(1);
t170 = Ifges(4,5) * qJD(3) + qJD(1) * t259;
t169 = Ifges(4,6) * qJD(3) + qJD(1) * t257;
t157 = Ifges(5,4) * t161;
t146 = t214 * t296 - t299;
t145 = t214 * t298 + t297;
t144 = t214 * t297 + t298;
t143 = -t214 * t299 + t296;
t124 = mrSges(5,1) * t161 + mrSges(5,2) * t165;
t122 = -qJDD(3) * mrSges(5,2) - mrSges(5,3) * t130;
t119 = t165 * Ifges(5,1) + Ifges(5,5) * qJD(3) - t157;
t118 = -t161 * Ifges(5,2) + Ifges(5,6) * qJD(3) + t327;
t115 = t250 * t175;
t114 = t250 * t176;
t112 = t178 * t176;
t92 = -pkin(5) * t304 + t132;
t70 = -pkin(5) * t308 + t106;
t69 = -pkin(5) * t306 + t104;
t68 = t136 * Ifges(6,1) + Ifges(6,4) * t266 + t161 * Ifges(6,5);
t67 = t136 * Ifges(6,4) + Ifges(6,2) * t266 + t161 * Ifges(6,6);
t61 = -pkin(5) * t266 + t84;
t60 = mrSges(7,1) * t158 - mrSges(7,3) * t77;
t59 = -mrSges(7,2) * t158 + mrSges(7,3) * t397;
t57 = t163 * t178 - t167 * t175;
t55 = qJD(6) * t113 + t163 * t250;
t36 = t111 * Ifges(6,1) + t110 * Ifges(6,4) + t130 * Ifges(6,5);
t31 = pkin(8) * t306 + t40;
t27 = -t110 * pkin(5) + t46;
t23 = -pkin(5) * t162 + pkin(8) * t305 + t39;
t19 = -mrSges(7,2) * t129 + mrSges(7,3) * t25;
t18 = mrSges(7,1) * t129 - mrSges(7,3) * t24;
t11 = pkin(8) * t110 + t15;
t8 = pkin(5) * t130 - pkin(8) * t111 + t14;
t4 = -qJD(6) * t17 + t23 * t233 - t230 * t31;
t3 = qJD(6) * t16 + t23 * t230 + t233 * t31;
t2 = -qJD(6) * t10 - t11 * t230 + t233 * t8;
t1 = qJD(6) * t9 + t11 * t233 + t230 * t8;
t5 = [(Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (Ifges(4,5) * t234 - Ifges(5,5) * t175 - Ifges(4,6) * t231 - Ifges(5,6) * t176) * qJDD(3) - t231 * (Ifges(4,4) * t183 + Ifges(4,2) * t184) / 0.2e1 + (-t146 * mrSges(7,1) + t145 * mrSges(7,2) - m(6) * (t219 + t377) - m(4) * t219 - m(3) * t273 + t395 * (t273 + t377) + (-m(4) * t236 + m(6) * pkin(1) + t366) * t232 + t367 * t235) * g(1) + (-m(6) * t276 - t144 * mrSges(7,1) - t143 * mrSges(7,2) + t353 * t291 + t395 * (-t228 * t235 + t276) + (-m(4) * pkin(7) + m(6) * t228 - t366) * t235 + t367 * t232) * g(2) + t369 * mrSges(5,3) + t368 * t236 + (-m(5) * t94 + m(6) * t84 - t384) * t104 + (-m(5) * t49 + m(6) * t46 + t385) * t132 + t266 * (-Ifges(6,6) * t162 + t163 * t256) / 0.2e1 + m(4) * t267 + m(7) * (t1 * t17 + t10 * t3 + t16 * t2 + t27 * t92 + t4 * t9 + t61 * t69) - pkin(1) * t285 - t169 * t287 / 0.2e1 - t170 * t288 / 0.2e1 + t260 * t320 + t47 * (-mrSges(6,1) * t162 + mrSges(6,3) * t305) + t67 * t306 / 0.2e1 + t48 * (mrSges(6,2) * t162 + mrSges(6,3) * t306) - t36 * t303 / 0.2e1 + t14 * (mrSges(6,1) * t176 + mrSges(6,3) * t303) + t15 * (-mrSges(6,2) * t176 + mrSges(6,3) * t304) - t68 * t305 / 0.2e1 + t40 * t97 + t39 * t98 + t92 * t7 + t371 * t282 + t211 * mrSges(3,2) + t196 * t124 + qJ(2) * (-mrSges(4,1) * t184 + mrSges(4,2) * t183) + t186 * (-mrSges(5,1) * t162 - mrSges(5,2) * t163) + t1 * (-mrSges(7,2) * t176 + mrSges(7,3) * t113) + t2 * (mrSges(7,1) * t176 - mrSges(7,3) * t115) + t141 * (mrSges(5,1) * t176 - mrSges(5,2) * t175) + t131 * (-Ifges(5,1) * t175 - Ifges(5,4) * t176) + qJD(2) * t180 + qJD(3) * (-Ifges(5,5) * t163 + Ifges(5,6) * t162) / 0.2e1 - t163 * t119 / 0.2e1 + t162 * t118 / 0.2e1 + t9 * (-mrSges(7,1) * t162 - mrSges(7,3) * t55) + t10 * (mrSges(7,2) * t162 + mrSges(7,3) * t57) + t105 * t147 + t133 * t122 + t244 * t270 + t236 * t370 + t206 * t269 + m(3) * (-pkin(1) * t211 + t267) + t184 * t257 / 0.2e1 + t136 * (-Ifges(6,5) * t162 + t163 * t258) / 0.2e1 + t183 * t259 / 0.2e1 + (Ifges(5,4) * t175 + Ifges(5,2) * t176 + Ifges(6,3) * t340) * t130 - t393 * t162 / 0.2e1 + (Ifges(6,5) * t111 + Ifges(6,6) * t110 + t279) * t340 + qJD(3) ^ 2 * t255 / 0.2e1 + (Ifges(7,1) * t55 + Ifges(7,4) * t57 - Ifges(7,5) * t162) * t354 + (Ifges(7,4) * t55 + Ifges(7,2) * t57 - Ifges(7,6) * t162) * t356 + t304 * t358 + t55 * t359 + t57 * t360 + (Ifges(7,4) * t115 + Ifges(7,2) * t113 + Ifges(7,6) * t176) * t361 + (Ifges(7,1) * t115 + Ifges(7,4) * t113 + Ifges(7,5) * t176) * t362 + t115 * t363 + t113 * t364 + (-Ifges(5,1) * t163 + Ifges(5,4) * t162) * t341 + (-Ifges(6,3) * t162 + t163 * t254) * t344 + (-Ifges(5,4) * t163 + Ifges(5,2) * t162) * t345 + (Ifges(7,5) * t55 + Ifges(7,6) * t57 - Ifges(7,3) * t162) * t346 + (Ifges(6,3) * t176 + t175 * t254) * t348 + (Ifges(7,5) * t115 + Ifges(7,6) * t113 + Ifges(7,3) * t176) * t349 + (Ifges(6,5) * t176 + t175 * t258) * t350 + (Ifges(6,6) * t176 + t175 * t256) * t351 + m(6) * (t14 * t64 + t15 * t65 + t39 * t47 + t40 * t48) + m(5) * (t105 * t95 + t133 * t50 + t141 * t206 + t186 * t196) + t69 * t34 + t3 * t59 + t4 * t60 + t61 * (-mrSges(7,1) * t57 + mrSges(7,2) * t55) + t64 * t63 + t65 * t62 - t251 * mrSges(4,3) + t16 * t18 + t17 * t19 + (t263 + 0.2e1 * mrSges(3,3)) * t195 + (Ifges(4,1) * t183 + Ifges(4,4) * t184) * t399 + t163 * t246 + t27 * (-mrSges(7,1) * t113 + mrSges(7,2) * t115); t285 - t112 * t18 - t114 * t19 + t387 * t60 + t386 * t59 + t370 + (qJ(2) * t353 - mrSges(3,3)) * t365 + (t122 + t373) * t176 + (t7 + t385) * t175 - t278 * t163 + t247 * t162 + (-m(5) * t186 - t225 * t97 - t227 * t98 - t124 - t180) * qJD(1) - m(5) * t369 + m(3) * t211 + t372 * (-t284 - t353) + (-t1 * t114 + t10 * t386 - t112 * t2 + t163 * t61 + t175 * t27 + t387 * t9) * m(7) + (t252 * t162 + t163 * t84 + t253 * t176 + t320 - (t225 * t48 + t227 * t47) * qJD(1)) * m(6) + t368; ((t226 * t50 + t311 * t49) * pkin(3) + t106 * t94 - t107 * t95 - t186 * t277) * m(5) + (-Ifges(5,2) * t165 + t119 - t157) * t344 + t384 * t106 + (mrSges(7,1) * t379 - mrSges(7,2) * t378) * t61 - t266 * (Ifges(6,6) * t165 + t161 * t256) / 0.2e1 + (-Ifges(7,1) * t167 - Ifges(7,4) * t168) * t354 + (-Ifges(7,4) * t167 - Ifges(7,2) * t168) * t356 + (-Ifges(7,5) * t167 - Ifges(7,6) * t168) * t346 - t94 * t331 - t161 * t246 - t102 * t30 / 0.2e1 + t170 * t290 / 0.2e1 + t169 * t289 / 0.2e1 + t95 * t330 + t191 * t300 - t124 * t277 + (-g(1) * t200 - qJD(5) * t252 - t106 * t84 + t202 * t253 + t205 * t46 - t47 * t51 - t48 * t52) * m(6) + t68 * t307 / 0.2e1 - t67 * t308 / 0.2e1 - t52 * t97 - t51 * t98 - t101 * t29 / 0.2e1 + (t244 / 0.2e1 - t371) * t365 + t205 * t58 + Ifges(4,5) * t183 + Ifges(4,6) * t184 - t186 * (mrSges(5,1) * t165 - mrSges(5,2) * t161) + t189 * t7 - qJD(3) * (-Ifges(5,5) * t161 - Ifges(5,6) * t165) / 0.2e1 - t139 * mrSges(4,2) - t107 * t147 + t138 * mrSges(4,1) - Ifges(5,6) * t130 + Ifges(5,5) * t131 + t255 * t270 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) - t136 * (Ifges(6,5) * t165 + t161 * t258) / 0.2e1 + t46 * t261 + t391 * t60 + t392 * t59 + (-t61 * t70 + t1 * t121 + t120 * t2 + t189 * t27 - g(1) * (-t229 * t383 + t200) + t391 * t9 + t392 * t10) * m(7) - (-Ifges(5,1) * t161 - t327 + t393) * t165 / 0.2e1 + (Ifges(7,1) * t102 + Ifges(7,4) * t101 + Ifges(7,5) * t165) * t355 + (Ifges(7,4) * t102 + Ifges(7,2) * t101 + Ifges(7,6) * t165) * t357 + t227 * t358 - t167 * t359 - t168 * t360 + t178 * t363 + t118 * t341 + (Ifges(6,3) * t165 + t161 * t254) * t345 + (Ifges(7,5) * t102 + Ifges(7,6) * t101 + Ifges(7,3) * t165) * t347 + (Ifges(6,5) * t225 + Ifges(6,6) * t227) * t348 + (Ifges(6,1) * t225 + t325) * t350 + (Ifges(6,2) * t227 + t326) * t351 + t122 * t338 - t70 * t34 + (-t307 * t47 - t308 * t48 + t253) * mrSges(6,3) + t373 * t202 + t374 * qJD(5) + t49 * mrSges(5,1) - t50 * mrSges(5,2) + (m(5) * t220 - m(7) * (-t220 - t295) - m(6) * (qJ(5) * t216 - t220) - t216 * mrSges(6,3) + t376 * t214 + t398) * g(3) + (t394 * t337 + (-m(7) * t229 + mrSges(7,3) + t274) * t214 + t396) * t334 + (-t214 * t274 - t396) * t335 + (-g(1) * t383 - t1 * t250 - t10 * t379 - t178 * t2 + t378 * t9) * mrSges(7,3) + t27 * (mrSges(7,1) * t250 + mrSges(7,2) * t178) + (Ifges(7,4) * t178 - Ifges(7,2) * t250) * t361 + (Ifges(7,1) * t178 - Ifges(7,4) * t250) * t362 + (Ifges(7,5) * t178 - Ifges(7,6) * t250) * t349 - t250 * t364 + t225 * t36 / 0.2e1 + t120 * t18 + t121 * t19 - t190 * t185 + t123 * t275 + t48 * mrSges(6,2) * t165 - t47 * mrSges(6,1) * t165 - t9 * mrSges(7,1) * t165 + t10 * mrSges(7,2) * t165; -t250 * t18 + t178 * t19 + t225 * t62 + t227 * t63 - t379 * t60 - t378 * t59 + t278 * t165 - t247 * t161 + t269 + (g(1) * t235 + g(2) * t232) * t284 + (t1 * t178 - t10 * t378 - t165 * t61 - t2 * t250 - t379 * t9) * m(7) + (t14 * t227 + t15 * t225 - t161 * t252 - t165 * t84) * m(6) + (t161 * t95 + t165 * t94 + t141) * m(5); -t394 * t214 * g(3) + t136 * t98 - t266 * t97 - t397 * t59 + t77 * t60 + t58 + t7 + (-t10 * t397 + t77 * t9 + t27 - t382) * m(7) + (t136 * t47 - t266 * t48 - t382 + t46) * m(6); -t1 * mrSges(7,2) + t2 * mrSges(7,1) - t61 * (mrSges(7,1) * t77 + mrSges(7,2) * t397) + (Ifges(7,1) * t397 - t339) * t355 + t29 * t354 + (Ifges(7,5) * t397 - Ifges(7,6) * t77) * t347 - t9 * t59 + t10 * t60 - g(1) * (mrSges(7,1) * t143 - mrSges(7,2) * t144) - g(2) * (mrSges(7,1) * t145 + mrSges(7,2) * t146) - g(3) * (-mrSges(7,1) * t213 - mrSges(7,2) * t215) * t216 + (t10 * t77 + t397 * t9) * mrSges(7,3) + t279 + (-Ifges(7,2) * t77 + t30 + t73) * t357;];
tau  = t5;
