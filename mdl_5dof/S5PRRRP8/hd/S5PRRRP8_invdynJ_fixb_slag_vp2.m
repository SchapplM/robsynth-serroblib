% Calculate vector of inverse dynamics joint torques for
% S5PRRRP8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP8_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP8_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:15
% EndTime: 2019-12-05 16:58:49
% DurationCPUTime: 13.81s
% Computational Cost: add. (3074->524), mult. (7306->717), div. (0->0), fcn. (5217->10), ass. (0->243)
t382 = Ifges(5,1) + Ifges(6,1);
t381 = Ifges(6,4) + Ifges(5,5);
t380 = Ifges(5,6) - Ifges(6,6);
t174 = sin(qJ(4));
t177 = cos(qJ(4));
t270 = t177 * qJD(3);
t175 = sin(qJ(3));
t280 = qJD(2) * t175;
t134 = t174 * t280 - t270;
t332 = t134 / 0.2e1;
t135 = qJD(3) * t174 + t177 * t280;
t384 = -t135 / 0.2e1;
t178 = cos(qJ(3));
t278 = qJD(2) * t178;
t165 = qJD(4) - t278;
t328 = t165 / 0.2e1;
t383 = qJD(3) / 0.2e1;
t379 = -Ifges(5,3) - Ifges(6,2);
t351 = -t174 * t380 + t177 * t381;
t305 = Ifges(6,5) * t174;
t308 = Ifges(5,4) * t174;
t349 = t177 * t382 + t305 - t308;
t269 = qJD(2) * qJD(3);
t140 = qJDD(2) * t175 + t178 * t269;
t275 = qJD(4) * t134;
t59 = qJDD(3) * t174 + t140 * t177 - t275;
t337 = t59 / 0.2e1;
t60 = qJD(4) * t135 - t177 * qJDD(3) + t140 * t174;
t336 = -t60 / 0.2e1;
t139 = qJDD(2) * t178 - t175 * t269;
t132 = qJDD(4) - t139;
t334 = t132 / 0.2e1;
t179 = cos(qJ(2));
t172 = sin(pkin(5));
t281 = qJD(2) * t172;
t248 = qJD(1) * t281;
t156 = t179 * t248;
t176 = sin(qJ(2));
t268 = qJDD(1) * t172;
t110 = t176 * t268 + t156;
t173 = cos(pkin(5));
t283 = qJD(1) * t173;
t378 = qJDD(2) * pkin(7) + qJD(3) * t283 + t110;
t168 = Ifges(4,4) * t278;
t131 = Ifges(5,4) * t134;
t306 = Ifges(6,5) * t134;
t364 = t135 * t382 + t165 * t381 - t131 + t306;
t130 = Ifges(6,5) * t135;
t48 = t165 * Ifges(6,6) + t134 * Ifges(6,3) + t130;
t377 = Ifges(4,1) * t280 + Ifges(4,5) * qJD(3) + t174 * t48 + t177 * t364 + t168;
t311 = Ifges(4,4) * t175;
t362 = Ifges(4,2) * t178;
t211 = t311 + t362;
t376 = t379 * t328 + t381 * t384 + t380 * t332 + Ifges(4,6) * t383 + qJD(2) * t211 / 0.2e1;
t284 = qJD(1) * t172;
t255 = t176 * t284;
t141 = qJD(2) * pkin(7) + t255;
t253 = t175 * t283;
t90 = t141 * t178 + t253;
t78 = qJD(3) * pkin(8) + t90;
t325 = pkin(3) * t178;
t147 = -pkin(8) * t175 - pkin(2) - t325;
t254 = t179 * t284;
t91 = qJD(2) * t147 - t254;
t23 = t174 * t91 + t177 * t78;
t15 = qJ(5) * t165 + t23;
t375 = -t15 * mrSges(6,2) - t23 * mrSges(5,3);
t374 = -m(5) - m(6);
t373 = (-Ifges(5,4) + Ifges(6,5)) * t60 + t382 * t59 + t381 * t132;
t372 = t139 / 0.2e1;
t371 = t140 / 0.2e1;
t370 = mrSges(3,2) - mrSges(4,3);
t369 = mrSges(5,3) + mrSges(6,2);
t29 = mrSges(5,1) * t132 - mrSges(5,3) * t59;
t30 = -t132 * mrSges(6,1) + t59 * mrSges(6,2);
t367 = t30 - t29;
t31 = -mrSges(5,2) * t132 - mrSges(5,3) * t60;
t32 = -mrSges(6,2) * t60 + mrSges(6,3) * t132;
t366 = t31 + t32;
t313 = mrSges(5,3) * t134;
t85 = -mrSges(5,2) * t165 - t313;
t88 = -mrSges(6,2) * t134 + mrSges(6,3) * t165;
t315 = t85 + t88;
t312 = mrSges(5,3) * t135;
t86 = mrSges(5,1) * t165 - t312;
t87 = -mrSges(6,1) * t165 + mrSges(6,2) * t135;
t363 = t86 - t87;
t220 = mrSges(4,1) * t178 - mrSges(4,2) * t175;
t361 = mrSges(3,1) + t220;
t11 = mrSges(5,1) * t60 + mrSges(5,2) * t59;
t360 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t140 + t11;
t200 = pkin(4) * t174 - qJ(5) * t177;
t359 = qJD(4) * t200 - qJD(5) * t174 - t253 - (qJD(2) * t200 + t141) * t178;
t261 = mrSges(4,3) * t280;
t358 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t134 + mrSges(5,2) * t135 + t261;
t22 = -t174 * t78 + t177 * t91;
t357 = qJD(5) - t22;
t355 = -t175 * t379 + t178 * t351;
t354 = t175 * t381 + t178 * t349;
t216 = mrSges(6,1) * t174 - mrSges(6,3) * t177;
t218 = mrSges(5,1) * t174 + mrSges(5,2) * t177;
t89 = -t175 * t141 + t178 * t283;
t77 = -qJD(3) * pkin(3) - t89;
t24 = pkin(4) * t134 - qJ(5) * t135 + t77;
t353 = t24 * t216 + t77 * t218;
t352 = t174 * t381 + t177 * t380;
t304 = Ifges(6,5) * t177;
t307 = Ifges(5,4) * t177;
t350 = t174 * t382 - t304 + t307;
t348 = -t132 * t379 - t380 * t60 + t381 * t59;
t267 = qJDD(1) * t173;
t277 = qJD(3) * t175;
t19 = -t141 * t277 + t175 * t267 + t178 * t378;
t276 = qJD(3) * t178;
t20 = -t141 * t276 - t175 * t378 + t178 * t267;
t347 = -t175 * t20 + t178 * t19;
t17 = qJDD(3) * pkin(8) + t19;
t272 = qJD(4) * t177;
t274 = qJD(4) * t174;
t155 = t176 * t248;
t109 = t179 * t268 - t155;
t94 = -qJDD(2) * pkin(2) - t109;
t37 = -pkin(3) * t139 - pkin(8) * t140 + t94;
t3 = t177 * t17 + t174 * t37 + t91 * t272 - t274 * t78;
t4 = -qJD(4) * t23 - t17 * t174 + t177 * t37;
t346 = -t174 * t4 + t177 * t3;
t1 = qJ(5) * t132 + qJD(5) * t165 + t3;
t2 = -pkin(4) * t132 + qJDD(5) - t4;
t345 = t1 * t177 + t174 * t2;
t344 = -t59 * Ifges(6,5) / 0.2e1 - t132 * Ifges(6,6) / 0.2e1 + Ifges(5,4) * t337 + Ifges(5,6) * t334 + (Ifges(6,3) + Ifges(5,2)) * t336;
t201 = pkin(4) * t177 + qJ(5) * t174;
t217 = -t177 * mrSges(6,1) - t174 * mrSges(6,3);
t219 = mrSges(5,1) * t177 - mrSges(5,2) * t174;
t343 = -m(6) * t201 - mrSges(4,1) + t217 - t219;
t221 = pkin(3) * t175 - pkin(8) * t178;
t138 = t221 * qJD(3);
t271 = qJD(4) * t178;
t40 = pkin(7) * (t174 * t277 - t177 * t271) + t138 * t177 - t147 * t274;
t229 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t341 = -m(5) * t77 - t358;
t227 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t292 = t172 * t178;
t121 = t173 * t175 + t176 * t292;
t299 = cos(pkin(9));
t232 = t299 * t176;
t171 = sin(pkin(9));
t294 = t171 * t179;
t116 = t173 * t232 + t294;
t233 = t172 * t299;
t67 = t116 * t178 - t175 * t233;
t231 = t299 * t179;
t295 = t171 * t176;
t118 = -t173 * t295 + t231;
t69 = t171 * t172 * t175 + t118 * t178;
t340 = -g(1) * t69 - g(2) * t67 - g(3) * t121;
t339 = -t4 * mrSges(5,1) + t2 * mrSges(6,1) + t3 * mrSges(5,2) - t1 * mrSges(6,3);
t180 = qJD(2) ^ 2;
t309 = Ifges(5,4) * t135;
t51 = -t134 * Ifges(5,2) + t165 * Ifges(5,6) + t309;
t338 = -t51 / 0.2e1;
t335 = t60 / 0.2e1;
t333 = -t134 / 0.2e1;
t330 = t135 / 0.2e1;
t327 = t174 / 0.2e1;
t317 = -qJD(2) / 0.2e1;
t316 = qJD(4) / 0.2e1;
t310 = Ifges(4,4) * t178;
t18 = -qJDD(3) * pkin(3) - t20;
t302 = t175 * t18;
t137 = t221 * qJD(2);
t46 = t174 * t137 + t177 * t89;
t115 = -t173 * t231 + t295;
t298 = t115 * t175;
t117 = t173 * t294 + t232;
t297 = t117 * t175;
t296 = t147 * t177;
t293 = t172 * t176;
t291 = t172 * t179;
t289 = t174 * t178;
t287 = t177 * t178;
t286 = t178 * t179;
t285 = pkin(2) * t291 + pkin(7) * t293;
t97 = pkin(7) * t287 + t174 * t147;
t279 = qJD(2) * t176;
t273 = qJD(4) * t175;
t260 = mrSges(4,3) * t278;
t259 = t175 * t291;
t258 = t174 * t291;
t252 = t172 * t279;
t251 = t179 * t281;
t250 = t174 * t276;
t237 = t272 / 0.2e1;
t236 = -t115 * pkin(2) + pkin(7) * t116;
t235 = -t117 * pkin(2) + pkin(7) * t118;
t230 = t269 / 0.2e1;
t225 = t178 * t254;
t210 = -Ifges(5,2) * t174 + t307;
t209 = Ifges(5,2) * t177 + t308;
t206 = Ifges(4,5) * t178 - Ifges(4,6) * t175;
t203 = Ifges(6,3) * t174 + t304;
t202 = -Ifges(6,3) * t177 + t305;
t45 = t137 * t177 - t174 * t89;
t196 = pkin(7) + t200;
t72 = t121 * t174 + t177 * t291;
t120 = -t173 * t178 + t175 * t293;
t142 = -qJD(2) * pkin(2) - t254;
t192 = t142 * (mrSges(4,1) * t175 + mrSges(4,2) * t178);
t191 = t175 * (Ifges(4,1) * t178 - t311);
t93 = (t174 * t176 + t177 * t286) * t172;
t184 = Ifges(5,6) * t175 + t178 * t210;
t183 = Ifges(6,6) * t175 + t178 * t203;
t39 = t174 * t138 + t147 * t272 + (-t174 * t271 - t175 * t270) * pkin(7);
t151 = -qJD(3) * mrSges(4,2) + t260;
t143 = -pkin(3) - t201;
t136 = t220 * qJD(2);
t111 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t139;
t106 = t196 * t175;
t96 = -pkin(7) * t289 + t296;
t84 = qJD(1) * t93;
t83 = t174 * t225 - t177 * t255;
t80 = -t296 + (pkin(7) * t174 + pkin(4)) * t178;
t79 = -qJ(5) * t178 + t97;
t74 = -mrSges(4,1) * t139 + mrSges(4,2) * t140;
t73 = t121 * t177 - t258;
t71 = -qJD(3) * t120 + t178 * t251;
t70 = qJD(3) * t121 + t175 * t251;
t68 = -t118 * t175 + t171 * t292;
t66 = -t116 * t175 - t178 * t233;
t64 = mrSges(6,1) * t134 - mrSges(6,3) * t135;
t63 = pkin(4) * t135 + qJ(5) * t134;
t38 = (qJD(4) * t201 - qJD(5) * t177) * t175 + t196 * t276;
t36 = -pkin(4) * t280 - t45;
t35 = qJ(5) * t280 + t46;
t33 = -pkin(4) * t277 - t40;
t27 = -t117 * t177 + t174 * t69;
t25 = -t115 * t177 + t174 * t67;
t21 = qJ(5) * t277 - qJD(5) * t178 + t39;
t14 = -pkin(4) * t165 + t357;
t13 = -qJD(4) * t72 + t174 * t252 + t177 * t71;
t12 = -qJD(4) * t258 + t121 * t272 + t174 * t71 - t177 * t252;
t10 = mrSges(6,1) * t60 - mrSges(6,3) * t59;
t5 = pkin(4) * t60 - qJ(5) * t59 - qJD(5) * t135 + t18;
t6 = [m(2) * qJDD(1) + t121 * t111 + t71 * t151 + t366 * t73 + t367 * t72 + t315 * t13 - t363 * t12 + (t64 + t358) * t70 + (t10 + t360) * t120 + (-m(2) - m(3) - m(4) + t374) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t180 - t74) * t179 + (-mrSges(3,1) * t180 - mrSges(3,2) * qJDD(2) - qJD(2) * t136) * t176) * t172 + m(4) * (-t120 * t20 + t121 * t19 - t70 * t89 + t71 * t90 + (t142 * t279 - t179 * t94) * t172) + m(3) * (qJDD(1) * t173 ^ 2 + (t109 * t179 + t110 * t176) * t172) + m(6) * (t1 * t73 + t12 * t14 + t120 * t5 + t13 * t15 + t2 * t72 + t24 * t70) + m(5) * (-t12 * t22 + t120 * t18 + t13 * t23 + t3 * t73 - t4 * t72 + t70 * t77); (-t276 * t89 + t347) * mrSges(4,3) + t363 * t83 - t315 * t84 + (t22 * mrSges(5,1) - t14 * mrSges(6,1) - t23 * mrSges(5,2) - t90 * mrSges(4,3) + t15 * mrSges(6,3) - t376) * t277 + qJDD(3) * (Ifges(4,5) * t175 + Ifges(4,6) * t178) + t310 * t371 + t211 * t372 + (t3 * t97 + t4 * t96 + (-t84 + t39) * t23 + (t83 + t40) * t22) * m(5) + (-(t142 * t176 + (-t175 * t89 + t178 * t90) * t179) * t284 - pkin(2) * t94) * m(4) + (t358 * t276 + t360 * t175 + t111 * t178 - t277 * t151 + (t276 * t77 + t302) * m(5) + ((-t175 * t90 - t178 * t89) * qJD(3) + t347) * m(4)) * pkin(7) + t96 * t29 + t97 * t31 + t106 * t10 + t80 * t30 + t39 * t85 + t40 * t86 + t33 * t87 + t21 * t88 - pkin(2) * t74 + t79 * t32 + t38 * t64 + (-t1 * mrSges(6,2) - t3 * mrSges(5,3) - t344) * t174 * t175 - t94 * t220 + (Ifges(4,1) * t140 + Ifges(4,4) * t372 + t203 * t335 + t210 * t336 + t5 * t216 - t230 * t362 + t48 * t237 + t334 * t351 + t337 * t349) * t175 + (-m(4) * t235 + t374 * (-pkin(8) * t297 - t117 * t325 + t235) - t229 * (-t117 * t287 + t118 * t174) + t227 * (-t117 * t289 - t118 * t177) + t369 * t297 + t370 * t118 + t361 * t117) * g(1) + (-m(4) * t236 + t374 * (-pkin(8) * t298 - t115 * t325 + t236) - t229 * (-t115 * t287 + t116 * t174) + t227 * (-t115 * t289 - t116 * t177) + t369 * t298 + t370 * t116 + t361 * t115) * g(2) + (-m(4) * t285 + t374 * (t172 * pkin(3) * t286 + pkin(8) * t259 + t285) - t229 * t93 + t227 * (-t177 * t293 + t178 * t258) - t369 * t259 + (t176 * t370 - t179 * t361) * t172) * g(3) + (-t110 + t156) * mrSges(3,2) + (-m(6) * t24 + t341 - t64) * t175 * t254 + t377 * t276 / 0.2e1 + (t373 / 0.2e1 + mrSges(6,2) * t2 - mrSges(5,3) * t4) * t175 * t177 + (mrSges(5,2) * t77 + mrSges(6,2) * t14 - mrSges(5,3) * t22 - mrSges(6,3) * t24) * (-t174 * t273 + t178 * t270) - (t174 * t364 + t177 * t51) * t273 / 0.2e1 + (Ifges(4,4) * t371 + Ifges(4,2) * t372 - Ifges(5,6) * t336 - Ifges(6,6) * t335 + t310 * t230 + t334 * t379 - t337 * t381 + t339) * t178 + (t1 * t79 + t106 * t5 + t2 * t80 + t24 * t38 + (-t84 + t21) * t15 + (-t83 + t33) * t14) * m(6) + (t109 + t155) * mrSges(3,1) + (-t202 * t332 - t209 * t333 - t328 * t352 - t330 * t350) * t273 - t348 * t178 / 0.2e1 + (mrSges(5,1) * t77 + mrSges(6,1) * t24 + t375) * (t175 * t272 + t250) + t191 * t230 + t218 * t302 - t151 * t225 + (t183 * t332 + t184 * t333 + t206 * t383 + t355 * t328 + t354 * t330 + t192) * qJD(3) + Ifges(3,3) * qJDD(2) + t136 * t255 + t250 * t338; t364 * t237 + (-t151 + t260) * t89 + t376 * t280 + (t316 * t351 + t317 * t355) * t165 + (t14 * t272 + t340 + t345) * mrSges(6,2) - t363 * pkin(8) * t272 + t366 * pkin(8) * t177 + t367 * pkin(8) * t174 + Ifges(4,5) * t140 + t143 * t10 + Ifges(4,6) * t139 - t180 * t191 / 0.2e1 - t46 * t85 - t45 * t86 - t36 * t87 - t35 * t88 - t18 * t219 + t5 * t217 + t373 * t327 + (mrSges(4,2) * t121 + t374 * (-t120 * pkin(3) + pkin(8) * t121) - t343 * t120) * g(3) + (mrSges(4,2) * t69 + t374 * (t68 * pkin(3) + pkin(8) * t69) + t343 * t68) * g(1) + (mrSges(4,2) * t67 + t374 * (t66 * pkin(3) + pkin(8) * t67) + t343 * t66) * g(2) + (t316 * t349 + t317 * t354) * t135 - t19 * mrSges(4,2) + t20 * mrSges(4,1) - pkin(3) * t11 + t359 * t64 + (t143 * t5 + ((t14 * t177 - t15 * t174) * qJD(4) + t345) * pkin(8) - t14 * t36 - t15 * t35 + t359 * t24) * m(6) - (-Ifges(4,2) * t280 + t168 + t377) * t278 / 0.2e1 + (-t22 * t272 + t340 + t346) * mrSges(5,3) + t344 * t177 + (t261 + t341) * t90 + t352 * t334 + (t327 * t51 - t353) * t278 + t353 * qJD(4) + (-t210 / 0.2e1 + t203 / 0.2e1) * t275 + (-t22 * t45 - t23 * t46 - pkin(3) * t18 + ((-t23 * t174 - t22 * t177) * qJD(4) + t346) * pkin(8)) * m(5) + t350 * t337 + (-t315 * pkin(8) + t338 + t48 / 0.2e1 + t375) * t274 + t202 * t335 + Ifges(4,3) * qJDD(3) + ((t184 / 0.2e1 - t183 / 0.2e1) * t134 - t192 - t22 * (mrSges(5,1) * t175 - mrSges(5,3) * t287) - t14 * (-mrSges(6,1) * t175 + mrSges(6,2) * t287) - t23 * (-mrSges(5,2) * t175 - mrSges(5,3) * t289) - t15 * (-mrSges(6,2) * t289 + mrSges(6,3) * t175)) * qJD(2) + t209 * t336 - t206 * t269 / 0.2e1; (t134 * t14 + t135 * t15) * mrSges(6,2) - t24 * (mrSges(6,1) * t135 + mrSges(6,3) * t134) - t77 * (mrSges(5,1) * t135 - mrSges(5,2) * t134) + qJD(5) * t88 - t63 * t64 + (t312 + t363) * t23 + (-t313 - t315) * t22 - pkin(4) * t30 + qJ(5) * t32 - (-t134 * t381 - t135 * t380) * t165 / 0.2e1 + (-t134 * t382 + t130 - t309 + t48) * t384 + (-Ifges(5,2) * t135 - t131 + t364) * t332 - t339 + t51 * t330 + (Ifges(6,3) * t135 - t306) * t333 + (-pkin(4) * t2 + qJ(5) * t1 - t14 * t23 + t15 * t357 - t24 * t63) * m(6) + (t227 * t73 + t229 * t72) * g(3) + (t227 * (t115 * t174 + t177 * t67) + t229 * t25) * g(2) + (t227 * (t117 * t174 + t177 * t69) + t229 * t27) * g(1) + t348; t135 * t64 - t165 * t88 + (-g(1) * t27 - g(2) * t25 - g(3) * t72 + t24 * t135 - t15 * t165 + t2) * m(6) + t30;];
tau = t6;
