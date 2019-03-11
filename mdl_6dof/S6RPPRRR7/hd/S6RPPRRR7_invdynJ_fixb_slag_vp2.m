% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:08
% EndTime: 2019-03-09 02:33:26
% DurationCPUTime: 10.81s
% Computational Cost: add. (11069->575), mult. (23264->751), div. (0->0), fcn. (17024->14), ass. (0->265)
t215 = sin(pkin(10));
t216 = cos(pkin(10));
t300 = t215 ^ 2 + t216 ^ 2;
t270 = t300 * mrSges(4,3);
t218 = -pkin(1) - qJ(3);
t398 = -qJD(1) * qJD(3) + qJDD(1) * t218;
t210 = pkin(10) + qJ(4);
t201 = qJ(5) + t210;
t190 = sin(t201);
t191 = cos(t201);
t397 = mrSges(6,1) * t190 + (mrSges(6,2) - mrSges(7,3)) * t191;
t221 = sin(qJ(4));
t225 = cos(qJ(4));
t308 = t216 * t225;
t243 = t215 * t221 - t308;
t219 = sin(qJ(6));
t340 = mrSges(7,2) * t219;
t223 = cos(qJ(6));
t342 = mrSges(7,1) * t223;
t396 = t340 - t342;
t395 = -m(4) - m(5);
t394 = -m(6) - m(7);
t271 = -pkin(5) * t190 + t191 * pkin(9);
t393 = m(7) * t271;
t211 = qJD(4) + qJD(5);
t220 = sin(qJ(5));
t224 = cos(qJ(5));
t177 = qJD(1) * t218 + qJD(2);
t268 = -pkin(7) * qJD(1) + t177;
t144 = t268 * t215;
t145 = t268 * t216;
t103 = t144 * t225 + t145 * t221;
t160 = -t225 * t215 - t221 * t216;
t150 = t160 * qJD(1);
t89 = pkin(8) * t150 + t103;
t319 = t224 * t89;
t102 = -t144 * t221 + t225 * t145;
t299 = qJD(1) * t215;
t151 = qJD(1) * t308 - t221 * t299;
t88 = -pkin(8) * t151 + t102;
t85 = qJD(4) * pkin(4) + t88;
t50 = t220 * t85 + t319;
t45 = pkin(9) * t211 + t50;
t197 = qJD(1) * qJ(2) + qJD(3);
t169 = pkin(3) * t299 + t197;
t124 = -pkin(4) * t150 + t169;
t245 = t150 * t220 + t224 * t151;
t263 = t224 * t150 - t151 * t220;
t63 = -pkin(5) * t263 - pkin(9) * t245 + t124;
t20 = -t219 * t45 + t223 * t63;
t392 = t20 * mrSges(7,1);
t21 = t219 * t63 + t223 * t45;
t391 = t21 * mrSges(7,2);
t336 = mrSges(6,3) * t245;
t93 = t211 * t223 - t219 * t245;
t94 = t211 * t219 + t223 * t245;
t344 = -mrSges(6,1) * t211 - mrSges(7,1) * t93 + mrSges(7,2) * t94 + t336;
t105 = qJD(6) - t263;
t70 = -mrSges(7,2) * t105 + mrSges(7,3) * t93;
t71 = mrSges(7,1) * t105 - mrSges(7,3) * t94;
t250 = -t219 * t71 + t223 * t70;
t337 = mrSges(6,3) * t263;
t98 = -mrSges(6,2) * t211 + t337;
t239 = -t250 - t98;
t121 = t160 * t220 - t224 * t243;
t345 = -pkin(7) + t218;
t166 = t345 * t215;
t167 = t345 * t216;
t123 = t225 * t166 + t221 * t167;
t226 = cos(qJ(1));
t287 = t191 * t340;
t341 = mrSges(6,2) * t190;
t390 = (-t287 - t341) * t226;
t222 = sin(qJ(1));
t284 = t191 * t342;
t310 = t191 * t222;
t311 = t190 * t222;
t389 = -mrSges(6,1) * t310 - mrSges(7,3) * t311 - (t284 - t287) * t222;
t213 = qJD(1) * qJD(2);
t179 = -qJDD(1) * qJ(2) - t213;
t388 = -t190 * t396 + t397;
t293 = qJD(6) * t223;
t297 = qJD(4) * t225;
t298 = qJD(4) * t221;
t152 = -t215 * t298 + t216 * t297;
t153 = t160 * qJD(4);
t244 = t224 * t160 + t220 * t243;
t78 = qJD(5) * t244 - t152 * t220 + t224 * t153;
t237 = t121 * t293 + t219 * t78;
t206 = qJDD(4) + qJDD(5);
t116 = qJD(1) * t153 - qJDD(1) * t243;
t117 = qJD(1) * qJD(4) * t243 + qJDD(1) * t160;
t61 = qJD(5) * t263 + t116 * t224 + t117 * t220;
t34 = qJD(6) * t93 + t206 * t219 + t223 * t61;
t62 = -qJD(5) * t245 - t116 * t220 + t117 * t224;
t60 = qJDD(6) - t62;
t17 = mrSges(7,1) * t60 - mrSges(7,3) * t34;
t35 = -qJD(6) * t94 + t206 * t223 - t219 * t61;
t18 = -mrSges(7,2) * t60 + mrSges(7,3) * t35;
t387 = -t219 * t17 + t223 * t18;
t386 = mrSges(6,1) * t191 + t190 * mrSges(7,3) + t284;
t385 = -g(1) * t222 + g(2) * t226;
t165 = qJDD(2) + t398;
t264 = -pkin(7) * qJDD(1) + t165;
t129 = t264 * t215;
t130 = t264 * t216;
t67 = -qJD(4) * t103 - t129 * t221 + t225 * t130;
t47 = qJDD(4) * pkin(4) - pkin(8) * t116 + t67;
t66 = t225 * t129 + t221 * t130 - t144 * t298 + t145 * t297;
t48 = pkin(8) * t117 + t66;
t15 = -qJD(5) * t50 - t220 * t48 + t224 * t47;
t12 = -pkin(5) * t206 - t15;
t16 = -mrSges(7,1) * t35 + mrSges(7,2) * t34;
t382 = -m(7) * t12 - t16;
t323 = t220 * t89;
t49 = t224 * t85 - t323;
t44 = -pkin(5) * t211 - t49;
t381 = -m(7) * t44 - t344;
t380 = -t102 * t153 - t103 * t152 + t160 * t66 + t243 * t67;
t379 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) - mrSges(3,2);
t295 = qJD(5) * t224;
t296 = qJD(5) * t220;
t14 = t220 * t47 + t224 * t48 + t85 * t295 - t296 * t89;
t11 = pkin(9) * t206 + t14;
t176 = qJDD(3) - t179;
t291 = qJDD(1) * t215;
t154 = pkin(3) * t291 + t176;
t92 = -pkin(4) * t117 + t154;
t19 = -pkin(5) * t62 - pkin(9) * t61 + t92;
t2 = qJD(6) * t20 + t11 * t223 + t19 * t219;
t3 = -qJD(6) * t21 - t11 * t219 + t19 * t223;
t378 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t76 = pkin(5) * t245 - pkin(9) * t263;
t202 = t215 * pkin(3);
t199 = sin(t210);
t200 = cos(t210);
t258 = -mrSges(5,1) * t199 - mrSges(5,2) * t200;
t260 = mrSges(4,1) * t215 + mrSges(4,2) * t216;
t377 = -m(5) * t202 + mrSges(2,2) - mrSges(3,3) + t258 - t260 + t393 - t397;
t252 = t20 * t223 + t21 * t219;
t231 = -qJD(6) * t252 + t2 * t223 - t219 * t3;
t294 = qJD(6) * t219;
t376 = m(7) * t231 - t71 * t293 - t70 * t294 + t387;
t375 = m(4) * t197 + m(5) * t169 - mrSges(5,1) * t150 + mrSges(5,2) * t151 + t260 * qJD(1);
t359 = -t211 / 0.2e1;
t364 = -t263 / 0.2e1;
t365 = -t105 / 0.2e1;
t367 = -t94 / 0.2e1;
t368 = -t93 / 0.2e1;
t374 = -t124 * mrSges(6,1) + Ifges(7,5) * t367 - Ifges(6,2) * t364 - Ifges(6,6) * t359 + Ifges(7,6) * t368 + Ifges(7,3) * t365 + t391 - t392;
t256 = mrSges(7,1) * t219 + mrSges(7,2) * t223;
t235 = t44 * t256;
t253 = Ifges(7,5) * t223 - Ifges(7,6) * t219;
t330 = Ifges(7,4) * t223;
t254 = -Ifges(7,2) * t219 + t330;
t331 = Ifges(7,4) * t219;
t255 = Ifges(7,1) * t223 - t331;
t91 = Ifges(7,4) * t93;
t42 = t94 * Ifges(7,1) + t105 * Ifges(7,5) + t91;
t321 = t223 * t42;
t334 = mrSges(7,3) * t223;
t335 = mrSges(7,3) * t219;
t358 = t219 / 0.2e1;
t363 = -t245 / 0.2e1;
t347 = t94 * Ifges(7,4);
t41 = t93 * Ifges(7,2) + t105 * Ifges(7,6) + t347;
t373 = -t124 * mrSges(6,2) + Ifges(6,1) * t363 + Ifges(6,5) * t359 + t20 * t334 + t21 * t335 + t253 * t365 + t254 * t368 + t255 * t367 - t235 - t321 / 0.2e1 + t41 * t358;
t372 = t34 / 0.2e1;
t371 = t35 / 0.2e1;
t369 = t60 / 0.2e1;
t366 = t94 / 0.2e1;
t362 = t245 / 0.2e1;
t360 = t151 / 0.2e1;
t357 = pkin(4) * t151;
t356 = pkin(4) * t199;
t355 = pkin(4) * t200;
t354 = pkin(4) * t220;
t353 = pkin(4) * t224;
t349 = t49 * mrSges(6,3);
t348 = t50 * mrSges(6,3);
t217 = -pkin(7) - qJ(3);
t339 = mrSges(5,3) * t150;
t338 = mrSges(5,3) * t151;
t333 = Ifges(5,4) * t151;
t332 = Ifges(6,4) * t245;
t318 = pkin(1) * qJDD(1);
t315 = t121 * t219;
t314 = t121 * t223;
t307 = t219 * t222;
t306 = t219 * t226;
t305 = t222 * t223;
t304 = t223 * t226;
t185 = qJ(2) + t202;
t303 = pkin(5) * t310 + pkin(9) * t311;
t290 = qJDD(1) * t216;
t302 = mrSges(4,1) * t291 + mrSges(4,2) * t290;
t301 = t226 * pkin(1) + t222 * qJ(2);
t289 = m(6) * t355;
t288 = Ifges(7,5) * t34 + Ifges(7,6) * t35 + Ifges(7,3) * t60;
t283 = t394 + t395;
t278 = t321 / 0.2e1;
t274 = -t62 * mrSges(6,1) + t61 * mrSges(6,2);
t273 = -t294 / 0.2e1;
t204 = t226 * qJ(2);
t272 = -pkin(1) * t222 + t204;
t269 = -t117 * mrSges(5,1) + t116 * mrSges(5,2);
t135 = pkin(4) * t152 + qJD(2);
t267 = t300 * t177;
t266 = t300 * t165;
t122 = -t166 * t221 + t225 * t167;
t134 = -pkin(4) * t160 + t185;
t261 = -pkin(5) * t191 - pkin(9) * t190;
t251 = -t20 * t219 + t21 * t223;
t100 = pkin(8) * t243 + t122;
t101 = pkin(8) * t160 + t123;
t69 = t100 * t220 + t101 * t224;
t72 = -pkin(5) * t244 - pkin(9) * t121 + t134;
t29 = t219 * t72 + t223 * t69;
t28 = -t219 * t69 + t223 * t72;
t249 = -t219 * t70 - t223 * t71;
t247 = t224 * t100 - t101 * t220;
t236 = t121 * t294 - t223 * t78;
t96 = qJD(3) * t243 - qJD(4) * t123;
t95 = qJD(3) * t160 - t166 * t298 + t167 * t297;
t230 = qJD(5) * t121 + t224 * t152 + t153 * t220;
t229 = -pkin(8) * t153 + t96;
t8 = t34 * Ifges(7,4) + t35 * Ifges(7,2) + t60 * Ifges(7,6);
t9 = t34 * Ifges(7,1) + t35 * Ifges(7,4) + t60 * Ifges(7,5);
t228 = -t14 * mrSges(6,2) - t3 * t335 + t2 * t334 + t12 * t396 + t15 * mrSges(6,1) + Ifges(6,3) * t206 + (Ifges(7,1) * t219 + t330) * t372 + (Ifges(7,2) * t223 + t331) * t371 + t41 * t273 + (Ifges(7,5) * t219 + Ifges(7,6) * t223) * t369 + Ifges(6,6) * t62 + Ifges(6,5) * t61 + t9 * t358 + t223 * t8 / 0.2e1 + (-t20 * t293 - t21 * t294) * mrSges(7,3) + (t278 + t235) * qJD(6) + (t105 * t253 + t254 * t93 + t255 * t94) * qJD(6) / 0.2e1;
t227 = qJD(1) ^ 2;
t207 = -pkin(8) + t217;
t198 = qJDD(2) - t318;
t195 = -pkin(5) - t353;
t168 = t202 + t356;
t146 = Ifges(5,4) * t150;
t143 = t190 * t304 - t307;
t142 = t190 * t306 + t305;
t141 = t190 * t305 + t306;
t140 = -t190 * t307 + t304;
t133 = qJD(4) * mrSges(5,1) - t338;
t132 = -qJD(4) * mrSges(5,2) + t339;
t109 = t151 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t146;
t108 = t150 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t333;
t107 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t117;
t106 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t116;
t104 = Ifges(6,4) * t263;
t83 = -pkin(8) * t152 + t95;
t75 = -mrSges(6,1) * t263 + mrSges(6,2) * t245;
t74 = Ifges(6,1) * t245 + t211 * Ifges(6,5) + t104;
t73 = Ifges(6,2) * t263 + t211 * Ifges(6,6) + t332;
t65 = t357 + t76;
t55 = -mrSges(6,2) * t206 + mrSges(6,3) * t62;
t54 = mrSges(6,1) * t206 - mrSges(6,3) * t61;
t52 = t224 * t88 - t323;
t51 = t220 * t88 + t319;
t40 = t94 * Ifges(7,5) + t93 * Ifges(7,6) + t105 * Ifges(7,3);
t36 = pkin(5) * t230 - pkin(9) * t78 + t135;
t27 = t219 * t76 + t223 * t49;
t26 = -t219 * t49 + t223 * t76;
t25 = t219 * t65 + t223 * t52;
t24 = -t219 * t52 + t223 * t65;
t22 = qJD(5) * t247 + t220 * t229 + t224 * t83;
t5 = -qJD(6) * t29 - t219 * t22 + t223 * t36;
t4 = qJD(6) * t28 + t219 * t36 + t22 * t223;
t1 = [-t230 * t391 + (-Ifges(7,1) * t236 - Ifges(7,4) * t237 + Ifges(7,5) * t230) * t366 + (Ifges(6,1) * t78 - Ifges(6,4) * t230) * t362 + t211 * (Ifges(6,5) * t78 - Ifges(6,6) * t230) / 0.2e1 + t124 * (mrSges(6,1) * t230 + mrSges(6,2) * t78) - t230 * t348 + t230 * t40 / 0.2e1 - t230 * t73 / 0.2e1 + t105 * (-Ifges(7,5) * t236 - Ifges(7,6) * t237 + Ifges(7,3) * t230) / 0.2e1 + t93 * (-Ifges(7,4) * t236 - Ifges(7,2) * t237 + Ifges(7,6) * t230) / 0.2e1 + t263 * (Ifges(6,4) * t78 - Ifges(6,2) * t230) / 0.2e1 + (Ifges(2,3) + Ifges(3,1)) * qJDD(1) + m(5) * (t102 * t96 + t103 * t95 + t122 * t67 + t123 * t66 + t154 * t185) + m(7) * (t2 * t29 + t20 * t5 + t21 * t4 + t28 * t3) + m(6) * (t124 * t135 + t134 * t92 + t14 * t69 + t22 * t50) - (-m(6) * t15 - t382 - t54) * t247 + t375 * qJD(2) + (-mrSges(5,1) * t154 + Ifges(5,4) * t116 + Ifges(5,2) * t117 + Ifges(5,6) * qJDD(4)) * t160 + (Ifges(4,1) * t216 - Ifges(4,4) * t215) * t290 + m(4) * (qJ(2) * t176 - qJD(3) * t267 + t218 * t266) + (Ifges(5,1) * t153 - Ifges(5,4) * t152) * t360 + (-t141 * mrSges(7,1) - t140 * mrSges(7,2) + t394 * (t222 * t168 - t207 * t226 + t301) + (-m(3) + t395) * t301 + (-m(4) * qJ(3) + m(5) * t217 - t379) * t226 + t377 * t222) * g(2) + (-m(3) * t272 - t143 * mrSges(7,1) + t142 * mrSges(7,2) + t394 * (t226 * t168 + t222 * t207 + t272) + t395 * t204 + (-m(4) * t218 - m(5) * (-pkin(1) + t217) + t379) * t222 + t377 * t226) * g(1) + (t198 - t318) * mrSges(3,2) + (t92 * mrSges(6,2) - t15 * mrSges(6,3) + Ifges(6,1) * t61 + Ifges(6,4) * t62 + Ifges(6,5) * t206 + t12 * t256 + t253 * t369 + t254 * t371 + t255 * t372 + t273 * t42) * t121 - 0.2e1 * t179 * mrSges(3,3) + t169 * (mrSges(5,1) * t152 + mrSges(5,2) * t153) + t150 * (Ifges(5,4) * t153 - Ifges(5,2) * t152) / 0.2e1 + qJD(4) * (Ifges(5,5) * t153 - Ifges(5,6) * t152) / 0.2e1 + t153 * t109 / 0.2e1 - t152 * t108 / 0.2e1 + t135 * t75 + t95 * t132 + t96 * t133 + t122 * t106 + t123 * t107 - t78 * t349 + t22 * t98 + t78 * t74 / 0.2e1 + t69 * t55 + t4 * t70 + t5 * t71 + t9 * t314 / 0.2e1 - t8 * t315 / 0.2e1 + t28 * t17 + t29 * t18 + qJ(2) * t302 + (-t2 * t315 + t20 * t236 - t21 * t237 - t3 * t314) * mrSges(7,3) + m(3) * (-pkin(1) * t198 + (-t179 + t213) * qJ(2)) - (Ifges(4,4) * t216 - Ifges(4,2) * t215) * t291 + t134 * t274 + t185 * t269 + t176 * t260 - t237 * t41 / 0.2e1 - (Ifges(7,3) * t369 + Ifges(7,6) * t371 + Ifges(7,5) * t372 - Ifges(6,6) * t206 - Ifges(6,4) * t61 - Ifges(6,2) * t62 + t92 * mrSges(6,1) + t288 / 0.2e1 - t14 * mrSges(6,3) + t378) * t244 - (mrSges(5,2) * t154 + Ifges(5,1) * t116 + Ifges(5,4) * t117 + Ifges(5,5) * qJDD(4)) * t243 + t78 * t278 + (-t165 - t398) * t270 + t44 * (mrSges(7,1) * t237 - mrSges(7,2) * t236) + t230 * t392 + t380 * mrSges(5,3) + (-m(6) * t49 - t381) * (qJD(5) * t69 + t220 * t83 - t224 * t229); -t243 * t106 - t160 * t107 + t152 * t132 + t153 * t133 - t344 * t78 + (-m(3) * qJ(2) - mrSges(3,3)) * t227 - (t16 - t54) * t121 - t239 * t230 + (mrSges(3,2) - t270) * qJDD(1) - (qJD(6) * t249 + t387 + t55) * t244 + m(6) * (t121 * t15 - t14 * t244 + t230 * t50 + t49 * t78) - m(5) * t380 + m(3) * t198 + m(7) * (-t12 * t121 + t230 * t251 - t231 * t244 - t44 * t78) + m(4) * t266 + (-m(6) * t124 - m(7) * t252 + t249 - t375 - t75) * qJD(1) + t385 * (m(3) - t283); t250 * qJD(6) - t344 * t245 + t239 * t263 - t150 * t132 + t151 * t133 + t223 * t17 + t219 * t18 - t227 * t270 + t269 + t274 + t302 + (g(1) * t226 + g(2) * t222) * t283 + (t105 * t251 + t2 * t219 + t3 * t223 - t245 * t44) * m(7) + (t245 * t49 - t263 * t50 + t92) * m(6) + (t102 * t151 - t103 * t150 + t154) * m(5) + (qJD(1) * t267 + t176) * m(4); -(-Ifges(5,2) * t151 + t109 + t146) * t150 / 0.2e1 + t385 * (mrSges(5,1) * t200 - mrSges(5,2) * t199) + (-Ifges(6,4) * t363 + t348 - t40 / 0.2e1 + t73 / 0.2e1 + t374) * t245 + (Ifges(6,4) * t364 + t373 + t349 - t74 / 0.2e1) * t263 + (t339 - t132) * t102 + t376 * (pkin(9) + t354) + (t338 + t133) * t103 - t344 * t51 + t228 + t108 * t360 + t54 * t353 + t55 * t354 + t195 * t16 - t169 * (mrSges(5,1) * t151 + mrSges(5,2) * t150) - qJD(4) * (Ifges(5,5) * t150 - Ifges(5,6) * t151) / 0.2e1 + Ifges(5,6) * t117 - t75 * t357 + Ifges(5,5) * t116 - t52 * t98 - t151 * (Ifges(5,1) * t150 - t333) / 0.2e1 + t67 * mrSges(5,1) - t25 * t70 - t24 * t71 - t66 * mrSges(5,2) + Ifges(5,3) * qJDD(4) + (-t124 * t357 + t49 * t51 - t50 * t52) * m(6) + (t12 * t195 - t20 * t24 - t21 * t25 - t44 * t51) * m(7) + (-t239 * t295 + t344 * t296 + (t14 * t220 + t15 * t224 + (-t220 * t49 + t224 * t50) * qJD(5)) * m(6) + (t220 * t44 + t224 * t251) * qJD(5) * m(7)) * pkin(4) + (m(6) * t356 - m(7) * (t271 - t356) - t258 + t388) * g(3) + (-m(7) * (t222 * t355 + t303) - (t289 - t341) * t222 + t389) * g(1) + ((-m(7) * (t261 - t355) + t289 + t386) * t226 + t390) * g(2); t228 + t73 * t362 - m(7) * (t20 * t26 + t21 * t27) - t27 * t70 - t26 * t71 + (t337 - t98) * t49 + (t104 + t74) * t364 + (t40 - t332) * t363 + (t388 - t393) * g(3) + ((-m(7) * t261 + t386) * t226 + t390) * g(2) + (-m(7) * t303 + mrSges(6,2) * t311 + t389) * g(1) + t382 * pkin(5) + (t336 + t381) * t50 + t374 * t245 + t373 * t263 + t376 * pkin(9); -t44 * (mrSges(7,1) * t94 + mrSges(7,2) * t93) + (Ifges(7,1) * t93 - t347) * t367 + t41 * t366 + (Ifges(7,5) * t93 - Ifges(7,6) * t94) * t365 - t20 * t70 + t21 * t71 - g(1) * (mrSges(7,1) * t140 - mrSges(7,2) * t141) - g(2) * (mrSges(7,1) * t142 + mrSges(7,2) * t143) + g(3) * t256 * t191 + (t20 * t93 + t21 * t94) * mrSges(7,3) + t288 + (-Ifges(7,2) * t94 + t42 + t91) * t368 + t378;];
tau  = t1;
