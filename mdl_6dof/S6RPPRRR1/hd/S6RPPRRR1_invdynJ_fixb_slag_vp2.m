% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:02
% EndTime: 2019-03-09 02:18:17
% DurationCPUTime: 10.67s
% Computational Cost: add. (11094->569), mult. (25148->756), div. (0->0), fcn. (18969->18), ass. (0->266)
t204 = pkin(11) + qJ(4);
t197 = qJ(5) + t204;
t182 = sin(t197);
t183 = cos(t197);
t212 = sin(qJ(6));
t327 = mrSges(7,2) * t212;
t399 = -t182 * t327 + t183 * (-m(7) * pkin(9) - mrSges(7,3));
t208 = sin(pkin(10));
t179 = pkin(1) * t208 + qJ(3);
t161 = qJD(1) * qJD(3) + qJDD(1) * t179;
t377 = t183 * pkin(5) + t182 * pkin(9);
t398 = m(7) * t377;
t205 = qJD(4) + qJD(5);
t207 = sin(pkin(11));
t214 = sin(qJ(4));
t209 = cos(pkin(11));
t218 = cos(qJ(4));
t292 = t209 * t218;
t163 = -t207 * t214 + t292;
t153 = t163 * qJD(1);
t164 = t207 * t218 + t209 * t214;
t154 = t164 * qJD(1);
t213 = sin(qJ(5));
t217 = cos(qJ(5));
t241 = t153 * t213 + t217 * t154;
t216 = cos(qJ(6));
t98 = t205 * t216 - t212 * t241;
t99 = t205 * t212 + t216 * t241;
t306 = -mrSges(6,1) * t205 - mrSges(7,1) * t98 + mrSges(7,2) * t99 + mrSges(6,3) * t241;
t397 = -t183 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t182;
t396 = t207 ^ 2 + t209 ^ 2;
t200 = qJDD(4) + qJDD(5);
t155 = t163 * qJD(4);
t116 = qJD(1) * t155 + qJDD(1) * t164;
t156 = t164 * qJD(4);
t117 = -qJD(1) * t156 + qJDD(1) * t163;
t263 = t217 * t153 - t154 * t213;
t63 = qJD(5) * t263 + t116 * t217 + t117 * t213;
t34 = qJD(6) * t98 + t200 * t212 + t216 * t63;
t64 = -qJD(5) * t241 - t116 * t213 + t117 * t217;
t62 = qJDD(6) - t64;
t17 = mrSges(7,1) * t62 - mrSges(7,3) * t34;
t35 = -qJD(6) * t99 + t200 * t216 - t212 * t63;
t18 = -mrSges(7,2) * t62 + mrSges(7,3) * t35;
t248 = -t212 * t17 + t216 * t18;
t286 = qJD(6) * t216;
t287 = qJD(6) * t212;
t103 = qJD(6) - t263;
t71 = -mrSges(7,2) * t103 + mrSges(7,3) * t98;
t72 = mrSges(7,1) * t103 - mrSges(7,3) * t99;
t395 = -t72 * t286 - t71 * t287 + t248;
t330 = mrSges(7,1) * t216;
t394 = t327 - t330;
t210 = cos(pkin(10));
t185 = -pkin(1) * t210 - pkin(2);
t198 = t209 * pkin(3);
t376 = t185 - t198;
t151 = qJD(1) * t376 + qJD(3);
t115 = -pkin(4) * t153 + t151;
t315 = t205 * Ifges(6,6);
t321 = Ifges(6,4) * t241;
t173 = t179 * qJD(1);
t192 = t209 * qJD(2);
t318 = pkin(7) * qJD(1);
t128 = t192 + (-t173 - t318) * t207;
t140 = t207 * qJD(2) + t209 * t173;
t129 = t209 * t318 + t140;
t89 = t128 * t214 + t129 * t218;
t82 = pkin(8) * t153 + t89;
t307 = t217 * t82;
t302 = t129 * t214;
t88 = t218 * t128 - t302;
t81 = -pkin(8) * t154 + t88;
t79 = qJD(4) * pkin(4) + t81;
t45 = t213 * t79 + t307;
t336 = t45 * mrSges(6,3);
t40 = pkin(9) * t205 + t45;
t65 = -pkin(5) * t263 - pkin(9) * t241 + t115;
t21 = t212 * t65 + t216 * t40;
t382 = t21 * mrSges(7,2);
t20 = -t212 * t40 + t216 * t65;
t383 = t20 * mrSges(7,1);
t317 = t103 * Ifges(7,3);
t333 = t99 * Ifges(7,5);
t335 = t98 * Ifges(7,6);
t50 = t317 + t333 + t335;
t379 = t263 * Ifges(6,2);
t73 = t315 + t321 + t379;
t393 = -t115 * mrSges(6,1) + t336 + t73 / 0.2e1 - t50 / 0.2e1 + t382 - t383 + t321 / 0.2e1 + t315 / 0.2e1;
t102 = Ifges(6,4) * t263;
t253 = mrSges(7,1) * t212 + mrSges(7,2) * t216;
t311 = t213 * t82;
t44 = t217 * t79 - t311;
t39 = -pkin(5) * t205 - t44;
t234 = t39 * t253;
t97 = Ifges(7,4) * t98;
t52 = Ifges(7,1) * t99 + Ifges(7,5) * t103 + t97;
t309 = t216 * t52;
t316 = t205 * Ifges(6,5);
t337 = t44 * mrSges(6,3);
t348 = t212 / 0.2e1;
t334 = t99 * Ifges(7,4);
t51 = t98 * Ifges(7,2) + t103 * Ifges(7,6) + t334;
t380 = t241 * Ifges(6,1);
t74 = t102 + t316 + t380;
t392 = -t115 * mrSges(6,2) - t234 - t309 / 0.2e1 + t51 * t348 + t337 - t74 / 0.2e1 - t316 / 0.2e1 - t102 / 0.2e1;
t389 = m(4) + m(5);
t388 = -m(7) - m(6);
t16 = -mrSges(7,1) * t35 + mrSges(7,2) * t34;
t56 = mrSges(6,1) * t200 - mrSges(6,3) * t63;
t381 = t16 - t56;
t100 = -mrSges(6,2) * t205 + mrSges(6,3) * t263;
t245 = -t212 * t72 + t216 * t71;
t237 = -t100 - t245;
t331 = pkin(7) + t179;
t157 = t331 * t207;
t158 = t331 * t209;
t113 = -t214 * t157 + t218 * t158;
t375 = t183 * t394 + t397;
t119 = t163 * t213 + t164 * t217;
t240 = t217 * t163 - t164 * t213;
t83 = qJD(5) * t240 + t155 * t217 - t156 * t213;
t236 = t119 * t286 + t212 * t83;
t372 = -t20 * t286 - t21 * t287;
t190 = t209 * qJDD(2);
t130 = -t161 * t207 + t190;
t131 = t207 * qJDD(2) + t209 * t161;
t371 = -t130 * t207 + t131 * t209;
t369 = m(3) + t389;
t124 = t190 + (-pkin(7) * qJDD(1) - t161) * t207;
t283 = qJDD(1) * t209;
t125 = pkin(7) * t283 + t131;
t55 = -qJD(4) * t89 + t218 * t124 - t125 * t214;
t42 = qJDD(4) * pkin(4) - pkin(8) * t116 + t55;
t290 = qJD(4) * t218;
t54 = -qJD(4) * t302 + t214 * t124 + t218 * t125 + t128 * t290;
t43 = pkin(8) * t117 + t54;
t11 = -qJD(5) * t45 - t213 * t43 + t217 * t42;
t149 = qJDD(1) * t376 + qJDD(3);
t94 = -pkin(4) * t117 + t149;
t19 = -pkin(5) * t64 - pkin(9) * t63 + t94;
t288 = qJD(5) * t217;
t289 = qJD(5) * t213;
t10 = t213 * t42 + t217 * t43 + t79 * t288 - t289 * t82;
t7 = pkin(9) * t200 + t10;
t2 = qJD(6) * t20 + t19 * t212 + t216 * t7;
t3 = -qJD(6) * t21 + t19 * t216 - t212 * t7;
t367 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t77 = pkin(5) * t241 - pkin(9) * t263;
t184 = t198 + pkin(2);
t193 = sin(t204);
t195 = cos(t204);
t256 = mrSges(5,1) * t195 - mrSges(5,2) * t193;
t257 = -mrSges(4,1) * t209 + mrSges(4,2) * t207;
t366 = m(4) * pkin(2) + m(5) * t184 + mrSges(3,1) + t256 - t257 - t397;
t211 = -pkin(7) - qJ(3);
t364 = -m(4) * qJ(3) + m(5) * t211 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t362 = m(7) * pkin(5);
t361 = t34 / 0.2e1;
t360 = t35 / 0.2e1;
t357 = t62 / 0.2e1;
t354 = -t98 / 0.2e1;
t353 = -t99 / 0.2e1;
t352 = t99 / 0.2e1;
t351 = -t103 / 0.2e1;
t349 = t154 / 0.2e1;
t215 = sin(qJ(1));
t347 = pkin(1) * t215;
t346 = pkin(4) * t154;
t345 = pkin(4) * t156;
t344 = pkin(4) * t193;
t180 = pkin(4) * t195;
t343 = pkin(4) * t213;
t342 = pkin(4) * t217;
t339 = t2 * t216;
t219 = cos(qJ(1));
t199 = t219 * pkin(1);
t338 = t3 * t212;
t328 = mrSges(6,2) * t183;
t326 = mrSges(5,3) * t153;
t325 = mrSges(5,3) * t154;
t324 = mrSges(7,3) * t212;
t323 = mrSges(7,3) * t216;
t322 = Ifges(5,4) * t154;
t320 = Ifges(7,4) * t212;
t319 = Ifges(7,4) * t216;
t304 = t119 * t212;
t303 = t119 * t216;
t206 = qJ(1) + pkin(10);
t194 = sin(t206);
t296 = t194 * t212;
t295 = t194 * t216;
t196 = cos(t206);
t294 = t196 * t212;
t293 = t196 * t216;
t284 = qJDD(1) * t207;
t282 = Ifges(7,5) * t34 + Ifges(7,6) * t35 + Ifges(7,3) * t62;
t278 = -t388 + t389;
t271 = t309 / 0.2e1;
t267 = -t64 * mrSges(6,1) + t63 * mrSges(6,2);
t266 = -t287 / 0.2e1;
t265 = -t117 * mrSges(5,1) + t116 * mrSges(5,2);
t112 = -t218 * t157 - t158 * t214;
t262 = t399 * t194;
t261 = t399 * t196;
t258 = -mrSges(4,1) * t283 + mrSges(4,2) * t284;
t252 = Ifges(7,1) * t216 - t320;
t251 = -Ifges(7,2) * t212 + t319;
t250 = Ifges(7,5) * t216 - Ifges(7,6) * t212;
t247 = t20 * t216 + t21 * t212;
t246 = -t20 * t212 + t21 * t216;
t95 = -pkin(8) * t164 + t112;
t96 = pkin(8) * t163 + t113;
t67 = t213 * t95 + t217 * t96;
t127 = -pkin(4) * t163 + t376;
t70 = -pkin(5) * t240 - pkin(9) * t119 + t127;
t29 = t212 * t70 + t216 * t67;
t28 = -t212 * t67 + t216 * t70;
t66 = t213 * t96 - t217 * t95;
t242 = -(-t173 * t207 + t192) * t207 + t140 * t209;
t235 = t119 * t287 - t216 * t83;
t233 = t98 * t251;
t232 = t99 * t252;
t231 = t103 * t250;
t227 = -qJD(6) * t247 - t338;
t90 = -t157 * t290 + qJD(3) * t292 + (-qJD(3) * t207 - qJD(4) * t158) * t214;
t226 = m(7) * (-pkin(5) * t182 - t344) - t182 * t330;
t225 = t328 + (mrSges(6,1) + t330 + t362) * t182;
t224 = t227 + t339;
t91 = -t164 * qJD(3) - qJD(4) * t113;
t14 = t34 * Ifges(7,4) + t35 * Ifges(7,2) + t62 * Ifges(7,6);
t15 = t34 * Ifges(7,1) + t35 * Ifges(7,4) + t62 * Ifges(7,5);
t8 = -pkin(5) * t200 - t11;
t223 = -t10 * mrSges(6,2) + t2 * t323 + t15 * t348 + t216 * t14 / 0.2e1 + Ifges(6,3) * t200 + (Ifges(7,1) * t212 + t319) * t361 + (Ifges(7,2) * t216 + t320) * t360 + t51 * t266 + (Ifges(7,5) * t212 + Ifges(7,6) * t216) * t357 + t8 * t394 + Ifges(6,6) * t64 + Ifges(6,5) * t63 + t11 * mrSges(6,1) + (t234 + t271) * qJD(6) + (t233 + t232 + t231) * qJD(6) / 0.2e1;
t222 = -pkin(8) * t155 + t91;
t201 = -pkin(8) + t211;
t187 = -pkin(5) - t342;
t170 = qJDD(1) * t185 + qJDD(3);
t167 = t180 + t184;
t148 = Ifges(5,4) * t153;
t137 = qJD(4) * mrSges(5,1) - t325;
t136 = -qJD(4) * mrSges(5,2) + t326;
t135 = t183 * t293 + t296;
t134 = -t183 * t294 + t295;
t133 = -t183 * t295 + t294;
t132 = t183 * t296 + t293;
t107 = t154 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t148;
t106 = t153 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t322;
t105 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t117;
t104 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t116;
t85 = -pkin(8) * t156 + t90;
t84 = qJD(5) * t119 + t155 * t213 + t217 * t156;
t76 = -mrSges(6,1) * t263 + mrSges(6,2) * t241;
t69 = t346 + t77;
t57 = -mrSges(6,2) * t200 + mrSges(6,3) * t64;
t47 = t217 * t81 - t311;
t46 = t213 * t81 + t307;
t37 = pkin(5) * t84 - pkin(9) * t83 + t345;
t27 = t212 * t77 + t216 * t44;
t26 = -t212 * t44 + t216 * t77;
t25 = t212 * t69 + t216 * t47;
t24 = -t212 * t47 + t216 * t69;
t22 = -qJD(5) * t66 + t213 * t222 + t217 * t85;
t5 = -qJD(6) * t29 - t212 * t22 + t216 * t37;
t4 = qJD(6) * t28 + t212 * t37 + t216 * t22;
t1 = [(-m(6) * t44 + m(7) * t39 + t306) * (qJD(5) * t67 + t213 * t85 - t217 * t222) - t236 * t51 / 0.2e1 + m(4) * (t242 * qJD(3) + t170 * t185 + t179 * t371) + m(7) * (t2 * t29 + t20 * t5 + t21 * t4 + t28 * t3) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t210 - 0.2e1 * mrSges(3,2) * t208 + m(3) * (t208 ^ 2 + t210 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) - (-Ifges(6,4) * t63 - Ifges(6,2) * t64 - Ifges(6,6) * t200 - t10 * mrSges(6,3) + Ifges(7,3) * t357 + Ifges(7,6) * t360 + Ifges(7,5) * t361 + t94 * mrSges(6,1) + t282 / 0.2e1 + t367) * t240 + m(5) * (t112 * t55 + t113 * t54 + t149 * t376 + t88 * t91 + t89 * t90) + t376 * t265 + (mrSges(2,1) * t215 - t133 * mrSges(7,1) + mrSges(2,2) * t219 - t132 * mrSges(7,2) + t388 * (-t196 * t201 - t347) + t369 * t347 + t364 * t196 + (-m(7) * (-t167 - t377) + m(6) * t167 + t366) * t194) * g(1) + (t94 * mrSges(6,2) - t11 * mrSges(6,3) + Ifges(6,1) * t63 + Ifges(6,4) * t64 + Ifges(6,5) * t200 + t250 * t357 + t251 * t360 + t252 * t361 + t253 * t8 + t266 * t52) * t119 + t241 * (Ifges(6,1) * t83 - Ifges(6,4) * t84) / 0.2e1 + (-t155 * t88 - t156 * t89) * mrSges(5,3) + (Ifges(5,1) * t155 - Ifges(5,4) * t156) * t349 + qJD(4) * (Ifges(5,5) * t155 - Ifges(5,6) * t156) / 0.2e1 + t151 * (mrSges(5,1) * t156 + mrSges(5,2) * t155) + t153 * (Ifges(5,4) * t155 - Ifges(5,2) * t156) / 0.2e1 + m(6) * (t10 * t67 + t115 * t345 + t127 * t94 + t22 * t45) + t263 * (Ifges(6,4) * t83 - Ifges(6,2) * t84) / 0.2e1 - t84 * t382 + (-m(6) * t11 + m(7) * t8 + t381) * t66 + (t161 * t396 + t371) * mrSges(4,3) + (mrSges(5,2) * t149 - mrSges(5,3) * t55 + Ifges(5,1) * t116 + Ifges(5,4) * t117 + Ifges(5,5) * qJDD(4)) * t164 + (-t2 * t304 + t20 * t235 - t21 * t236 - t3 * t303) * mrSges(7,3) + t84 * t383 + (Ifges(4,4) * t207 + Ifges(4,2) * t209) * t283 + (Ifges(4,1) * t207 + Ifges(4,4) * t209) * t284 + t83 * t271 + (-mrSges(5,1) * t149 + mrSges(5,3) * t54 + Ifges(5,4) * t116 + Ifges(5,2) * t117 + Ifges(5,6) * qJDD(4)) * t163 + t76 * t345 - t84 * t336 - t83 * t337 + t205 * (Ifges(6,5) * t83 - Ifges(6,6) * t84) / 0.2e1 + t15 * t303 / 0.2e1 - t14 * t304 / 0.2e1 + (-mrSges(2,1) * t219 - t135 * mrSges(7,1) + mrSges(2,2) * t215 - t134 * mrSges(7,2) + t388 * (t196 * t167 - t194 * t201 + t199) - t369 * t199 + t364 * t194 + (-t366 - t398) * t196) * g(2) + t39 * (mrSges(7,1) * t236 - mrSges(7,2) * t235) + t98 * (-Ifges(7,4) * t235 - Ifges(7,2) * t236 + Ifges(7,6) * t84) / 0.2e1 + t103 * (-Ifges(7,5) * t235 - Ifges(7,6) * t236 + Ifges(7,3) * t84) / 0.2e1 + (-Ifges(7,1) * t235 - Ifges(7,4) * t236 + Ifges(7,5) * t84) * t352 + t170 * t257 + t185 * t258 + t28 * t17 + t29 * t18 + t67 * t57 + t4 * t71 + t5 * t72 + t83 * t74 / 0.2e1 + t84 * t50 / 0.2e1 - t84 * t73 / 0.2e1 + t22 * t100 + t127 * t267 + t112 * t104 + t113 * t105 + t115 * (mrSges(6,1) * t84 + mrSges(6,2) * t83) + t90 * t136 + t91 * t137 + t155 * t107 / 0.2e1 - t156 * t106 / 0.2e1; m(3) * qJDD(2) + t163 * t104 + t164 * t105 + t155 * t136 - t156 * t137 + t306 * t84 - t381 * t240 - t237 * t83 + (t57 + (-t212 * t71 - t216 * t72) * qJD(6) + t248) * t119 + (-m(3) - t278) * g(3) + m(6) * (t10 * t119 + t11 * t240 - t44 * t84 + t45 * t83) + m(5) * (t155 * t89 - t156 * t88 + t163 * t55 + t164 * t54) + m(4) * (t130 * t209 + t131 * t207) + m(7) * (t119 * t224 - t240 * t8 + t246 * t83 + t39 * t84); -t396 * qJD(1) ^ 2 * mrSges(4,3) + t267 + t265 + t258 - t306 * t241 + t216 * t17 + t212 * t18 + t237 * t263 + t245 * qJD(6) - t153 * t136 + t154 * t137 + (-g(1) * t194 + g(2) * t196) * t278 + (t103 * t246 + t2 * t212 + t216 * t3 - t241 * t39) * m(7) + (t241 * t44 - t263 * t45 + t94) * m(6) + (-t153 * t89 + t154 * t88 + t149) * m(5) + (-qJD(1) * t242 + t170) * m(4); t372 * mrSges(7,3) + (t20 * t323 + t21 * t324 + t250 * t351 + t252 * t353 + t251 * t354 - t380 / 0.2e1 + t392) * t263 + (Ifges(7,3) * t351 + Ifges(7,5) * t353 + Ifges(7,6) * t354 + t379 / 0.2e1 + t393) * t241 + (m(7) * t224 + t395) * (pkin(9) + t343) - t306 * (-pkin(4) * t289 + t46) + (t325 + t137) * t89 + (m(6) * t344 + mrSges(5,1) * t193 + mrSges(6,1) * t182 + mrSges(5,2) * t195 + t328) * (g(1) * t196 + g(2) * t194) + (-m(6) * t180 - t256 - m(7) * (t180 + t377) + t375) * g(3) - (-Ifges(5,2) * t154 + t107 + t148) * t153 / 0.2e1 + (t326 - t136) * t88 - t237 * pkin(4) * t288 + (t187 * t8 + (t213 * t39 + t217 * t246) * qJD(5) * pkin(4) - t20 * t24 - t21 * t25 - t39 * t46) * m(7) + t56 * t342 + t57 * t343 + t223 + ((t10 * t213 + t11 * t217 + (-t213 * t44 + t217 * t45) * qJD(5)) * pkin(4) - t115 * t346 + t44 * t46 - t45 * t47) * m(6) - t76 * t346 - t3 * t324 - t154 * (Ifges(5,1) * t153 - t322) / 0.2e1 + Ifges(5,3) * qJDD(4) + t187 * t16 + t106 * t349 - g(1) * (t196 * t226 - t261) - g(2) * (t194 * t226 - t262) - t54 * mrSges(5,2) + t55 * mrSges(5,1) - t25 * t71 - t24 * t72 - t47 * t100 + Ifges(5,5) * t116 + Ifges(5,6) * t117 - qJD(4) * (Ifges(5,5) * t153 - Ifges(5,6) * t154) / 0.2e1 - t151 * (mrSges(5,1) * t154 + mrSges(5,2) * t153); -m(7) * (t20 * t26 + t21 * t27 + t39 * t45) + (-t335 / 0.2e1 - t333 / 0.2e1 - t317 / 0.2e1 + t393) * t241 - t8 * t362 + (-t233 / 0.2e1 - t232 / 0.2e1 - t231 / 0.2e1 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t241 + t247 * mrSges(7,3) + t392) * t263 + t227 * mrSges(7,3) + (t196 * t225 + t261) * g(1) + (t194 * t225 + t262) * g(2) + t223 - t306 * t45 + (t375 - t398) * g(3) + (m(7) * (-t338 + t339 + t372) + t395) * pkin(9) - pkin(5) * t16 - t27 * t71 - t26 * t72 - t44 * t100; -t39 * (mrSges(7,1) * t99 + mrSges(7,2) * t98) + (Ifges(7,1) * t98 - t334) * t353 + t51 * t352 + (Ifges(7,5) * t98 - Ifges(7,6) * t99) * t351 - t20 * t71 + t21 * t72 - g(1) * (mrSges(7,1) * t134 - mrSges(7,2) * t135) - g(2) * (-mrSges(7,1) * t132 + mrSges(7,2) * t133) + g(3) * t253 * t182 + (t20 * t98 + t21 * t99) * mrSges(7,3) + t282 + (-Ifges(7,2) * t99 + t52 + t97) * t354 + t367;];
tau  = t1;
