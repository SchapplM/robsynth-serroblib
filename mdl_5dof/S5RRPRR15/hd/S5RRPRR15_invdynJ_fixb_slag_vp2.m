% Calculate vector of inverse dynamics joint torques for
% S5RRPRR15
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR15_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:12
% EndTime: 2019-12-31 20:41:45
% DurationCPUTime: 18.28s
% Computational Cost: add. (5223->627), mult. (10814->861), div. (0->0), fcn. (6514->10), ass. (0->301)
t210 = sin(qJ(5));
t211 = sin(qJ(4));
t290 = qJD(4) + qJD(5);
t297 = qJD(5) * t210;
t301 = qJD(4) * t211;
t214 = cos(qJ(5));
t215 = cos(qJ(4));
t312 = t214 * t215;
t79 = -t210 * t301 - t211 * t297 + t290 * t312;
t212 = sin(qJ(2));
t195 = t212 * qJD(1);
t283 = t215 * t195;
t318 = t210 * t211;
t98 = -t195 * t318 + t214 * t283;
t416 = t79 + t98;
t242 = t210 * t215 + t214 * t211;
t80 = t290 * t242;
t231 = t242 * t212;
t99 = qJD(1) * t231;
t436 = t80 + t99;
t374 = pkin(2) + pkin(7);
t189 = pkin(6) * t195;
t151 = -pkin(3) * t195 - t189;
t399 = qJD(3) - t151;
t102 = -qJD(2) * t374 + t399;
t300 = qJD(4) * t215;
t216 = cos(qJ(2));
t295 = qJD(1) * qJD(2);
t155 = -t216 * qJDD(1) + t212 * t295;
t156 = qJDD(1) * t212 + t216 * t295;
t324 = qJDD(1) * pkin(1);
t226 = -qJ(3) * t156 - qJD(3) * t195 - t324;
t50 = t155 * t374 + t226;
t139 = t156 * pkin(6);
t274 = qJDD(3) + t139;
t70 = pkin(3) * t156 - qJDD(2) * t374 + t274;
t198 = t212 * qJ(3);
t273 = -pkin(1) - t198;
t93 = (-t216 * t374 + t273) * qJD(1);
t11 = t102 * t300 + t211 * t70 + t215 * t50 - t301 * t93;
t303 = qJD(2) * t215;
t305 = qJD(1) * t216;
t232 = t211 * t305 - t303;
t69 = qJD(4) * t232 - qJDD(2) * t211 + t155 * t215;
t10 = pkin(8) * t69 + t11;
t143 = -qJD(2) * t211 - t215 * t305;
t52 = t102 * t211 + t215 * t93;
t43 = pkin(8) * t143 + t52;
t337 = t210 * t43;
t180 = t195 + qJD(4);
t51 = t215 * t102 - t211 * t93;
t42 = pkin(8) * t232 + t51;
t38 = pkin(4) * t180 + t42;
t13 = t214 * t38 - t337;
t12 = -qJD(4) * t52 - t211 * t50 + t215 * t70;
t142 = qJDD(4) + t156;
t68 = qJD(4) * t143 + qJDD(2) * t215 + t155 * t211;
t9 = pkin(4) * t142 - pkin(8) * t68 + t12;
t2 = qJD(5) * t13 + t10 * t214 + t210 * t9;
t435 = t2 * mrSges(6,2);
t333 = t214 * t43;
t14 = t210 * t38 + t333;
t3 = -qJD(5) * t14 - t10 * t210 + t214 * t9;
t434 = t3 * mrSges(6,1);
t420 = -Ifges(4,4) + Ifges(3,5);
t419 = Ifges(4,5) - Ifges(3,6);
t317 = t211 * t212;
t240 = pkin(4) * t216 - pkin(8) * t317;
t351 = pkin(8) + t374;
t190 = pkin(2) * t195;
t325 = qJ(3) * t216;
t246 = pkin(7) * t212 - t325;
t108 = qJD(1) * t246 + t190;
t191 = pkin(6) * t305;
t152 = pkin(3) * t305 + t191;
t65 = -t108 * t211 + t215 * t152;
t433 = -qJD(1) * t240 + t351 * t301 - t65;
t159 = t351 * t215;
t66 = t215 * t108 + t211 * t152;
t432 = pkin(8) * t283 + qJD(4) * t159 + t66;
t209 = qJ(4) + qJ(5);
t196 = sin(t209);
t197 = cos(t209);
t259 = mrSges(5,1) * t211 + mrSges(5,2) * t215;
t363 = pkin(4) * t211;
t431 = -m(6) * t363 - t196 * mrSges(6,1) - t197 * mrSges(6,2) - t259;
t218 = -pkin(8) - pkin(7);
t430 = -m(6) * (-pkin(2) + t218) + mrSges(6,3) + m(5) * t374 + mrSges(5,3);
t268 = t214 * t143 + t210 * t232;
t133 = qJDD(5) + t142;
t21 = qJD(5) * t268 + t210 * t69 + t214 * t68;
t76 = t143 * t210 - t214 * t232;
t22 = -qJD(5) * t76 - t210 * t68 + t214 * t69;
t289 = Ifges(6,5) * t21 + Ifges(6,6) * t22 + Ifges(6,3) * t133;
t365 = Ifges(6,4) * t76;
t171 = t195 + t290;
t368 = -t171 / 0.2e1;
t376 = -t76 / 0.2e1;
t208 = qJD(2) * qJ(3);
t122 = t208 + t152;
t82 = -pkin(4) * t143 + t122;
t429 = t434 - t435 + t289 + (Ifges(6,5) * t268 - Ifges(6,6) * t76) * t368 + (t13 * t268 + t14 * t76) * mrSges(6,3) - t82 * (mrSges(6,1) * t76 + mrSges(6,2) * t268) + (Ifges(6,1) * t268 - t365) * t376;
t340 = Ifges(4,6) * t216;
t248 = -t212 * Ifges(4,2) - t340;
t428 = t14 * mrSges(6,2) + Ifges(4,4) * qJD(2) / 0.2e1 + qJD(1) * t248 / 0.2e1 - t13 * mrSges(6,1);
t213 = sin(qJ(1));
t427 = g(2) * t213;
t258 = t216 * mrSges(4,2) - t212 * mrSges(4,3);
t262 = mrSges(3,1) * t216 - mrSges(3,2) * t212;
t426 = t258 - t262;
t241 = -t312 + t318;
t425 = t436 * t13 - t14 * t416 - t2 * t242 + t241 * t3;
t424 = qJD(3) + t189;
t71 = Ifges(6,4) * t268;
t423 = -Ifges(6,2) * t76 + t71;
t347 = mrSges(5,3) * t143;
t88 = -mrSges(5,2) * t180 + t347;
t346 = mrSges(5,3) * t232;
t89 = mrSges(5,1) * t180 + t346;
t244 = -t211 * t89 + t215 * t88;
t47 = mrSges(5,1) * t142 - mrSges(5,3) * t68;
t48 = -mrSges(5,2) * t142 + mrSges(5,3) * t69;
t422 = t244 * qJD(4) + t211 * t48 + t215 * t47;
t386 = t21 / 0.2e1;
t385 = t22 / 0.2e1;
t372 = t133 / 0.2e1;
t158 = t351 * t211;
t84 = -t158 * t214 - t159 * t210;
t418 = -qJD(5) * t84 + t210 * t432 + t214 * t433;
t83 = t158 * t210 - t159 * t214;
t417 = qJD(5) * t83 + t210 * t433 - t214 * t432;
t387 = m(6) * pkin(4);
t410 = -mrSges(5,1) - t387;
t284 = mrSges(4,1) * t305;
t165 = -qJD(2) * mrSges(4,3) - t284;
t81 = -mrSges(5,1) * t143 - mrSges(5,2) * t232;
t409 = t81 - t165;
t361 = pkin(4) * t215;
t185 = pkin(3) + t361;
t408 = pkin(4) * t300 + t185 * t195 + t424;
t202 = t216 * pkin(2);
t307 = t202 + t198;
t264 = pkin(7) * t216 + t307;
t135 = -pkin(1) - t264;
t373 = pkin(3) + pkin(6);
t168 = t373 * t212;
t147 = t211 * t168;
t78 = t215 * t135 + t147;
t407 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t305 - t165;
t285 = mrSges(4,1) * t195;
t406 = -mrSges(3,3) * t195 - t285 + (mrSges(3,1) - mrSges(4,2)) * qJD(2);
t341 = Ifges(4,6) * t212;
t405 = t212 * (-Ifges(4,2) * t216 + t341) + t216 * (Ifges(4,3) * t212 - t340);
t404 = t212 * t419 + t216 * t420;
t138 = t155 * pkin(6);
t403 = -t138 * t216 + t139 * t212;
t100 = -qJDD(2) * qJ(3) - qJD(2) * qJD(3) + t138;
t107 = -qJDD(2) * pkin(2) + t274;
t402 = -t100 * t216 + t107 * t212;
t217 = cos(qJ(1));
t400 = g(1) * t217 + t427;
t296 = m(4) + m(5) + m(6);
t188 = Ifges(3,4) * t305;
t398 = Ifges(3,1) * t195 + Ifges(3,5) * qJD(2) - Ifges(5,5) * t232 + t76 * Ifges(6,5) + t143 * Ifges(5,6) + Ifges(6,6) * t268 + t180 * Ifges(5,3) + t171 * Ifges(6,3) + t188;
t397 = -mrSges(2,1) + t426;
t328 = t216 * mrSges(4,3);
t396 = -t328 + t431 * t216 + (m(4) * pkin(2) - mrSges(4,2) + t430) * t212;
t247 = -t216 * Ifges(4,3) - t341;
t338 = t232 * Ifges(5,4);
t58 = t143 * Ifges(5,2) + t180 * Ifges(5,6) - t338;
t136 = Ifges(5,4) * t143;
t59 = -Ifges(5,1) * t232 + t180 * Ifges(5,5) + t136;
t395 = Ifges(4,5) * qJD(2) + qJD(1) * t247 + t211 * t59 + t215 * t58;
t394 = t216 * t290;
t391 = -m(5) * pkin(3) - m(6) * t185 - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t389 = Ifges(6,4) * t386 + Ifges(6,2) * t385 + Ifges(6,6) * t372;
t388 = Ifges(6,1) * t386 + Ifges(6,4) * t385 + Ifges(6,5) * t372;
t384 = -t68 * Ifges(5,4) / 0.2e1 - t69 * Ifges(5,2) / 0.2e1 - t142 * Ifges(5,6) / 0.2e1;
t32 = Ifges(6,2) * t268 + Ifges(6,6) * t171 + t365;
t383 = -t32 / 0.2e1;
t33 = Ifges(6,1) * t76 + Ifges(6,5) * t171 + t71;
t382 = -t33 / 0.2e1;
t381 = t33 / 0.2e1;
t380 = t68 / 0.2e1;
t379 = t69 / 0.2e1;
t378 = -t268 / 0.2e1;
t377 = t268 / 0.2e1;
t375 = t76 / 0.2e1;
t371 = t142 / 0.2e1;
t369 = -t232 / 0.2e1;
t367 = t171 / 0.2e1;
t364 = pkin(4) * t232;
t360 = pkin(6) * t212;
t200 = t216 * pkin(6);
t314 = t212 * t217;
t103 = -t196 * t213 + t197 * t314;
t104 = t196 * t314 + t197 * t213;
t350 = t103 * mrSges(6,1) - t104 * mrSges(6,2);
t315 = t212 * t213;
t105 = t196 * t217 + t197 * t315;
t106 = -t196 * t315 + t197 * t217;
t349 = t105 * mrSges(6,1) + t106 * mrSges(6,2);
t348 = mrSges(6,1) * t197;
t345 = Ifges(3,4) * t212;
t344 = Ifges(3,4) * t216;
t343 = Ifges(5,4) * t211;
t342 = Ifges(5,4) * t215;
t339 = t11 * t211;
t332 = t215 * mrSges(5,3);
t327 = t216 * mrSges(6,3);
t316 = t211 * t216;
t313 = t213 * t215;
t311 = t215 * t216;
t310 = t215 * t217;
t309 = t216 * t217;
t308 = t216 * t218;
t169 = t216 * pkin(3) + t200;
t306 = t217 * pkin(1) + t213 * pkin(6);
t304 = qJD(2) * t212;
t302 = qJD(2) * t216;
t299 = qJD(4) * t216;
t288 = Ifges(5,5) * t68 + Ifges(5,6) * t69 + Ifges(5,3) * t142;
t282 = t211 * t299;
t276 = -t300 / 0.2e1;
t182 = qJ(3) + t363;
t272 = pkin(8) * t216 - t135;
t271 = -t295 / 0.2e1;
t113 = t156 * mrSges(4,1) + qJDD(2) * mrSges(4,2);
t154 = t373 * t302;
t267 = pkin(2) * t304 - qJD(3) * t212;
t90 = qJD(2) * t246 + t267;
t269 = t215 * t154 - t211 * t90;
t261 = mrSges(3,1) * t212 + mrSges(3,2) * t216;
t260 = mrSges(5,1) * t215 - mrSges(5,2) * t211;
t257 = Ifges(5,1) * t215 - t343;
t256 = Ifges(5,1) * t211 + t342;
t255 = t216 * Ifges(3,2) + t345;
t253 = -Ifges(5,2) * t211 + t342;
t252 = Ifges(5,2) * t215 + t343;
t250 = Ifges(5,5) * t215 - Ifges(5,6) * t211;
t249 = Ifges(5,5) * t211 + Ifges(5,6) * t215;
t148 = t215 * t168;
t56 = pkin(4) * t212 + t211 * t272 + t148;
t62 = -pkin(8) * t311 + t78;
t29 = -t210 * t62 + t214 * t56;
t30 = t210 * t56 + t214 * t62;
t245 = t51 * t211 - t52 * t215;
t157 = -qJD(2) * pkin(2) + t424;
t163 = -t191 - t208;
t243 = t157 * t216 + t163 * t212;
t239 = t273 - t202;
t237 = pkin(1) * t261;
t124 = -t211 * t213 + t212 * t310;
t126 = t211 * t217 + t212 * t313;
t123 = t239 * qJD(1);
t236 = t123 * (-mrSges(4,2) * t212 - t328);
t235 = t212 * (Ifges(3,1) * t216 - t345);
t35 = -t135 * t301 + t211 * t154 + t168 * t300 + t215 * t90;
t229 = t212 * t303 + t282;
t228 = t211 * t304 - t215 * t299;
t225 = Ifges(5,5) * t216 + t212 * t256;
t224 = Ifges(5,6) * t216 + t212 * t252;
t223 = Ifges(5,3) * t216 + t212 * t249;
t72 = -pkin(3) * t155 - t100;
t222 = -qJD(4) * t245 + t12 * t215 + t339;
t170 = t216 * t196 * mrSges(6,2);
t160 = -pkin(1) - t307;
t153 = t373 * t304;
t150 = -qJ(3) * t305 + t190;
t149 = t258 * qJD(1);
t134 = t260 * t216;
t127 = -t211 * t315 + t310;
t125 = t211 * t314 + t313;
t118 = Ifges(3,6) * qJD(2) + qJD(1) * t255;
t117 = pkin(4) * t311 + t169;
t114 = -qJ(3) * t302 + t267;
t112 = mrSges(4,1) * t155 - qJDD(2) * mrSges(4,3);
t110 = t242 * t216;
t109 = t241 * t216;
t85 = -pkin(4) * t282 + (-pkin(6) - t185) * t304;
t77 = -t135 * t211 + t148;
t67 = pkin(2) * t155 + t226;
t55 = mrSges(6,1) * t171 - mrSges(6,3) * t76;
t54 = -mrSges(6,2) * t171 + mrSges(6,3) * t268;
t45 = -t241 * t304 + t242 * t394;
t44 = qJD(2) * t231 + t241 * t394;
t39 = -mrSges(6,1) * t268 + mrSges(6,2) * t76;
t37 = -pkin(4) * t69 + t72;
t36 = -qJD(4) * t78 + t269;
t34 = -mrSges(5,1) * t69 + mrSges(5,2) * t68;
t28 = t68 * Ifges(5,1) + t69 * Ifges(5,4) + t142 * Ifges(5,5);
t26 = pkin(8) * t229 + t35;
t25 = t240 * qJD(2) + (t215 * t272 - t147) * qJD(4) + t269;
t18 = t214 * t42 - t337;
t17 = -t210 * t42 - t333;
t16 = -mrSges(6,2) * t133 + mrSges(6,3) * t22;
t15 = mrSges(6,1) * t133 - mrSges(6,3) * t21;
t8 = -mrSges(6,1) * t22 + mrSges(6,2) * t21;
t5 = -qJD(5) * t30 - t210 * t26 + t214 * t25;
t4 = qJD(5) * t29 + t210 * t25 + t214 * t26;
t1 = [(Ifges(6,1) * t44 + Ifges(6,4) * t45) * t375 + t52 * mrSges(5,3) * t229 + (-qJDD(2) * mrSges(3,2) - t112) * t200 + (Ifges(6,4) * t44 + Ifges(6,2) * t45) * t377 + m(4) * (t114 * t123 + t160 * t67 + (qJD(2) * t243 + t402) * pkin(6)) + (-t155 * t200 + t156 * t360 + t403) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t403) + t404 * qJD(2) ^ 2 / 0.2e1 + t405 * t271 + (-m(3) * t306 - t125 * mrSges(5,1) - t104 * mrSges(6,1) - t124 * mrSges(5,2) - t103 * mrSges(6,2) + (-m(5) * pkin(7) - mrSges(5,3)) * t309 - t296 * (pkin(2) * t309 + qJ(3) * t314 + t306) + t391 * t213 + (-m(6) * (pkin(4) * t317 - t308) - t327 + t397) * t217) * g(2) + t212 * t434 - t51 * mrSges(5,3) * t228 + (t212 * t420 - t216 * t419) * qJDD(2) / 0.2e1 + t156 * t344 / 0.2e1 + (-Ifges(3,4) * t155 + Ifges(3,5) * qJDD(2) + t288 + t289) * t212 / 0.2e1 + (-qJDD(2) * mrSges(3,1) + t113) * t360 + t262 * t324 + t402 * mrSges(4,1) + t156 * t212 * Ifges(3,1) + (t157 * mrSges(4,1) + t398 / 0.2e1 + Ifges(6,3) * t367 + Ifges(6,5) * t375 + Ifges(6,6) * t377 - t52 * mrSges(5,2) + t51 * mrSges(5,1) - t406 * pkin(6) - t428) * t302 + (t109 * t2 + t110 * t3 - t13 * t44 + t14 * t45) * mrSges(6,3) + t122 * (-mrSges(5,1) * t229 + mrSges(5,2) * t228) + m(6) * (t117 * t37 + t13 * t5 + t14 * t4 + t2 * t30 + t29 * t3 + t82 * t85) + (t216 * (-Ifges(3,2) * t212 + t344) + t235) * t295 / 0.2e1 + t58 * t282 / 0.2e1 - t28 * t316 / 0.2e1 + t12 * (t212 * mrSges(5,1) + mrSges(5,3) * t316) + t216 * t59 * t276 + (Ifges(5,6) * t212 - t216 * t252) * t379 + (Ifges(5,5) * t212 - t216 * t256) * t380 + t44 * t381 + t311 * t384 - t110 * t388 + t109 * t389 + (qJD(2) * t225 - t257 * t299) * t369 + (Ifges(5,3) * t212 - t216 * t249) * t371 + t11 * (-mrSges(5,2) * t212 - mrSges(5,3) * t311) + t143 * (qJD(2) * t224 - t253 * t299) / 0.2e1 + t180 * (qJD(2) * t223 - t250 * t299) / 0.2e1 - t237 * t295 + (Ifges(6,5) * t44 + Ifges(6,6) * t45) * t367 + (-Ifges(6,4) * t110 + Ifges(6,2) * t109 + Ifges(6,6) * t212) * t385 + (-Ifges(6,1) * t110 + Ifges(6,4) * t109 + Ifges(6,5) * t212) * t386 + (-Ifges(6,5) * t110 + Ifges(6,6) * t109 + Ifges(6,3) * t212) * t372 + t37 * (-mrSges(6,1) * t109 - mrSges(6,2) * t110) + m(5) * (t11 * t78 + t12 * t77 - t122 * t153 + t169 * t72 + t35 * t52 + t36 * t51) + (t163 * mrSges(4,1) + t395 / 0.2e1 - t118 / 0.2e1 - t407 * pkin(6)) * t304 + t29 * t15 + t30 * t16 + (-t127 * mrSges(5,1) - t106 * mrSges(6,1) + t126 * mrSges(5,2) + t105 * mrSges(6,2) + (-m(6) * (-t182 * t212 - pkin(1)) - m(5) * t273 - m(4) * t239 + m(3) * pkin(1) + t430 * t216 - t397) * t213 + ((-m(3) - t296) * pkin(6) + t391) * t217) * g(1) + t45 * t32 / 0.2e1 + t4 * t54 + t5 * t55 + t77 * t47 + t78 * t48 + t82 * (-mrSges(6,1) * t45 + mrSges(6,2) * t44) + t85 * t39 + t35 * t88 + t36 * t89 - t212 * t435 + t117 * t8 + t72 * t134 + t155 * t247 / 0.2e1 - t156 * t248 / 0.2e1 + qJD(2) * t236 - t155 * t255 / 0.2e1 + t67 * t258 + t114 * t149 - t153 * t81 - pkin(1) * (mrSges(3,1) * t155 + mrSges(3,2) * t156) + t160 * (-mrSges(4,2) * t155 - mrSges(4,3) * t156) + t169 * t34 - t212 * (Ifges(4,4) * qJDD(2) - Ifges(4,2) * t156 + Ifges(4,6) * t155) / 0.2e1 + t216 * (Ifges(3,4) * t156 - Ifges(3,2) * t155 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t216 * (Ifges(4,5) * qJDD(2) - Ifges(4,6) * t156 + Ifges(4,3) * t155) / 0.2e1 + Ifges(2,3) * qJDD(1); (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + t404 * t271 + t406 * t191 + t407 * t189 + t408 * t39 + (-qJ(3) * t296 * t309 + t217 * t396) * g(1) + t400 * t261 + t425 * mrSges(6,3) + (-t236 - t51 * (t216 * mrSges(5,1) - mrSges(5,3) * t317) - t52 * (-mrSges(5,2) * t216 + t212 * t332) - m(4) * t243 * pkin(6)) * qJD(1) + (t72 * qJ(3) + t122 * t399 - t51 * t65 - t52 * t66) * m(5) + t409 * qJD(3) + (t34 - t112) * qJ(3) + t416 * t383 + t417 * t54 + (t418 * t13 + t417 * t14 + t182 * t37 + t2 * t84 + t3 * t83 + t408 * t82) * m(6) + t418 * t55 + t419 * t155 + t420 * t156 + (mrSges(6,1) * t416 - mrSges(6,2) * t436) * t82 - t395 * t195 / 0.2e1 - (-Ifges(3,2) * t195 + t188 + t398) * t305 / 0.2e1 + (-t300 * t52 + t301 * t51 - t339) * mrSges(5,3) + (-m(5) * t222 - t422) * t374 + (Ifges(6,4) * t99 + Ifges(6,2) * t98) * t378 - (t143 * t252 + t180 * t249 - t232 * t256) * qJD(4) / 0.2e1 - (t143 * t224 + t180 * t223 - t225 * t232) * qJD(1) / 0.2e1 - t242 * t389 + (-Ifges(6,4) * t241 - Ifges(6,2) * t242) * t385 + (-Ifges(6,1) * t241 - Ifges(6,4) * t242) * t386 + (-Ifges(6,5) * t241 - Ifges(6,6) * t242) * t372 + t37 * (mrSges(6,1) * t242 - mrSges(6,2) * t241) - t241 * t388 + (-Ifges(6,5) * t80 - Ifges(6,6) * t79) * t367 + (-Ifges(6,1) * t80 - Ifges(6,4) * t79) * t375 + (-Ifges(6,4) * t80 - Ifges(6,2) * t79) * t377 + (t405 / 0.2e1 - t235 / 0.2e1 + t237) * qJD(1) ^ 2 + (Ifges(6,1) * t99 + Ifges(6,4) * t98) * t376 + (Ifges(6,5) * t376 + Ifges(6,6) * t378 + Ifges(6,3) * t368 + t428) * t305 - t157 * t284 - t163 * t285 + t253 * t379 + t257 * t380 - t80 * t381 + t99 * t382 + t211 * t384 + t250 * t371 + (-pkin(2) * t107 - qJ(3) * t100 - qJD(3) * t163 - t123 * t150) * m(4) - t59 * t301 / 0.2e1 + (Ifges(6,5) * t99 + Ifges(6,6) * t98) * t368 + t118 * t195 / 0.2e1 + t180 * t260 * t122 - t12 * t332 + (-t296 * t325 + t396) * t427 + (-m(6) * (t307 - t308) - t327 - m(4) * t307 - m(5) * t264 - t216 * mrSges(5,3) + t431 * t212 + t426) * g(3) + t83 * t15 + t84 * t16 - t66 * t88 - t65 * t89 - t100 * mrSges(4,3) + t107 * mrSges(4,2) - pkin(2) * t113 + t58 * t276 + t138 * mrSges(3,2) - t139 * mrSges(3,1) - t150 * t149 + t72 * t259 - t151 * t81 + t182 * t8 + t215 * t28 / 0.2e1; t242 * t16 - t241 * t15 - t436 * t55 + t416 * t54 + (-t39 - t409) * qJD(2) + t296 * t216 * g(3) + ((t149 + t244) * qJD(1) - t400 * t296) * t212 + t113 + (-qJD(2) * t82 - t425) * m(6) + (-qJD(2) * t122 - t195 * t245 + t222) * m(5) + (qJD(2) * t163 + t123 * t195 + t107) * m(4) + t422; t423 * t378 + (mrSges(5,2) * t125 + t124 * t410 - t350) * g(1) + (-mrSges(5,2) * t127 + t126 * t410 - t349) * g(2) - t76 * t383 + t232 * (Ifges(5,1) * t143 + t338) / 0.2e1 - (Ifges(5,2) * t232 + t136 + t59) * t143 / 0.2e1 - t122 * (-mrSges(5,1) * t232 + mrSges(5,2) * t143) - t180 * (Ifges(5,5) * t143 + Ifges(5,6) * t232) / 0.2e1 - t11 * mrSges(5,2) + t12 * mrSges(5,1) + t288 + (t2 * t210 + t214 * t3 + (-t13 * t210 + t14 * t214) * qJD(5)) * t387 + t58 * t369 + (-t346 + t89) * t52 - m(6) * (t13 * t17 + t14 * t18 - t364 * t82) + t429 + t39 * t364 + t268 * t382 + (t347 - t88) * t51 + (-t170 - (-m(6) * t361 - t348) * t216 + t134) * g(3) - t18 * t54 - t17 * t55 + (t16 * t210 - t297 * t55 + (qJD(5) * t54 + t15) * t214) * pkin(4); t32 * t375 - t13 * t54 + t14 * t55 - g(1) * t350 - g(2) * t349 - g(3) * (-t216 * t348 + t170) + (t33 + t423) * t378 + t429;];
tau = t1;
