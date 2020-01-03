% Calculate vector of inverse dynamics joint torques for
% S5RRPPP1
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:21
% EndTime: 2019-12-31 19:23:58
% DurationCPUTime: 21.91s
% Computational Cost: add. (3750->584), mult. (10065->773), div. (0->0), fcn. (7356->8), ass. (0->271)
t364 = Ifges(4,5) + Ifges(6,5);
t400 = -Ifges(5,4) + t364;
t365 = Ifges(6,4) + Ifges(5,5);
t399 = -Ifges(4,6) + t365;
t219 = sin(qJ(2));
t287 = qJD(1) * qJD(2);
t264 = t219 * t287;
t221 = cos(qJ(2));
t286 = qJDD(1) * t221;
t180 = -t264 + t286;
t181 = qJDD(1) * t219 + t221 * t287;
t215 = sin(pkin(8));
t217 = cos(pkin(8));
t216 = sin(pkin(5));
t285 = qJDD(2) * t216;
t218 = cos(pkin(5));
t308 = t217 * t218;
t91 = -t180 * t308 + t181 * t215 - t217 * t285;
t343 = -t91 / 0.2e1;
t243 = t180 * t218 + t285;
t92 = t181 * t217 + t215 * t243;
t341 = -t92 / 0.2e1;
t340 = t92 / 0.2e1;
t133 = qJDD(2) * t218 - t180 * t216;
t398 = -t133 / 0.2e1;
t335 = t133 / 0.2e1;
t367 = Ifges(4,1) + Ifges(6,3);
t366 = Ifges(4,4) - Ifges(6,6);
t363 = Ifges(6,2) + Ifges(5,3);
t362 = -Ifges(5,6) + Ifges(6,6);
t351 = Ifges(5,1) + Ifges(6,1) + Ifges(4,3);
t369 = -mrSges(5,2) + mrSges(4,1);
t393 = m(6) * qJ(5) + mrSges(6,3) + t369;
t368 = -mrSges(4,3) - mrSges(5,1);
t371 = m(6) * pkin(4) + mrSges(6,1) - t368;
t294 = qJD(1) * t221;
t171 = -t218 * qJD(2) + t216 * t294;
t208 = pkin(7) * t286;
t306 = t218 * t221;
t269 = qJD(3) * t306;
t292 = qJD(2) * t219;
t281 = pkin(7) * t292;
t290 = qJD(3) * t216;
t75 = qJD(2) * t290 + t208 + (t269 - t281) * qJD(1) + t243 * qJ(3);
t289 = qJD(3) * t219;
t235 = -qJ(3) * t181 - qJD(1) * t289;
t320 = qJDD(1) * pkin(1);
t90 = -pkin(2) * t180 + t216 * t235 - t320;
t170 = t181 * pkin(7);
t93 = qJDD(2) * pkin(2) + t218 * t235 - t170;
t7 = -t215 * t75 + t217 * (t216 * t90 + t218 * t93);
t226 = qJDD(4) - t7;
t324 = pkin(3) + qJ(5);
t2 = pkin(4) * t92 + qJD(5) * t171 - t133 * t324 + t226;
t6 = -pkin(3) * t133 + t226;
t392 = t6 * mrSges(5,1) + t2 * mrSges(6,1) - t7 * mrSges(4,3) + Ifges(5,4) * t398 + t364 * t335 + (Ifges(5,6) + t366) * t343 + (Ifges(5,2) + t367) * t340;
t314 = t215 * t218;
t315 = t215 * t216;
t8 = t217 * t75 + t93 * t314 + t90 * t315;
t5 = -qJ(4) * t133 + qJD(4) * t171 - t8;
t3 = -pkin(4) * t91 + qJDD(5) - t5;
t391 = -t5 * mrSges(5,1) + t3 * mrSges(6,1) + t8 * mrSges(4,3) + Ifges(4,4) * t340 + Ifges(4,6) * t335 + t362 * t341 + t365 * t398 + (Ifges(4,2) + t363) * t343;
t342 = t91 / 0.2e1;
t390 = -t7 * mrSges(4,1) + t8 * mrSges(4,2) - t6 * mrSges(5,2) - t3 * mrSges(6,2) + t5 * mrSges(5,3) + t2 * mrSges(6,3) - t340 * t364 - t342 * t365 + (-t335 + t398) * t351 + (-Ifges(4,6) + t399) * t343 + (-Ifges(5,4) + t400) * t341;
t389 = -m(6) - m(5);
t388 = t180 / 0.2e1;
t387 = t181 / 0.2e1;
t386 = -mrSges(3,3) + mrSges(2,2);
t295 = qJD(1) * t219;
t275 = t215 * t295;
t293 = qJD(2) * t216;
t172 = t218 * t294 + t293;
t357 = t172 * t217;
t109 = t275 - t357;
t153 = t215 * t306 + t217 * t219;
t110 = qJD(1) * t153 + t215 * t293;
t383 = t363 * t109 + t362 * t110 - t365 * t171;
t382 = -t366 * t109 + t367 * t110 - t364 * t171;
t381 = t221 * Ifges(3,2);
t253 = t221 * mrSges(3,1) - t219 * mrSges(3,2);
t380 = -t253 - mrSges(2,1);
t222 = cos(qJ(1));
t213 = t222 * pkin(7);
t220 = sin(qJ(1));
t300 = t222 * t218;
t312 = t216 * t219;
t201 = qJ(3) * t312;
t297 = t221 * pkin(2) + t201;
t361 = -pkin(1) - t297;
t379 = qJ(3) * t300 + t361 * t220 + t213;
t169 = -pkin(7) * t264 + t208;
t378 = t169 * t221 + t170 * t219;
t352 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t374 = t109 * t399 + t110 * t400 - t351 * t171;
t81 = -mrSges(4,1) * t171 - mrSges(4,3) * t110;
t85 = mrSges(5,1) * t110 - mrSges(5,2) * t171;
t360 = -t81 + t85;
t83 = mrSges(5,1) * t109 + mrSges(5,3) * t171;
t84 = -mrSges(6,1) * t109 - mrSges(6,2) * t171;
t323 = -t83 + t84;
t307 = t218 * t219;
t154 = t215 * t221 + t217 * t307;
t138 = t154 * qJD(1);
t139 = t217 * t294 - t218 * t275;
t310 = t216 * t221;
t325 = pkin(2) * t219;
t240 = -qJ(3) * t310 + t325;
t164 = t240 * qJD(1);
t261 = qJ(3) * t218 + pkin(7);
t248 = qJD(1) * t261;
t166 = t221 * t248;
t96 = t218 * t164 + t166 * t216;
t238 = -qJ(4) * t139 + t96;
t288 = qJD(4) * t215;
t359 = -t138 * t324 + (-qJD(5) * t217 - t288) * t216 - t238;
t165 = t219 * t248;
t143 = t215 * t165;
t256 = t166 * t308 - t143;
t265 = t324 * t219;
t271 = t215 * t290;
t318 = t164 * t217;
t358 = -qJD(5) * t218 + t271 - pkin(4) * t139 - (-qJD(1) * t265 - t318) * t216 - t256;
t174 = t261 * t219;
t356 = t217 * (-t174 * t218 + t216 * t361);
t252 = mrSges(3,1) * t219 + mrSges(3,2) * t221;
t355 = -m(6) * (pkin(4) * t310 - t325) - mrSges(6,1) * t310 + t252;
t353 = t216 ^ 2 + t218 ^ 2;
t122 = qJD(2) * t240 - t216 * t289;
t247 = qJD(2) * t261;
t123 = -t219 * t247 + t269;
t124 = -t218 * t289 - t221 * t247;
t40 = t122 * t315 + t217 * t123 + t124 * t314;
t26 = -t216 * (qJ(4) * t292 - qJD(4) * t221) - t40;
t282 = pkin(7) * t294;
t128 = qJ(3) * t172 + t282;
t147 = qJD(2) * pkin(2) - t165;
t148 = t361 * qJD(1);
t43 = -t215 * t128 + t217 * (t147 * t218 + t148 * t216);
t339 = -t109 / 0.2e1;
t338 = t109 / 0.2e1;
t337 = -t110 / 0.2e1;
t336 = t110 / 0.2e1;
t322 = Ifges(3,4) * t219;
t321 = Ifges(3,4) * t221;
t319 = t122 * t217;
t313 = t216 * t217;
t311 = t216 * t220;
t309 = t216 * t222;
t305 = t219 * t220;
t304 = t219 * t222;
t303 = t220 * t218;
t302 = t220 * t221;
t301 = t221 * t222;
t175 = t261 * t221;
t158 = t215 * t175;
t299 = pkin(3) * t310 + t158;
t157 = pkin(2) * t314 + qJ(3) * t313;
t296 = t222 * pkin(1) + t220 * pkin(7);
t291 = qJD(2) * t221;
t284 = pkin(2) * t305;
t283 = pkin(2) * t304;
t280 = t216 * t302;
t279 = t216 * t301;
t278 = t215 * t307;
t277 = t217 * t306;
t44 = t217 * t128 + t147 * t314 + t148 * t315;
t65 = t164 * t315 - t217 * t165 - t166 * t314;
t71 = -t174 * t314 + t217 * t175 + t315 * t361;
t276 = -pkin(2) * t217 - pkin(3);
t274 = t216 * t295;
t272 = t216 * t292;
t270 = t217 * t290;
t34 = -t92 * mrSges(6,2) + t91 * mrSges(6,3);
t62 = t92 * mrSges(5,1) + t133 * mrSges(5,2);
t61 = -t91 * mrSges(6,1) + t133 * mrSges(6,2);
t59 = t92 * mrSges(6,1) - t133 * mrSges(6,3);
t262 = -qJ(4) * t215 - pkin(2);
t260 = t287 / 0.2e1;
t74 = t218 * t122 - t124 * t216;
t97 = t174 * t216 + t218 * t361;
t111 = t215 * t123;
t257 = -t124 * t308 + t111;
t129 = t215 * t305 - t220 * t277;
t130 = t153 * t220;
t187 = qJ(3) * t280;
t255 = -t130 * pkin(3) - qJ(4) * t129 + t187;
t131 = t215 * t304 - t222 * t277;
t132 = t153 * t222;
t188 = qJ(3) * t279;
t254 = -t132 * pkin(3) - t131 * qJ(4) + t188;
t134 = -qJ(4) * t218 - t157;
t30 = -t216 * t93 + t218 * t90 + qJDD(3);
t79 = -t147 * t216 + t218 * t148 + qJD(3);
t251 = t322 + t381;
t250 = Ifges(3,5) * t221 - Ifges(3,6) * t219;
t246 = pkin(2) * t301 + qJ(3) * t303 + t222 * t201 + t296;
t33 = qJ(4) * t171 - t44;
t239 = pkin(1) * t252;
t237 = -qJ(4) * t153 + t97;
t236 = t219 * (Ifges(3,1) * t221 - t322);
t63 = qJ(4) * t310 - t71;
t101 = t215 * t302 + (t219 * t303 + t309) * t217;
t103 = t154 * t222 - t217 * t311;
t152 = t215 * t219 - t277;
t232 = -g(1) * t103 - g(2) * t101 - g(3) * t152;
t231 = -qJ(4) * t110 + t79;
t54 = -qJ(4) * t274 - t65;
t160 = t216 * t305 - t300;
t161 = t216 * t304 + t303;
t230 = -g(1) * t161 - g(2) * t160 + g(3) * t310;
t141 = -qJD(2) * t278 + t217 * t291;
t227 = -qJ(4) * t141 - qJD(4) * t153 + t74;
t225 = qJD(4) - t43;
t224 = -qJ(4) * t92 - qJD(4) * t110 + t30;
t209 = Ifges(3,4) * t294;
t199 = qJ(3) * t315;
t191 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t294;
t190 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t295;
t168 = qJD(4) * t218 + t270;
t163 = Ifges(3,1) * t295 + Ifges(3,5) * qJD(2) + t209;
t162 = Ifges(3,6) * qJD(2) + qJD(1) * t251;
t156 = pkin(2) * t308 - t199;
t155 = t217 * t221 - t278;
t140 = t154 * qJD(2);
t136 = (-pkin(3) * t217 + t262) * t216;
t135 = t218 * t276 + t199;
t108 = (-t217 * t324 + t262) * t216;
t107 = pkin(4) * t313 - t134;
t104 = t215 * t311 + t217 * t301 - t222 * t278;
t102 = -t215 * t309 + t217 * t302 - t220 * t278;
t98 = pkin(4) * t315 + t199 + (-qJ(5) + t276) * t218;
t95 = -t275 * t353 + t357;
t94 = t217 * t295 * t353 + t172 * t215;
t89 = t92 * mrSges(5,3);
t87 = t92 * mrSges(4,2);
t82 = mrSges(6,1) * t110 + mrSges(6,3) * t171;
t80 = mrSges(4,2) * t171 - mrSges(4,3) * t109;
t70 = -t158 + t356;
t69 = -mrSges(5,2) * t109 - mrSges(5,3) * t110;
t68 = mrSges(4,1) * t109 + mrSges(4,2) * t110;
t67 = -mrSges(6,2) * t110 + mrSges(6,3) * t109;
t66 = t299 - t356;
t64 = t143 + (t164 * t216 - t166 * t218) * t217;
t60 = mrSges(5,1) * t91 - mrSges(5,3) * t133;
t58 = mrSges(4,1) * t133 - mrSges(4,3) * t92;
t57 = -mrSges(4,2) * t133 - mrSges(4,3) * t91;
t56 = pkin(3) * t152 + t237;
t55 = (-pkin(3) * t295 - t318) * t216 + t256;
t50 = t110 * Ifges(4,4) - t109 * Ifges(4,2) - t171 * Ifges(4,6);
t49 = -t171 * Ifges(5,4) - t110 * Ifges(5,2) + t109 * Ifges(5,6);
t42 = pkin(3) * t138 + t238;
t41 = -pkin(4) * t152 - t63;
t39 = -t111 + (t122 * t216 + t124 * t218) * t217;
t38 = t152 * t324 + t237;
t37 = t174 * t308 + pkin(4) * t153 + (qJ(5) * t221 - t217 * t361) * t216 + t299;
t36 = -t91 * mrSges(5,2) - t89;
t35 = t91 * mrSges(4,1) + t87;
t32 = pkin(3) * t171 + t225;
t31 = -pkin(4) * t138 - t54;
t29 = (-pkin(3) * t292 - t319) * t216 + t257;
t25 = pkin(3) * t109 + t231;
t15 = pkin(3) * t140 + t227;
t14 = -pkin(4) * t109 + qJD(5) - t33;
t13 = -pkin(4) * t140 - t26;
t12 = pkin(4) * t141 + (-qJD(2) * t265 + qJD(5) * t221 - t319) * t216 + t257;
t11 = t109 * t324 + t231;
t10 = pkin(4) * t110 + t171 * t324 + t225;
t9 = qJD(5) * t152 + t140 * t324 + t227;
t4 = pkin(3) * t91 + t224;
t1 = qJD(5) * t109 + t324 * t91 + t224;
t16 = [(t30 * mrSges(4,1) - t4 * mrSges(5,2) + t1 * mrSges(6,3) - Ifges(4,2) * t343 + Ifges(5,6) * t341 + t335 * t399 - t366 * t340 + t363 * t342 - t391) * t152 - (t140 * t399 + t141 * t400 + t272 * t351) * t171 / 0.2e1 + (t30 * mrSges(4,2) - t1 * mrSges(6,2) - t4 * mrSges(5,3) + Ifges(4,4) * t343 - Ifges(5,2) * t341 + t335 * t400 + t367 * t340 + t362 * t342 + t392) * t153 + (Ifges(3,6) * qJDD(2) + t321 * t260 + pkin(7) * (-qJDD(2) * mrSges(3,2) + mrSges(3,3) * t180) + Ifges(3,4) * t387 + Ifges(3,2) * t388) * t221 + (-t260 * t381 + Ifges(3,1) * t181 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t181) + Ifges(3,4) * t388 + Ifges(3,5) * qJDD(2)) * t219 + (t140 * t363 + t141 * t362 + t272 * t365) * t338 + (-t140 * t366 + t141 * t367 + t272 * t364) * t336 + t382 * t141 / 0.2e1 + t383 * t140 / 0.2e1 + (-m(3) * t213 - t379 * m(4) + t389 * (-t102 * pkin(3) - qJ(4) * t101 + t379) + t386 * t222 + (m(3) * pkin(1) - t380) * t220 + t371 * t160 + t393 * t102 - t352 * t101) * g(1) + (-m(3) * t296 - m(4) * t246 + t389 * (t104 * pkin(3) + qJ(4) * t103 + t246) + t380 * t222 + t386 * t220 - t371 * t161 - t393 * t104 + t352 * t103) * g(2) + t378 * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t378) + t253 * t320 - t162 * t292 / 0.2e1 + m(4) * (t30 * t97 + t39 * t43 + t40 * t44 + t7 * t70 + t71 * t8 + t74 * t79) + m(5) * (t15 * t25 + t26 * t33 + t29 * t32 + t4 * t56 + t5 * t63 + t6 * t66) + m(6) * (t1 * t38 + t10 * t12 + t11 * t9 + t13 * t14 + t2 * t37 + t3 * t41) - t239 * t287 - t191 * t281 + (t163 / 0.2e1 - pkin(7) * t190) * t291 + t10 * (mrSges(6,1) * t141 - mrSges(6,3) * t272) + t33 * (mrSges(5,1) * t140 - mrSges(5,3) * t272) + t44 * (-mrSges(4,2) * t272 - mrSges(4,3) * t140) + t14 * (-mrSges(6,1) * t140 + mrSges(6,2) * t272) + t32 * (mrSges(5,1) * t141 + mrSges(5,2) * t272) + t43 * (mrSges(4,1) * t272 - mrSges(4,3) * t141) + t374 * t272 / 0.2e1 + qJD(2) ^ 2 * t250 / 0.2e1 + t321 * t387 + t251 * t388 + t236 * t260 - t140 * t50 / 0.2e1 + t11 * (-mrSges(6,2) * t141 + mrSges(6,3) * t140) + t79 * (mrSges(4,1) * t140 + mrSges(4,2) * t141) + t25 * (-mrSges(5,2) * t140 - mrSges(5,3) * t141) - t141 * t49 / 0.2e1 + t390 * t310 - pkin(1) * (-mrSges(3,1) * t180 + mrSges(3,2) * t181) + (Ifges(5,4) * t272 - Ifges(5,2) * t141 + Ifges(5,6) * t140) * t337 + (Ifges(4,4) * t141 - Ifges(4,2) * t140 + Ifges(4,6) * t272) * t339 + t38 * t34 + t56 * t36 + t37 * t59 + t41 * t61 + t63 * t60 + t66 * t62 + t9 * t67 + t15 * t69 + t70 * t58 + t71 * t57 + t74 * t68 + t40 * t80 + t39 * t81 + t12 * t82 + t26 * t83 + t13 * t84 + t29 * t85 + t97 * t35 + Ifges(2,3) * qJDD(1); (t134 * t5 + t135 * t6 + t136 * t4 - t25 * t42 - t32 * t55 + (-t168 - t54) * t33) * m(5) + (t162 / 0.2e1 + pkin(7) * t191) * t295 + (t270 - t65) * t80 + (t1 * (-mrSges(6,2) * t215 - mrSges(6,3) * t217) + t30 * (-mrSges(4,1) * t217 + mrSges(4,2) * t215) + t4 * (mrSges(5,2) * t217 - mrSges(5,3) * t215) - t69 * t288 - pkin(2) * t35 + (-Ifges(5,2) * t215 - Ifges(5,6) * t217) * t341 + (Ifges(4,4) * t215 + Ifges(4,2) * t217) * t343 + m(4) * (-pkin(2) * t30 + (-t215 * t43 + t217 * t44) * qJD(3)) + (t215 * t362 - t217 * t363) * t342 + (t215 * t367 + t217 * t366) * t340 + (t215 * t400 - t217 * t399) * t335) * t216 + (t138 * t399 + t139 * t400 + t274 * t351) * t171 / 0.2e1 - t382 * t139 / 0.2e1 - t383 * t138 / 0.2e1 + t358 * t82 + (t1 * t108 + t107 * t3 + t2 * t98 + (t168 - t31) * t14 + t359 * t11 + t358 * t10) * m(6) + (-m(4) * (t187 - t284) - m(5) * (t255 - t284) - m(6) * t255 + t368 * t280 + t355 * t220 + t393 * t130 - t352 * t129) * g(2) + (-m(4) * (t188 - t283) - m(5) * (t254 - t283) - m(6) * t254 + t368 * t279 + t355 * t222 + t393 * t132 - t352 * t131) * g(1) + (-m(4) * t297 - t253 + t389 * (t155 * pkin(3) + qJ(4) * t154 + t297) - t371 * t312 - t393 * t155 + t352 * t154) * g(3) - t250 * t287 / 0.2e1 + (-t138 * t366 + t139 * t367 + t274 * t364) * t337 - t10 * (mrSges(6,1) * t139 - mrSges(6,3) * t274) - t33 * (mrSges(5,1) * t138 - mrSges(5,3) * t274) - t44 * (-mrSges(4,2) * t274 - mrSges(4,3) * t138) - t14 * (-mrSges(6,1) * t138 + mrSges(6,2) * t274) - t32 * (mrSges(5,1) * t139 + mrSges(5,2) * t274) - t43 * (mrSges(4,1) * t274 - mrSges(4,3) * t139) - t374 * t274 / 0.2e1 + t190 * t282 + (m(5) * (qJD(3) * t32 - qJD(4) * t25) + t392) * t315 + t323 * t168 + t360 * t271 + (t138 * t363 + t139 * t362 + t274 * t365) * t339 - (-Ifges(3,2) * t295 + t163 + t209) * t294 / 0.2e1 + (t239 - t236 / 0.2e1) * qJD(1) ^ 2 + (t156 * t7 + t157 * t8 - t43 * t64 - t44 * t65 - t79 * t96) * m(4) + t138 * t50 / 0.2e1 - t79 * (mrSges(4,1) * t138 + mrSges(4,2) * t139) - t11 * (-mrSges(6,2) * t139 + mrSges(6,3) * t138) + t139 * t49 / 0.2e1 - t25 * (-mrSges(5,2) * t138 - mrSges(5,3) * t139) + t156 * t58 + t157 * t57 - t169 * mrSges(3,2) - t170 * mrSges(3,1) - t390 * t218 + t391 * t313 + Ifges(3,6) * t180 + Ifges(3,5) * t181 + (Ifges(5,4) * t274 - Ifges(5,2) * t139 + Ifges(5,6) * t138) * t336 + (Ifges(4,4) * t139 - Ifges(4,2) * t138 + Ifges(4,6) * t274) * t338 + Ifges(3,3) * qJDD(2) + t359 * t67 - t42 * t69 - t64 * t81 - t54 * t83 - t31 * t84 - t55 * t85 - t96 * t68 + t98 * t59 + t107 * t61 + t108 * t34 + t134 * t60 + t135 * t62 + t136 * t36; t87 - t89 + t369 * t91 + (-t80 - t323) * t95 + (-t82 - t360) * t94 + t34 + (-t10 * t94 - t14 * t95 + t1 + t230) * m(6) + (-t32 * t94 + t33 * t95 + t230 + t4) * m(5) + (t43 * t94 - t44 * t95 + t230 + t30) * m(4); t323 * t171 + (t67 + t69) * t110 + t59 + t62 + (t11 * t110 + t14 * t171 + t2 + t232) * m(6) + (t110 * t25 - t171 * t33 + t232 + t6) * m(5); -t109 * t67 - t171 * t82 + (-g(1) * t104 - g(2) * t102 - g(3) * t153 - t10 * t171 - t11 * t109 + t3) * m(6) + t61;];
tau = t16;
