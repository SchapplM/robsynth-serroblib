% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:20:46
% EndTime: 2018-11-23 17:20:55
% DurationCPUTime: 8.67s
% Computational Cost: add. (23912->575), mult. (63430->773), div. (0->0), fcn. (49911->10), ass. (0->278)
t251 = sin(pkin(11));
t252 = cos(pkin(11));
t256 = sin(qJ(2));
t260 = cos(qJ(2));
t226 = -t251 * t256 + t252 * t260;
t215 = t226 * qJD(1);
t316 = qJD(1) * t260;
t317 = qJD(1) * t256;
t216 = -t251 * t316 - t252 * t317;
t255 = sin(qJ(4));
t259 = cos(qJ(4));
t175 = t215 * t255 - t216 * t259;
t254 = sin(qJ(5));
t258 = cos(qJ(5));
t297 = t259 * t215 + t216 * t255;
t394 = -t175 * t254 + t258 * t297;
t116 = qJD(6) - t394;
t393 = t258 * t175 + t254 * t297;
t81 = pkin(5) * t393 - pkin(10) * t394;
t169 = Ifges(5,4) * t297;
t250 = qJD(2) + qJD(4);
t111 = t175 * Ifges(5,1) + t250 * Ifges(5,5) + t169;
t227 = t251 * t260 + t252 * t256;
t217 = t227 * qJD(2);
t205 = qJD(1) * t217;
t218 = t226 * qJD(2);
t206 = qJD(1) * t218;
t117 = qJD(4) * t297 - t205 * t255 + t206 * t259;
t118 = -qJD(4) * t175 - t205 * t259 - t206 * t255;
t245 = -pkin(2) * t260 - pkin(1);
t318 = qJD(1) * t245;
t233 = qJD(3) + t318;
t184 = -t215 * pkin(3) + t233;
t249 = qJD(5) + t250;
t253 = sin(qJ(6));
t257 = cos(qJ(6));
t101 = t249 * t257 - t253 * t393;
t61 = qJD(5) * t394 + t117 * t258 + t118 * t254;
t35 = qJD(6) * t101 + t257 * t61;
t102 = t249 * t253 + t257 * t393;
t36 = -qJD(6) * t102 - t253 * t61;
t62 = qJD(5) * t393 + t117 * t254 - t258 * t118;
t12 = Ifges(7,4) * t35 + Ifges(7,2) * t36 + Ifges(7,6) * t62;
t13 = Ifges(7,1) * t35 + Ifges(7,4) * t36 + Ifges(7,5) * t62;
t287 = Ifges(7,5) * t253 + Ifges(7,6) * t257;
t341 = Ifges(7,4) * t253;
t289 = Ifges(7,2) * t257 + t341;
t340 = Ifges(7,4) * t257;
t291 = Ifges(7,1) * t253 + t340;
t294 = mrSges(7,1) * t257 - mrSges(7,2) * t253;
t247 = pkin(2) * t317;
t241 = qJD(2) * t247;
t185 = pkin(3) * t205 + t241;
t93 = -pkin(4) * t118 + t185;
t17 = pkin(5) * t62 - pkin(10) * t61 + t93;
t388 = pkin(9) * t297;
t347 = -qJ(3) - pkin(7);
t237 = t347 * t260;
t232 = qJD(1) * t237;
t219 = t251 * t232;
t236 = t347 * t256;
t231 = qJD(1) * t236;
t225 = qJD(2) * pkin(2) + t231;
t178 = t252 * t225 + t219;
t352 = pkin(8) * t216;
t142 = qJD(2) * pkin(3) + t178 + t352;
t319 = t252 * t232;
t179 = t251 * t225 - t319;
t353 = pkin(8) * t215;
t147 = t179 + t353;
t92 = t142 * t255 + t147 * t259;
t85 = t92 + t388;
t328 = t258 * t85;
t397 = pkin(9) * t175;
t91 = t259 * t142 - t147 * t255;
t84 = t91 - t397;
t82 = pkin(4) * t250 + t84;
t43 = t254 * t82 + t328;
t41 = pkin(10) * t249 + t43;
t133 = -pkin(4) * t297 + t184;
t67 = -pkin(5) * t394 - pkin(10) * t393 + t133;
t20 = -t253 * t41 + t257 * t67;
t301 = qJD(2) * t347;
t213 = qJD(3) * t260 + t256 * t301;
t194 = t213 * qJD(1);
t214 = -t256 * qJD(3) + t260 * t301;
t195 = t214 * qJD(1);
t148 = -t194 * t251 + t195 * t252;
t275 = -pkin(8) * t206 + t148;
t149 = t194 * t252 + t195 * t251;
t276 = -pkin(8) * t205 + t149;
t64 = -qJD(4) * t92 - t255 * t276 + t259 * t275;
t264 = -t117 * pkin(9) + t64;
t313 = qJD(4) * t259;
t314 = qJD(4) * t255;
t63 = t142 * t313 - t147 * t314 + t255 * t275 + t259 * t276;
t38 = pkin(9) * t118 + t63;
t329 = t254 * t85;
t42 = t258 * t82 - t329;
t8 = qJD(5) * t42 + t254 * t264 + t258 * t38;
t2 = qJD(6) * t20 + t17 * t253 + t257 * t8;
t351 = t2 * t257;
t355 = t257 / 0.2e1;
t373 = t62 / 0.2e1;
t374 = t36 / 0.2e1;
t375 = t35 / 0.2e1;
t288 = Ifges(7,5) * t257 - Ifges(7,6) * t253;
t290 = -Ifges(7,2) * t253 + t340;
t292 = Ifges(7,1) * t257 - t341;
t293 = mrSges(7,1) * t253 + mrSges(7,2) * t257;
t356 = -t253 / 0.2e1;
t367 = t102 / 0.2e1;
t40 = -pkin(5) * t249 - t42;
t333 = t102 * Ifges(7,4);
t53 = t101 * Ifges(7,2) + t116 * Ifges(7,6) + t333;
t100 = Ifges(7,4) * t101;
t54 = t102 * Ifges(7,1) + t116 * Ifges(7,5) + t100;
t392 = t40 * t293 + t53 * t356 + t54 * t355 + t101 * t290 / 0.2e1 + t292 * t367 + t116 * t288 / 0.2e1;
t9 = qJD(5) * t43 + t254 * t38 - t258 * t264;
t269 = -t8 * mrSges(6,2) + mrSges(7,3) * t351 + t253 * t13 / 0.2e1 + t12 * t355 + t291 * t375 + t289 * t374 + t287 * t373 - Ifges(6,6) * t62 + Ifges(6,5) * t61 + (-mrSges(6,1) - t294) * t9 + t392 * qJD(6);
t343 = Ifges(5,4) * t175;
t363 = -t175 / 0.2e1;
t400 = t64 * mrSges(5,1) - t63 * mrSges(5,2) + Ifges(5,5) * t117 + Ifges(5,6) * t118 + t269 - (Ifges(5,5) * t297 - Ifges(5,6) * t175) * t250 / 0.2e1 - (-Ifges(5,2) * t175 + t111 + t169) * t297 / 0.2e1 - t184 * (mrSges(5,1) * t175 + mrSges(5,2) * t297) + (Ifges(5,1) * t297 - t343) * t363;
t399 = Ifges(6,2) / 0.2e1;
t398 = pkin(4) * t175;
t16 = -mrSges(7,1) * t36 + mrSges(7,2) * t35;
t396 = m(7) * t9 + t16;
t395 = t175 * t92 + t297 * t91;
t110 = Ifges(5,2) * t297 + t250 * Ifges(5,6) + t343;
t390 = t110 / 0.2e1;
t389 = t394 * t399;
t21 = t253 * t67 + t257 * t41;
t331 = t20 * t257;
t285 = t21 * t253 + t331;
t277 = t285 * mrSges(7,3);
t242 = pkin(2) * t252 + pkin(3);
t354 = pkin(2) * t251;
t211 = t259 * t242 - t255 * t354;
t209 = pkin(4) + t211;
t212 = t242 * t255 + t259 * t354;
t161 = t254 * t209 + t258 * t212;
t202 = t211 * qJD(4);
t203 = t212 * qJD(4);
t182 = -t231 * t251 + t319;
t150 = t182 - t353;
t183 = t252 * t231 + t219;
t151 = t183 + t352;
t96 = t259 * t150 - t151 * t255;
t280 = t96 - t388;
t97 = t255 * t150 + t259 * t151;
t87 = t97 - t397;
t383 = -qJD(5) * t161 + (-t203 - t280) * t258 + (-t202 + t87) * t254;
t326 = mrSges(6,1) * t249 + mrSges(7,1) * t101 - mrSges(7,2) * t102 - mrSges(6,3) * t393;
t382 = -t96 - t203;
t381 = -t97 + t202;
t186 = t252 * t236 + t237 * t251;
t167 = -pkin(8) * t227 + t186;
t187 = t251 * t236 - t252 * t237;
t168 = pkin(8) * t226 + t187;
t108 = t255 * t167 + t259 * t168;
t115 = Ifges(6,4) * t394;
t307 = -t115 / 0.2e1;
t346 = Ifges(6,1) * t393;
t379 = t307 - t346 / 0.2e1;
t378 = -t20 * t253 + t21 * t257;
t3 = -qJD(6) * t21 + t17 * t257 - t253 * t8;
t377 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t35 + Ifges(7,6) * t36;
t372 = pkin(1) * mrSges(3,1);
t371 = pkin(1) * mrSges(3,2);
t107 = t259 * t167 - t168 * t255;
t181 = t226 * t255 + t227 * t259;
t89 = -pkin(9) * t181 + t107;
t180 = t226 * t259 - t227 * t255;
t90 = pkin(9) * t180 + t108;
t65 = t254 * t90 - t258 * t89;
t370 = t65 * t9;
t369 = -t101 / 0.2e1;
t368 = -t102 / 0.2e1;
t366 = -t116 / 0.2e1;
t364 = t297 / 0.2e1;
t362 = t175 / 0.2e1;
t360 = -t216 / 0.2e1;
t359 = -t217 / 0.2e1;
t358 = t218 / 0.2e1;
t350 = t3 * t253;
t349 = t42 * mrSges(6,3);
t348 = t43 * mrSges(6,3);
t345 = Ifges(3,4) * t256;
t344 = Ifges(4,4) * t216;
t325 = Ifges(3,5) * qJD(2);
t324 = Ifges(3,6) * qJD(2);
t323 = qJD(2) * mrSges(3,1);
t322 = qJD(2) * mrSges(3,2);
t166 = t252 * t213 + t251 * t214;
t315 = qJD(2) * t256;
t312 = qJD(6) * t253;
t311 = qJD(6) * t257;
t304 = t325 / 0.2e1;
t303 = -t324 / 0.2e1;
t302 = t62 * mrSges(6,1) + t61 * mrSges(6,2);
t188 = -pkin(3) * t216 + t247;
t189 = pkin(2) * t315 + pkin(3) * t217;
t300 = t205 * mrSges(4,1) + t206 * mrSges(4,2);
t299 = -t118 * mrSges(5,1) + t117 * mrSges(5,2);
t165 = -t213 * t251 + t252 * t214;
t295 = -t2 * t253 - t257 * t3;
t18 = mrSges(7,1) * t62 - mrSges(7,3) * t35;
t19 = -mrSges(7,2) * t62 + mrSges(7,3) * t36;
t286 = -t253 * t18 + t257 * t19;
t66 = t254 * t89 + t258 * t90;
t130 = t180 * t254 + t181 * t258;
t196 = -t226 * pkin(3) + t245;
t143 = -t180 * pkin(4) + t196;
t281 = t258 * t180 - t181 * t254;
t74 = -pkin(5) * t281 - t130 * pkin(10) + t143;
t30 = t253 * t74 + t257 * t66;
t29 = -t253 * t66 + t257 * t74;
t76 = -mrSges(7,2) * t116 + mrSges(7,3) * t101;
t77 = mrSges(7,1) * t116 - mrSges(7,3) * t102;
t283 = -t253 * t77 + t257 * t76;
t160 = t209 * t258 - t212 * t254;
t128 = -qJD(4) * t181 - t217 * t259 - t218 * t255;
t99 = -pkin(4) * t128 + t189;
t134 = t188 + t398;
t103 = -mrSges(6,2) * t249 + mrSges(6,3) * t394;
t279 = -t103 - t283;
t139 = -pkin(8) * t218 + t165;
t140 = -pkin(8) * t217 + t166;
t71 = t255 * t139 + t259 * t140 + t167 * t313 - t168 * t314;
t271 = -qJD(6) * t285 - t350;
t72 = -qJD(4) * t108 + t259 * t139 - t140 * t255;
t127 = qJD(4) * t180 - t217 * t255 + t218 * t259;
t270 = -pkin(9) * t127 + t72;
t267 = -t133 * mrSges(6,1) - t20 * mrSges(7,1) + t21 * mrSges(7,2) + Ifges(6,4) * t393 - Ifges(7,5) * t102 + Ifges(6,6) * t249 - Ifges(7,6) * t101 - Ifges(7,3) * t116 + t389;
t266 = -t77 * t311 - t76 * t312 + m(7) * (-t20 * t311 - t21 * t312 - t350 + t351) + t286;
t265 = t389 + t267;
t263 = t133 * mrSges(6,2) + t115 / 0.2e1 + Ifges(6,5) * t249 + t346 / 0.2e1 + t392;
t262 = -t263 + t379;
t261 = -t263 + t277;
t246 = Ifges(3,4) * t316;
t235 = mrSges(3,3) * t316 - t322;
t234 = -mrSges(3,3) * t317 + t323;
t224 = Ifges(3,1) * t317 + t246 + t325;
t223 = t324 + (Ifges(3,2) * t260 + t345) * qJD(1);
t210 = Ifges(4,4) * t215;
t192 = qJD(2) * mrSges(4,1) + mrSges(4,3) * t216;
t191 = -qJD(2) * mrSges(4,2) + t215 * mrSges(4,3);
t177 = -mrSges(4,1) * t215 - mrSges(4,2) * t216;
t171 = -t216 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t210;
t170 = t215 * Ifges(4,2) + Ifges(4,6) * qJD(2) - t344;
t159 = pkin(10) + t161;
t158 = -pkin(5) - t160;
t153 = mrSges(5,1) * t250 - mrSges(5,3) * t175;
t152 = -mrSges(5,2) * t250 + mrSges(5,3) * t297;
t125 = -mrSges(5,1) * t297 + mrSges(5,2) * t175;
t105 = qJD(5) * t160 + t202 * t258 - t203 * t254;
t80 = -mrSges(6,1) * t394 + mrSges(6,2) * t393;
t73 = t81 + t398;
t70 = qJD(5) * t130 + t127 * t254 - t258 * t128;
t69 = qJD(5) * t281 + t127 * t258 + t128 * t254;
t68 = t134 + t81;
t57 = Ifges(7,3) * t62;
t51 = pkin(9) * t128 + t71;
t47 = t254 * t280 + t258 * t87;
t45 = t258 * t84 - t329;
t44 = t254 * t84 + t328;
t28 = t253 * t81 + t257 * t42;
t27 = -t253 * t42 + t257 * t81;
t26 = t253 * t73 + t257 * t45;
t25 = -t253 * t45 + t257 * t73;
t24 = t253 * t68 + t257 * t47;
t23 = -t253 * t47 + t257 * t68;
t22 = pkin(5) * t70 - pkin(10) * t69 + t99;
t15 = qJD(5) * t66 + t254 * t51 - t258 * t270;
t14 = -qJD(5) * t65 + t254 * t270 + t258 * t51;
t5 = -qJD(6) * t30 - t14 * t253 + t22 * t257;
t4 = qJD(6) * t29 + t14 * t257 + t22 * t253;
t1 = [-t265 * t70 - t326 * t15 + (-t226 * t205 + t215 * t359) * Ifges(4,2) + (-t227 * t205 + t226 * t206 + t215 * t358 - t217 * t360) * Ifges(4,4) + (-t148 * t227 + t149 * t226 - t178 * t218 - t179 * t217 - t186 * t206 - t187 * t205) * mrSges(4,3) + t128 * t390 + m(5) * (t107 * t64 + t108 * t63 + t184 * t189 + t185 * t196 + t71 * t92 + t72 * t91) + (-t42 * t69 - t43 * t70 + t61 * t65 - t62 * t66) * mrSges(6,3) + (t227 * t206 + t218 * t360) * Ifges(4,1) + (t181 * t117 + t127 * t362) * Ifges(5,1) + (t180 * t118 + t128 * t364) * Ifges(5,2) + (t180 * t117 + t181 * t118 + t127 * t364 + t128 * t362) * Ifges(5,4) + t143 * t302 + t196 * t299 + t245 * t300 + m(7) * (t15 * t40 + t2 * t30 + t20 * t5 + t21 * t4 + t29 * t3 + t370) + m(6) * (t133 * t99 + t14 * t43 + t143 * t93 - t15 * t42 + t66 * t8 + t370) + (Ifges(4,5) * t358 + Ifges(4,6) * t359 + (-pkin(7) * t234 + t224 / 0.2e1 + t304 + (-0.2e1 * t371 + 0.3e1 / 0.2e1 * Ifges(3,4) * t260) * qJD(1)) * t260) * qJD(2) + (-pkin(7) * t235 - t223 / 0.2e1 + t303 + (-0.2e1 * t372 - 0.3e1 / 0.2e1 * t345 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t260) * qJD(1) + (t177 + qJD(1) * (-mrSges(4,1) * t226 + mrSges(4,2) * t227) + m(4) * (t233 + t318)) * pkin(2)) * t315 + (t93 * mrSges(6,2) - Ifges(6,4) * t62 + Ifges(6,1) * t61 + t288 * t373 + t292 * t375 + t290 * t374 + t12 * t356 + t13 * t355 + (mrSges(6,3) + t293) * t9 + t295 * mrSges(7,3) + (-t257 * t53 / 0.2e1 + t54 * t356 + t287 * t366 + t289 * t369 + t291 * t368 + t40 * t294 - t378 * mrSges(7,3)) * qJD(6)) * t130 - (t93 * mrSges(6,1) - Ifges(6,4) * t61 + t57 / 0.2e1 - t8 * mrSges(6,3) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t62 + t377) * t281 + (-t107 * t117 + t108 * t118 - t127 * t91 + t128 * t92 + t180 * t63 - t181 * t64) * mrSges(5,3) + m(4) * (t148 * t186 + t149 * t187 + t165 * t178 + t166 * t179) + t127 * t111 / 0.2e1 + t14 * t103 + t99 * t80 + t5 * t77 + t4 * t76 + t65 * t16 + t30 * t19 + t29 * t18 + t171 * t358 + t170 * t359 + t71 * t152 + t72 * t153 + t184 * (-mrSges(5,1) * t128 + mrSges(5,2) * t127) + t185 * (-mrSges(5,1) * t180 + mrSges(5,2) * t181) + t189 * t125 + t166 * t191 + t165 * t192 + t250 * (Ifges(5,5) * t127 + Ifges(5,6) * t128) / 0.2e1 + t233 * (mrSges(4,1) * t217 + mrSges(4,2) * t218) + (-t277 - t262) * t69; t400 + (-t117 * t211 + t118 * t212 + t395) * mrSges(5,3) + t286 * t159 + m(4) * (t148 * t252 + t149 * t251) * pkin(2) + (t178 * t215 - t179 * t216 + (-t205 * t251 - t206 * t252) * pkin(2)) * mrSges(4,3) - mrSges(7,3) * t350 - (Ifges(4,2) * t216 + t171 + t210) * t215 / 0.2e1 + ((t304 - t224 / 0.2e1 - t246 / 0.2e1 + qJD(1) * t371 + (t234 - t323) * pkin(7)) * t260 + (t303 + t223 / 0.2e1 + (t372 + t345 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t260) * qJD(1) + (t235 + t322) * pkin(7) + (-m(4) * t233 - t177) * pkin(2)) * t256) * qJD(1) + t265 * t393 + (-t160 * t61 - t161 * t62 + t393 * t43 + t394 * t42) * mrSges(6,3) + (t261 + t379) * t394 - m(4) * (t178 * t182 + t179 * t183) + t175 * t390 - t134 * t80 + t381 * t152 + (-t184 * t188 + t211 * t64 + t212 * t63 + t381 * t92 + t382 * t91) * m(5) + t382 * t153 + (-t133 * t134 - t160 * t9 + t161 * t8 + (t105 - t47) * t43 + t383 * t42) * m(6) + t326 * t383 + (t158 * t9 - t20 * t23 - t21 * t24 + (t271 + t351) * t159 - t383 * t40 + t378 * t105) * m(7) - t279 * t105 - t47 * t103 - t23 * t77 - t24 * t76 + t216 * (Ifges(4,1) * t215 + t344) / 0.2e1 + t148 * mrSges(4,1) - t149 * mrSges(4,2) + t158 * t16 - t188 * t125 - t183 * t191 - t182 * t192 - Ifges(4,6) * t205 + Ifges(4,5) * t206 - qJD(2) * (Ifges(4,5) * t215 + Ifges(4,6) * t216) / 0.2e1 + t170 * t360 - t233 * (-mrSges(4,1) * t216 + mrSges(4,2) * t215) + ((-t253 * t76 - t257 * t77) * t159 - t277) * qJD(6); t283 * qJD(6) + t326 * t393 + t279 * t394 - t297 * t152 + t175 * t153 + t257 * t18 + t253 * t19 - t215 * t191 - t216 * t192 + t299 + t300 + t302 + (t116 * t378 - t393 * t40 - t295) * m(7) + (t393 * t42 - t394 * t43 + t93) * m(6) + (t175 * t91 - t297 * t92 + t185) * m(5) + (-t178 * t216 - t179 * t215 + t241) * m(4); t266 * (pkin(4) * t254 + pkin(10)) - m(6) * (-t42 * t44 + t43 * t45) + t326 * t44 + (-t116 * t331 + (-t116 * t21 - t3) * t253) * mrSges(7,3) + (t265 + t348) * t393 + (t262 + t349) * t394 + t396 * (-pkin(4) * t258 - pkin(5)) + (-t175 * t80 + (-t254 * t62 - t258 * t61) * mrSges(6,3) + (-t326 * t254 - t279 * t258 + m(7) * (t254 * t40 + t258 * t378)) * qJD(5) + (0.2e1 * t133 * t363 + t254 * t8 - t258 * t9 + (-t254 * t42 + t258 * t43) * qJD(5)) * m(6)) * pkin(4) + t395 * mrSges(5,3) - m(7) * (t20 * t25 + t21 * t26 + t40 * t44) - t45 * t103 - t25 * t77 - t26 * t76 - t91 * t152 + t92 * t153 + t110 * t362 + t400; t266 * pkin(10) + (t267 + t348) * t393 + (t307 + t349 + (-Ifges(6,1) / 0.2e1 + t399) * t393 + t261) * t394 - m(7) * (t20 * t27 + t21 * t28 + t40 * t43) + t269 + t271 * mrSges(7,3) + t326 * t43 - t42 * t103 - t27 * t77 - t28 * t76 - t396 * pkin(5); t57 - t40 * (mrSges(7,1) * t102 + mrSges(7,2) * t101) + (Ifges(7,1) * t101 - t333) * t368 + t53 * t367 + (Ifges(7,5) * t101 - Ifges(7,6) * t102) * t366 - t20 * t76 + t21 * t77 + (t101 * t20 + t102 * t21) * mrSges(7,3) + (-Ifges(7,2) * t102 + t100 + t54) * t369 + t377;];
tauc  = t1(:);
