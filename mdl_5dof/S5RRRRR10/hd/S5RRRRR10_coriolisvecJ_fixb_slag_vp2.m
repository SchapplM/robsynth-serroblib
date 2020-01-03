% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:32:00
% EndTime: 2019-12-31 22:32:43
% DurationCPUTime: 19.03s
% Computational Cost: add. (13286->677), mult. (34844->974), div. (0->0), fcn. (26946->10), ass. (0->312)
t251 = sin(qJ(4));
t252 = sin(qJ(3));
t255 = cos(qJ(4));
t256 = cos(qJ(3));
t224 = t251 * t252 - t255 * t256;
t402 = qJD(3) + qJD(4);
t170 = t402 * t224;
t257 = cos(qJ(2));
t248 = sin(pkin(5));
t318 = qJD(1) * t248;
t299 = t257 * t318;
t182 = t224 * t299;
t427 = t170 - t182;
t225 = t251 * t256 + t252 * t255;
t171 = t402 * t225;
t181 = t225 * t299;
t319 = t171 - t181;
t389 = -pkin(9) - pkin(8);
t301 = qJD(3) * t389;
t227 = t252 * t301;
t234 = t389 * t252;
t235 = t389 * t256;
t271 = t255 * t234 + t235 * t251;
t291 = t256 * t301;
t129 = qJD(4) * t271 + t255 * t227 + t251 * t291;
t253 = sin(qJ(2));
t300 = t253 * t318;
t249 = cos(pkin(5));
t317 = qJD(1) * t249;
t305 = pkin(1) * t317;
t213 = -pkin(7) * t300 + t257 * t305;
t270 = t248 * (pkin(2) * t253 - pkin(8) * t257);
t214 = qJD(1) * t270;
t154 = -t252 * t213 + t256 * t214;
t128 = (-pkin(9) * t256 * t257 + pkin(3) * t253) * t318 + t154;
t155 = t256 * t213 + t252 * t214;
t290 = t252 * t299;
t138 = -pkin(9) * t290 + t155;
t88 = t251 * t128 + t255 * t138;
t426 = -t88 + t129;
t240 = qJD(2) + t317;
t315 = qJD(2) * t257;
t297 = t252 * t315;
t313 = qJD(3) * t256;
t314 = qJD(3) * t252;
t165 = -t240 * t314 + (-t253 * t313 - t297) * t318;
t323 = t248 * t257;
t223 = t249 * t253 * pkin(1) + pkin(7) * t323;
t218 = t223 * qJD(2);
t206 = qJD(1) * t218;
t134 = -t165 * pkin(3) + t206;
t216 = pkin(7) * t299 + t253 * t305;
t180 = pkin(8) * t240 + t216;
t209 = (-pkin(2) * t257 - pkin(8) * t253 - pkin(1)) * t248;
t191 = qJD(1) * t209;
t135 = -t180 * t252 + t256 * t191;
t197 = t240 * t252 + t256 * t300;
t110 = -pkin(9) * t197 + t135;
t233 = qJD(3) - t299;
t101 = pkin(3) * t233 + t110;
t136 = t180 * t256 + t191 * t252;
t196 = t240 * t256 - t252 * t300;
t111 = pkin(9) * t196 + t136;
t311 = qJD(4) * t255;
t312 = qJD(4) * t251;
t296 = t256 * t315;
t164 = t240 * t313 + (-t253 * t314 + t296) * t318;
t316 = qJD(2) * t248;
t294 = qJD(1) * t316;
t289 = t253 * t294;
t215 = qJD(2) * t270;
t204 = qJD(1) * t215;
t324 = t248 * t253;
t241 = pkin(7) * t324;
t365 = pkin(1) * t257;
t222 = t249 * t365 - t241;
t217 = t222 * qJD(2);
t205 = qJD(1) * t217;
t95 = -qJD(3) * t136 + t256 * t204 - t205 * t252;
t62 = pkin(3) * t289 - pkin(9) * t164 + t95;
t94 = -t180 * t314 + t191 * t313 + t252 * t204 + t256 * t205;
t67 = pkin(9) * t165 + t94;
t15 = t101 * t311 - t111 * t312 + t251 * t62 + t255 * t67;
t12 = pkin(10) * t289 + t15;
t250 = sin(qJ(5));
t254 = cos(qJ(5));
t228 = qJD(4) + t233;
t322 = t255 * t111;
t56 = t101 * t251 + t322;
t54 = pkin(10) * t228 + t56;
t179 = -t240 * pkin(2) - t213;
t147 = -t196 * pkin(3) + t179;
t273 = t196 * t251 + t255 * t197;
t292 = t255 * t196 - t197 * t251;
t63 = -pkin(4) * t292 - pkin(10) * t273 + t147;
t24 = -t250 * t54 + t254 * t63;
t80 = qJD(4) * t292 + t164 * t255 + t165 * t251;
t81 = qJD(4) * t273 + t164 * t251 - t255 * t165;
t26 = t81 * pkin(4) - t80 * pkin(10) + t134;
t2 = qJD(5) * t24 + t12 * t254 + t250 * t26;
t25 = t250 * t63 + t254 * t54;
t3 = -qJD(5) * t25 - t12 * t250 + t254 * t26;
t287 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t360 = t80 * Ifges(5,4);
t119 = t228 * t250 + t254 * t273;
t44 = -qJD(5) * t119 - t250 * t80 + t254 * t289;
t41 = Ifges(6,6) * t44;
t118 = t228 * t254 - t250 * t273;
t43 = qJD(5) * t118 + t250 * t289 + t254 * t80;
t42 = Ifges(6,5) * t43;
t8 = Ifges(6,3) * t81 + t41 + t42;
t425 = t287 + t134 * mrSges(5,1) - t15 * mrSges(5,3) + t8 / 0.2e1 - t360 / 0.2e1;
t424 = t81 * Ifges(5,2);
t423 = -pkin(10) * t300 + t426;
t174 = pkin(3) * t290 + t216;
t422 = pkin(3) * t314 + t319 * pkin(4) + t427 * pkin(10) - t174;
t421 = -t147 * mrSges(5,1) + t56 * mrSges(5,3);
t142 = qJD(5) - t292;
t343 = Ifges(6,3) * t142;
t344 = Ifges(6,6) * t118;
t347 = Ifges(6,5) * t119;
t50 = t343 + t344 + t347;
t345 = Ifges(5,6) * t228;
t346 = Ifges(5,2) * t292;
t353 = Ifges(5,4) * t273;
t92 = t345 + t346 + t353;
t419 = t24 * mrSges(6,1) - t25 * mrSges(6,2) + t50 / 0.2e1 - t92 / 0.2e1;
t418 = t228 / 0.2e1;
t417 = t273 / 0.2e1;
t416 = t292 / 0.2e1;
t16 = -qJD(4) * t56 - t251 * t67 + t255 * t62;
t13 = -pkin(4) * t289 - t16;
t17 = -mrSges(6,1) * t44 + mrSges(6,2) * t43;
t415 = m(6) * t13 + t17;
t184 = t234 * t251 - t235 * t255;
t130 = qJD(4) * t184 + t227 * t251 - t255 * t291;
t87 = t128 * t255 - t138 * t251;
t414 = -t87 - t130;
t413 = Ifges(5,1) * t417;
t247 = -pkin(3) * t256 - pkin(2);
t168 = pkin(4) * t224 - pkin(10) * t225 + t247;
t121 = t168 * t254 - t184 * t250;
t412 = qJD(5) * t121 + t250 * t422 + t254 * t423;
t122 = t168 * t250 + t184 * t254;
t411 = -qJD(5) * t122 - t250 * t423 + t254 * t422;
t410 = pkin(4) * t300 - t414;
t409 = Ifges(5,4) * t416 + Ifges(5,5) * t418;
t280 = Ifges(6,5) * t254 - Ifges(6,6) * t250;
t350 = Ifges(6,4) * t254;
t282 = -Ifges(6,2) * t250 + t350;
t351 = Ifges(6,4) * t250;
t284 = Ifges(6,1) * t254 - t351;
t285 = mrSges(6,1) * t250 + mrSges(6,2) * t254;
t367 = t254 / 0.2e1;
t368 = -t250 / 0.2e1;
t380 = t142 / 0.2e1;
t382 = t119 / 0.2e1;
t384 = t118 / 0.2e1;
t352 = Ifges(6,4) * t119;
t51 = Ifges(6,2) * t118 + Ifges(6,6) * t142 + t352;
t116 = Ifges(6,4) * t118;
t52 = Ifges(6,1) * t119 + Ifges(6,5) * t142 + t116;
t326 = t111 * t251;
t55 = t101 * t255 - t326;
t53 = -pkin(4) * t228 - t55;
t408 = t280 * t380 + t282 * t384 + t284 * t382 + t53 * t285 + t52 * t367 + t51 * t368;
t390 = t413 + t409;
t332 = t24 * t254;
t278 = t25 * t250 + t332;
t407 = t278 * mrSges(6,3);
t208 = pkin(8) * t249 + t223;
t152 = -t252 * t208 + t256 * t209;
t220 = t249 * t252 + t256 * t324;
t115 = -pkin(3) * t323 - t220 * pkin(9) + t152;
t153 = t256 * t208 + t252 * t209;
t219 = t249 * t256 - t252 * t324;
t127 = pkin(9) * t219 + t153;
t406 = t251 * t115 + t255 * t127;
t405 = -t252 * t95 + t256 * t94;
t404 = -t24 * t250 + t25 * t254;
t403 = -t147 * mrSges(5,2) + t55 * mrSges(5,3);
t401 = t390 + t408 + t409;
t338 = t197 * Ifges(4,4);
t132 = t196 * Ifges(4,2) + t233 * Ifges(4,6) + t338;
t192 = Ifges(4,4) * t196;
t133 = t197 * Ifges(4,1) + t233 * Ifges(4,5) + t192;
t274 = t135 * t256 + t136 * t252;
t354 = Ifges(4,4) * t256;
t355 = Ifges(4,4) * t252;
t366 = t256 / 0.2e1;
t369 = t233 / 0.2e1;
t372 = t197 / 0.2e1;
t374 = t196 / 0.2e1;
t400 = -t274 * mrSges(4,3) + t133 * t366 - t252 * t132 / 0.2e1 + t179 * (mrSges(4,1) * t252 + mrSges(4,2) * t256) + (-Ifges(4,2) * t252 + t354) * t374 + (Ifges(4,1) * t256 - t355) * t372 + (Ifges(4,5) * t256 - Ifges(4,6) * t252) * t369;
t103 = -qJD(3) * t153 + t256 * t215 - t217 * t252;
t173 = qJD(3) * t219 + t248 * t296;
t298 = t253 * t316;
t82 = pkin(3) * t298 - pkin(9) * t173 + t103;
t102 = -t208 * t314 + t209 * t313 + t252 * t215 + t256 * t217;
t172 = -qJD(3) * t220 - t248 * t297;
t89 = pkin(9) * t172 + t102;
t21 = -qJD(4) * t406 - t251 * t89 + t255 * t82;
t99 = pkin(4) * t273 - pkin(10) * t292;
t398 = Ifges(5,2) / 0.2e1;
t10 = t43 * Ifges(6,1) + t44 * Ifges(6,4) + t81 * Ifges(6,5);
t397 = t10 / 0.2e1;
t396 = t43 / 0.2e1;
t395 = t44 / 0.2e1;
t393 = -t51 / 0.2e1;
t392 = t81 / 0.2e1;
t387 = pkin(1) * mrSges(3,1);
t386 = pkin(1) * mrSges(3,2);
t385 = -t118 / 0.2e1;
t383 = -t119 / 0.2e1;
t381 = -t142 / 0.2e1;
t272 = t255 * t219 - t220 * t251;
t379 = t272 / 0.2e1;
t160 = t219 * t251 + t220 * t255;
t378 = t160 / 0.2e1;
t377 = t164 / 0.2e1;
t376 = t165 / 0.2e1;
t373 = -t197 / 0.2e1;
t371 = t219 / 0.2e1;
t370 = t220 / 0.2e1;
t364 = t2 * t254;
t363 = t3 * t250;
t361 = t80 * Ifges(5,1);
t359 = t81 * Ifges(5,4);
t356 = Ifges(3,4) * t253;
t349 = Ifges(3,5) * t257;
t342 = t292 * Ifges(5,6);
t341 = t273 * Ifges(5,5);
t339 = t196 * Ifges(4,6);
t337 = t197 * Ifges(4,5);
t336 = t205 * mrSges(3,2);
t335 = t228 * Ifges(5,3);
t334 = t233 * Ifges(4,3);
t331 = t240 * Ifges(3,5);
t126 = mrSges(5,1) * t228 - mrSges(5,3) * t273;
t73 = -mrSges(6,1) * t118 + mrSges(6,2) * t119;
t327 = t126 - t73;
t321 = -mrSges(3,1) * t240 - mrSges(4,1) * t196 + mrSges(4,2) * t197 + mrSges(3,3) * t300;
t310 = qJD(5) * t250;
t309 = qJD(5) * t254;
t306 = Ifges(5,5) * t80 - Ifges(5,6) * t81 + Ifges(5,3) * t289;
t302 = Ifges(4,5) * t164 + Ifges(4,6) * t165 + Ifges(4,3) * t289;
t286 = mrSges(6,1) * t254 - mrSges(6,2) * t250;
t283 = Ifges(6,1) * t250 + t350;
t281 = Ifges(6,2) * t254 + t351;
t279 = Ifges(6,5) * t250 + Ifges(6,6) * t254;
t65 = -pkin(10) * t323 + t406;
t207 = t241 + (-pkin(2) - t365) * t249;
t163 = -t219 * pkin(3) + t207;
t90 = -pkin(4) * t272 - t160 * pkin(10) + t163;
t34 = t250 * t90 + t254 * t65;
t33 = -t250 * t65 + t254 * t90;
t71 = t255 * t115 - t251 * t127;
t139 = -t250 * t160 - t254 * t323;
t269 = -t254 * t160 + t250 * t323;
t20 = t115 * t311 - t127 * t312 + t251 * t82 + t255 * t89;
t150 = -t172 * pkin(3) + t218;
t9 = t43 * Ifges(6,4) + t44 * Ifges(6,2) + t81 * Ifges(6,6);
t264 = t16 * mrSges(5,1) - t15 * mrSges(5,2) + mrSges(6,3) * t364 + qJD(5) * t408 - t13 * t286 + t250 * t397 + t279 * t392 + t281 * t395 + t283 * t396 + t9 * t367 + t306;
t263 = -t353 / 0.2e1 + t347 / 0.2e1 - t345 / 0.2e1 + t344 / 0.2e1 + t343 / 0.2e1 + t419;
t22 = mrSges(6,1) * t81 - mrSges(6,3) * t43;
t23 = -mrSges(6,2) * t81 + mrSges(6,3) * t44;
t85 = -mrSges(6,2) * t142 + mrSges(6,3) * t118;
t86 = mrSges(6,1) * t142 - mrSges(6,3) * t119;
t262 = -t86 * t309 - t85 * t310 + m(6) * (-t24 * t309 - t25 * t310 - t363 + t364) + t254 * t23 - t250 * t22;
t261 = -t263 + t421;
t258 = t413 + t401;
t236 = Ifges(3,4) * t299;
t232 = t294 * t349;
t212 = -t240 * mrSges(3,2) + mrSges(3,3) * t299;
t177 = Ifges(3,1) * t300 + t236 + t331;
t176 = Ifges(3,6) * t240 + (Ifges(3,2) * t257 + t356) * t318;
t167 = mrSges(4,1) * t233 - mrSges(4,3) * t197;
t166 = -mrSges(4,2) * t233 + mrSges(4,3) * t196;
t157 = -t182 * t254 + t250 * t300;
t156 = t182 * t250 + t254 * t300;
t149 = -mrSges(4,2) * t289 + mrSges(4,3) * t165;
t148 = mrSges(4,1) * t289 - mrSges(4,3) * t164;
t131 = t334 + t337 + t339;
t125 = -mrSges(5,2) * t228 + mrSges(5,3) * t292;
t112 = -mrSges(4,1) * t165 + mrSges(4,2) * t164;
t105 = t164 * Ifges(4,1) + t165 * Ifges(4,4) + Ifges(4,5) * t289;
t104 = t164 * Ifges(4,4) + t165 * Ifges(4,2) + Ifges(4,6) * t289;
t98 = -mrSges(5,1) * t292 + mrSges(5,2) * t273;
t97 = qJD(4) * t160 - t255 * t172 + t173 * t251;
t96 = qJD(4) * t272 + t172 * t251 + t173 * t255;
t91 = t335 + t341 + t342;
t83 = pkin(3) * t197 + t99;
t69 = -mrSges(5,2) * t289 - mrSges(5,3) * t81;
t68 = mrSges(5,1) * t289 - mrSges(5,3) * t80;
t64 = pkin(4) * t323 - t71;
t60 = t110 * t255 - t326;
t59 = t110 * t251 + t322;
t58 = qJD(5) * t269 - t250 * t96 + t254 * t298;
t57 = qJD(5) * t139 + t250 * t298 + t254 * t96;
t36 = t97 * pkin(4) - t96 * pkin(10) + t150;
t35 = mrSges(5,1) * t81 + mrSges(5,2) * t80;
t32 = Ifges(5,5) * t289 - t359 + t361;
t31 = Ifges(5,6) * t289 + t360 - t424;
t30 = t250 * t99 + t254 * t55;
t29 = -t250 * t55 + t254 * t99;
t28 = t250 * t83 + t254 * t60;
t27 = -t250 * t60 + t254 * t83;
t19 = -pkin(4) * t298 - t21;
t18 = pkin(10) * t298 + t20;
t5 = -qJD(5) * t34 - t18 * t250 + t254 * t36;
t4 = qJD(5) * t33 + t18 * t254 + t250 * t36;
t1 = [(Ifges(5,5) * t96 + Ifges(5,3) * t298) * t418 + (t139 * t2 - t24 * t57 + t25 * t58 + t269 * t3) * mrSges(6,3) + (t205 * t257 + t206 * t253 + (-t213 * t257 - t216 * t253) * qJD(2)) * t248 * mrSges(3,3) + (Ifges(6,4) * t57 + Ifges(6,2) * t58) * t384 + (-Ifges(5,4) * t417 + Ifges(6,5) * t382 - Ifges(5,2) * t416 - Ifges(5,6) * t418 + Ifges(6,6) * t384 + Ifges(6,3) * t380 + t419 - t421) * t97 + (-Ifges(6,5) * t269 + Ifges(6,6) * t139) * t392 + (Ifges(5,4) * t96 + Ifges(5,6) * t298) * t416 + m(4) * (t102 * t136 + t103 * t135 + t152 * t95 + t153 * t94 + t179 * t218 + t206 * t207) + m(3) * (t205 * t223 - t206 * t222 - t213 * t218 + t216 * t217) + m(6) * (t13 * t64 + t19 * t53 + t2 * t34 + t24 * t5 + t25 * t4 + t3 * t33) + (Ifges(4,4) * t220 + Ifges(4,2) * t219 - Ifges(4,6) * t323) * t376 + (Ifges(4,1) * t220 + Ifges(4,4) * t219 - Ifges(4,5) * t323) * t377 + t32 * t378 + t31 * t379 + t96 * t390 + (Ifges(6,1) * t57 + Ifges(6,4) * t58) * t382 + (Ifges(5,1) * t96 + Ifges(5,5) * t298) * t417 - t81 * (Ifges(5,4) * t160 - Ifges(5,6) * t323) / 0.2e1 + (-Ifges(6,4) * t269 + Ifges(6,2) * t139) * t395 + (t134 * t160 + t147 * t96 + t15 * t323 - t298 * t56) * mrSges(5,2) + t4 * t85 + t5 * t86 + t71 * t68 + t19 * t73 + t64 * t17 + t57 * t52 / 0.2e1 + t53 * (-mrSges(6,1) * t58 + mrSges(6,2) * t57) + t58 * t51 / 0.2e1 - (t306 + t302) * t323 / 0.2e1 + (Ifges(6,5) * t57 + Ifges(6,6) * t58) * t380 + t34 * t23 + t33 * t22 + ((-t222 * mrSges(3,3) + Ifges(3,5) * t249 / 0.2e1 + (-0.2e1 * t386 + 0.3e1 / 0.2e1 * Ifges(3,4) * t257) * t248) * t257 + (-t223 * mrSges(3,3) + Ifges(4,5) * t370 + Ifges(4,6) * t371 + Ifges(5,5) * t378 + Ifges(5,6) * t379 - Ifges(3,6) * t249 + (-0.2e1 * t387 - 0.3e1 / 0.2e1 * t356) * t248 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3) / 0.2e1 - Ifges(5,3) / 0.2e1) * t323) * t253) * t294 + (t232 / 0.2e1 - t206 * mrSges(3,1) - t336) * t249 + t16 * (-mrSges(5,1) * t323 - t160 * mrSges(5,3)) + t95 * (-mrSges(4,1) * t323 - t220 * mrSges(4,3)) + t94 * (mrSges(4,2) * t323 + t219 * mrSges(4,3)) + (-Ifges(6,1) * t269 + Ifges(6,4) * t139) * t396 + t321 * t218 - t176 * t298 / 0.2e1 + t136 * (-mrSges(4,2) * t298 + mrSges(4,3) * t172) + t55 * (mrSges(5,1) * t298 - mrSges(5,3) * t96) + t135 * (mrSges(4,1) * t298 - mrSges(4,3) * t173) + t20 * t125 + t21 * t126 + t139 * t9 / 0.2e1 + t150 * t98 + t152 * t148 + t153 * t149 + t163 * t35 + t102 * t166 + t103 * t167 + t172 * t132 / 0.2e1 + t173 * t133 / 0.2e1 + t179 * (-mrSges(4,1) * t172 + mrSges(4,2) * t173) + t207 * t112 + t80 * (Ifges(5,1) * t160 - Ifges(5,5) * t323) / 0.2e1 + t217 * t212 + t206 * (-mrSges(4,1) * t219 + mrSges(4,2) * t220) + (-Ifges(6,6) * t395 - Ifges(6,5) * t396 - Ifges(6,3) * t392 - t424 / 0.2e1 - t425) * t272 + (Ifges(4,5) * t173 + Ifges(4,6) * t172 + Ifges(4,3) * t298) * t369 + t105 * t370 + t104 * t371 + (Ifges(4,1) * t173 + Ifges(4,4) * t172 + Ifges(4,5) * t298) * t372 + (Ifges(4,4) * t173 + Ifges(4,2) * t172 + Ifges(4,6) * t298) * t374 + ((t131 + t91) * t253 + t257 * t177 + t240 * (-Ifges(3,6) * t253 + t349)) * t316 / 0.2e1 + m(5) * (t134 * t163 + t147 * t150 + t15 * t406 + t16 * t71 + t20 * t56 + t21 * t55) + t406 * t69 - t269 * t397 + t13 * (-mrSges(6,1) * t139 - mrSges(6,2) * t269); t410 * t73 - t336 - m(4) * (t135 * t154 + t136 * t155 + t179 * t216) - m(5) * (t147 * t174 + t55 * t87 + t56 * t88) + t232 + (Ifges(4,2) * t256 + t355) * t376 + (Ifges(4,1) * t252 + t354) * t377 + (Ifges(6,5) * t157 + Ifges(6,6) * t156 + Ifges(6,3) * t181) * t381 + (Ifges(6,1) * t157 + Ifges(6,4) * t156 + Ifges(6,5) * t181) * t383 + (Ifges(6,4) * t157 + Ifges(6,2) * t156 + Ifges(6,6) * t181) * t385 + t156 * t393 + t411 * t86 + (t121 * t3 + t122 * t2 - t13 * t271 + t24 * t411 + t25 * t412 + t410 * t53) * m(6) + t412 * t85 + (t10 * t367 + t9 * t368 - t16 * mrSges(5,3) + t32 / 0.2e1 + t361 / 0.2e1 - t359 / 0.2e1 + t134 * mrSges(5,2) + t13 * t285 + t280 * t392 + t284 * t396 + t282 * t395 + (-t2 * t250 - t254 * t3) * mrSges(6,3) + (-mrSges(6,3) * t404 + t254 * t393 + t279 * t381 + t281 * t385 + t283 * t383 + t53 * t286 + t52 * t368) * qJD(5)) * t225 + m(4) * (-pkin(2) * t206 + pkin(8) * t405) + t405 * mrSges(4,3) + (-mrSges(4,1) * t256 + mrSges(4,2) * t252 - mrSges(3,1)) * t206 + (-t252 * t148 + t256 * t149) * pkin(8) + (mrSges(5,1) * t319 - mrSges(5,2) * t427) * t147 + (-t319 * t56 + t427 * t55) * mrSges(5,3) + ((m(5) * t147 + t98) * t252 * pkin(3) + (-m(4) * t274 - t252 * t166 - t256 * t167) * pkin(8) + t400) * qJD(3) + ((-t331 / 0.2e1 - t177 / 0.2e1 - t236 / 0.2e1 + t213 * mrSges(3,3) + t318 * t386 - t400) * t257 + (t136 * mrSges(4,2) - t135 * mrSges(4,1) - t339 / 0.2e1 - t337 / 0.2e1 - t334 / 0.2e1 + t176 / 0.2e1 - t131 / 0.2e1 - t91 / 0.2e1 + t216 * mrSges(3,3) + t56 * mrSges(5,2) - t55 * mrSges(5,1) - t342 / 0.2e1 - t341 / 0.2e1 - t335 / 0.2e1 + (t387 + t356 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t257) * t318 + (t240 / 0.2e1 - qJD(2)) * Ifges(3,6) + (Ifges(4,5) * t252 + Ifges(5,5) * t225 + Ifges(4,6) * t256 - Ifges(5,6) * t224) * qJD(2) / 0.2e1) * t253) * t318 + t414 * t126 - t228 * (-Ifges(5,5) * t182 - Ifges(5,6) * t181) / 0.2e1 - t292 * (-Ifges(5,4) * t182 - Ifges(5,2) * t181) / 0.2e1 - t273 * (-Ifges(5,1) * t182 - Ifges(5,4) * t181) / 0.2e1 + (-t346 / 0.2e1 + t263) * t171 - t321 * t216 - pkin(2) * t112 + t121 * t22 + t122 * t23 - t157 * t52 / 0.2e1 - t53 * (-mrSges(6,1) * t156 + mrSges(6,2) * t157) - t155 * t166 - t154 * t167 - t174 * t98 - t181 * t50 / 0.2e1 - t25 * (-mrSges(6,2) * t181 + mrSges(6,3) * t156) - t24 * (mrSges(6,1) * t181 - mrSges(6,3) * t157) + t181 * t92 / 0.2e1 + t184 * t69 - t213 * t212 + t247 * t35 + (t42 / 0.2e1 + t41 / 0.2e1 - t31 / 0.2e1 + (Ifges(6,3) / 0.2e1 + t398) * t81 + t425) * t224 + t104 * t366 + t252 * t105 / 0.2e1 - (t17 - t68) * t271 + m(5) * (t129 * t56 - t130 * t55 + t134 * t247 + t15 * t184 + t16 * t271) + t426 * t125 + t182 * t390 - (t258 - t407) * t170; t262 * (pkin(3) * t251 + pkin(10)) - t94 * mrSges(4,2) + t95 * mrSges(4,1) + (-t258 + t403) * t292 + t302 - m(6) * (t24 * t27 + t25 * t28 + t53 * t59) + (-t197 * t98 + t251 * t69 + t255 * t68 + ((m(6) * t53 - t327) * t251 + (m(6) * t404 - t250 * t86 + t254 * t85 + t125) * t255) * qJD(4) + (t15 * t251 + t16 * t255 + 0.2e1 * t147 * t373 + (-t251 * t55 + t255 * t56) * qJD(4)) * m(5)) * pkin(3) - m(5) * (-t55 * t59 + t56 * t60) + t264 - t28 * t85 - t27 * t86 + (t346 / 0.2e1 + t261) * t273 + t327 * t59 - t60 * t125 - t135 * t166 + t136 * t167 - t179 * (mrSges(4,1) * t197 + mrSges(4,2) * t196) - t233 * (Ifges(4,5) * t196 - Ifges(4,6) * t197) / 0.2e1 + t132 * t372 + (Ifges(4,1) * t196 - t338) * t373 + (-t142 * t332 + (-t142 * t25 - t3) * t250) * mrSges(6,3) + (t135 * t196 + t136 * t197) * mrSges(4,3) - (-Ifges(4,2) * t197 + t133 + t192) * t196 / 0.2e1 + t415 * (-pkin(3) * t255 - pkin(4)); t261 * t273 + t262 * pkin(10) + ((t398 - Ifges(5,1) / 0.2e1) * t273 + t407 - t401 + t403) * t292 + t264 - t30 * t85 - t29 * t86 - m(6) * (t24 * t29 + t25 * t30 + t53 * t56) + (-qJD(5) * t278 - t363) * mrSges(6,3) + t327 * t56 - t55 * t125 - t415 * pkin(4); -t53 * (mrSges(6,1) * t119 + mrSges(6,2) * t118) + (Ifges(6,1) * t118 - t352) * t383 + t51 * t382 + (Ifges(6,5) * t118 - Ifges(6,6) * t119) * t381 - t24 * t85 + t25 * t86 + (t118 * t24 + t119 * t25) * mrSges(6,3) + t287 + t8 + (-Ifges(6,2) * t119 + t116 + t52) * t385;];
tauc = t1(:);
