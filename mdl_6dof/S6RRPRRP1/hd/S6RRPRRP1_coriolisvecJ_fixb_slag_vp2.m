% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2018-11-23 17:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:10:35
% EndTime: 2018-11-23 17:10:44
% DurationCPUTime: 8.95s
% Computational Cost: add. (13408->584), mult. (35083->752), div. (0->0), fcn. (26536->8), ass. (0->279)
t423 = Ifges(6,4) + Ifges(7,4);
t424 = Ifges(6,1) + Ifges(7,1);
t418 = Ifges(6,5) + Ifges(7,5);
t422 = Ifges(6,2) + Ifges(7,2);
t417 = Ifges(6,6) + Ifges(7,6);
t257 = qJD(2) + qJD(4);
t260 = sin(qJ(5));
t263 = cos(qJ(5));
t258 = sin(pkin(10));
t259 = cos(pkin(10));
t262 = sin(qJ(2));
t265 = cos(qJ(2));
t229 = -t258 * t262 + t259 * t265;
t218 = t229 * qJD(1);
t330 = qJD(1) * t265;
t331 = qJD(1) * t262;
t219 = -t258 * t330 - t259 * t331;
t261 = sin(qJ(4));
t264 = cos(qJ(4));
t284 = t218 * t261 - t264 * t219;
t150 = t257 * t263 - t260 * t284;
t421 = t423 * t150;
t152 = t257 * t260 + t263 * t284;
t420 = t423 * t152;
t311 = t264 * t218 + t219 * t261;
t398 = t311 * t260;
t419 = pkin(5) * t398;
t230 = t258 * t265 + t259 * t262;
t220 = t230 * qJD(2);
t204 = qJD(1) * t220;
t221 = t229 * qJD(2);
t205 = qJD(1) * t221;
t123 = qJD(4) * t284 + t264 * t204 + t205 * t261;
t122 = qJD(4) * t311 - t204 * t261 + t205 * t264;
t82 = qJD(5) * t150 + t122 * t263;
t83 = -qJD(5) * t152 - t122 * t260;
t416 = t417 * t123 + t422 * t83 + t423 * t82;
t415 = t418 * t123 + t423 * t83 + t424 * t82;
t163 = qJD(5) - t311;
t414 = t422 * t150 + t417 * t163 + t420;
t406 = t424 * t152 + t418 * t163 + t421;
t249 = pkin(2) * t259 + pkin(3);
t374 = pkin(2) * t258;
t211 = t249 * t264 - t261 * t374;
t199 = t211 * qJD(4);
t367 = -qJ(3) - pkin(7);
t240 = t367 * t262;
t234 = qJD(1) * t240;
t242 = t367 * t265;
t235 = qJD(1) * t242;
t337 = t259 * t235;
t175 = -t234 * t258 + t337;
t372 = pkin(8) * t218;
t153 = t175 - t372;
t222 = t258 * t235;
t176 = t259 * t234 + t222;
t371 = pkin(8) * t219;
t154 = t176 + t371;
t98 = t153 * t261 + t154 * t264;
t127 = pkin(4) * t284 - pkin(9) * t311;
t253 = pkin(2) * t331;
t186 = -pkin(3) * t219 + t253;
t99 = t127 + t186;
t46 = t260 * t99 + t263 * t98;
t413 = t199 * t263 - t46;
t212 = t261 * t249 + t264 * t374;
t400 = -t212 * qJD(4) - t264 * t153 + t154 * t261;
t301 = mrSges(7,1) * t260 + mrSges(7,2) * t263;
t303 = mrSges(6,1) * t260 + mrSges(6,2) * t263;
t228 = qJD(2) * pkin(2) + t234;
t171 = t259 * t228 + t222;
t142 = qJD(2) * pkin(3) + t171 + t371;
t172 = t258 * t228 - t337;
t146 = t172 + t372;
t89 = t142 * t264 - t261 * t146;
t86 = -pkin(4) * t257 - t89;
t60 = -pkin(5) * t150 + qJD(6) + t86;
t412 = t60 * t301 + t86 * t303;
t411 = t418 * t260 + t417 * t263;
t359 = Ifges(7,4) * t260;
t362 = Ifges(6,4) * t260;
t410 = t422 * t263 + t359 + t362;
t358 = Ifges(7,4) * t263;
t361 = Ifges(6,4) * t263;
t409 = t424 * t260 + t358 + t361;
t408 = qJ(6) * t398 + t263 * qJD(6);
t387 = -t150 / 0.2e1;
t386 = -t152 / 0.2e1;
t383 = -t163 / 0.2e1;
t356 = t311 * Ifges(5,2);
t407 = t356 / 0.2e1;
t209 = pkin(9) + t212;
t336 = -qJ(6) - t209;
t310 = qJD(5) * t336;
t256 = t263 * qJ(6);
t394 = pkin(5) * t284 - t256 * t311;
t45 = -t260 * t98 + t263 * t99;
t405 = -t394 - t45 + (-qJD(6) - t199) * t260 + t263 * t310;
t366 = -qJ(6) - pkin(9);
t314 = qJD(5) * t366;
t51 = t263 * t127 - t260 * t89;
t404 = -qJD(6) * t260 + t263 * t314 - t394 - t51;
t403 = t260 * t310 + t408 + t413;
t52 = t260 * t127 + t263 * t89;
t402 = t260 * t314 + t408 - t52;
t328 = qJD(5) * t260;
t325 = pkin(5) * t328;
t401 = t325 - t400 - t419;
t399 = -t98 + t199;
t335 = -mrSges(5,1) * t257 - mrSges(6,1) * t150 + mrSges(6,2) * t152 + mrSges(5,3) * t284;
t179 = t259 * t240 + t242 * t258;
t160 = -pkin(8) * t230 + t179;
t180 = t258 * t240 - t259 * t242;
t161 = pkin(8) * t229 + t180;
t114 = t160 * t261 + t161 * t264;
t111 = t263 * t114;
t174 = t229 * t261 + t230 * t264;
t251 = -pkin(2) * t265 - pkin(1);
t194 = -t229 * pkin(3) + t251;
t283 = t264 * t229 - t230 * t261;
t112 = -pkin(4) * t283 - t174 * pkin(9) + t194;
t55 = t260 * t112 + t111;
t397 = t264 * t160 - t161 * t261;
t162 = Ifges(5,4) * t311;
t332 = qJD(1) * t251;
t236 = qJD(3) + t332;
t177 = -t218 * pkin(3) + t236;
t90 = t142 * t261 + t146 * t264;
t87 = pkin(9) * t257 + t90;
t94 = -pkin(4) * t311 - pkin(9) * t284 + t177;
t37 = -t260 * t87 + t263 * t94;
t25 = -qJ(6) * t152 + t37;
t22 = pkin(5) * t163 + t25;
t38 = t260 * t94 + t263 * t87;
t26 = qJ(6) * t150 + t38;
t290 = Ifges(7,5) * t263 - Ifges(7,6) * t260;
t274 = t163 * t290;
t292 = Ifges(6,5) * t263 - Ifges(6,6) * t260;
t275 = t163 * t292;
t298 = Ifges(7,1) * t263 - t359;
t276 = t152 * t298;
t300 = Ifges(6,1) * t263 - t362;
t277 = t152 * t300;
t294 = -Ifges(7,2) * t260 + t358;
t278 = t150 * t294;
t296 = -Ifges(6,2) * t260 + t361;
t279 = t150 * t296;
t286 = t260 * t38 + t263 * t37;
t376 = -t263 / 0.2e1;
t377 = t260 / 0.2e1;
t355 = t284 * Ifges(5,1);
t395 = t162 / 0.2e1 + t355 / 0.2e1;
t266 = t286 * mrSges(6,3) + (t22 * t263 + t26 * t260) * mrSges(7,3) + t89 * mrSges(5,3) - t257 * Ifges(5,5) - t279 / 0.2e1 - t278 / 0.2e1 - t277 / 0.2e1 - t276 / 0.2e1 - t275 / 0.2e1 - t274 / 0.2e1 - t177 * mrSges(5,2) + t414 * t377 + t406 * t376 - t395 - t412;
t396 = t266 - t162 / 0.2e1;
t322 = -Ifges(7,3) / 0.2e1 - Ifges(6,3) / 0.2e1;
t323 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t324 = Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t393 = t323 * t150 + t324 * t152 - t322 * t163 + t177 * mrSges(5,1) + t37 * mrSges(6,1) + t22 * mrSges(7,1) - t38 * mrSges(6,2) - t26 * mrSges(7,2) - t90 * mrSges(5,3) - Ifges(5,4) * t284 - t257 * Ifges(5,6) - t407 - t417 * t387 - t418 * t386 - (Ifges(6,3) + Ifges(7,3)) * t383;
t392 = t82 / 0.2e1;
t391 = t83 / 0.2e1;
t390 = pkin(1) * mrSges(3,1);
t389 = pkin(1) * mrSges(3,2);
t388 = t123 / 0.2e1;
t385 = t152 / 0.2e1;
t381 = -t219 / 0.2e1;
t380 = -t220 / 0.2e1;
t379 = t221 / 0.2e1;
t375 = t263 / 0.2e1;
t373 = pkin(5) * t263;
t370 = pkin(9) * t263;
t315 = qJD(2) * t367;
t215 = qJD(3) * t265 + t262 * t315;
t192 = t215 * qJD(1);
t216 = -t262 * qJD(3) + t265 * t315;
t193 = t216 * qJD(1);
t151 = t259 * t192 + t258 * t193;
t131 = -pkin(8) * t204 + t151;
t149 = -t192 * t258 + t193 * t259;
t273 = -pkin(8) * t205 + t149;
t32 = qJD(4) * t89 + t264 * t131 + t261 * t273;
t327 = qJD(5) * t263;
t248 = qJD(2) * t253;
t178 = pkin(3) * t204 + t248;
t53 = pkin(4) * t123 - pkin(9) * t122 + t178;
t7 = t260 * t53 + t263 * t32 + t94 * t327 - t328 * t87;
t369 = t263 * t7;
t365 = Ifges(3,4) * t262;
t33 = qJD(4) * t90 + t131 * t261 - t264 * t273;
t357 = t397 * t33;
t354 = t219 * Ifges(4,4);
t351 = t260 * t37;
t349 = Ifges(3,5) * qJD(2);
t348 = Ifges(3,6) * qJD(2);
t347 = qJD(2) * mrSges(3,1);
t346 = qJD(2) * mrSges(3,2);
t341 = t174 * t260;
t159 = t259 * t215 + t258 * t216;
t329 = qJD(2) * t262;
t158 = -t215 * t258 + t259 * t216;
t137 = -pkin(8) * t221 + t158;
t138 = -pkin(8) * t220 + t159;
t47 = qJD(4) * t397 + t137 * t261 + t138 * t264;
t128 = qJD(4) * t283 - t220 * t261 + t221 * t264;
t129 = qJD(4) * t174 + t264 * t220 + t221 * t261;
t187 = pkin(2) * t329 + pkin(3) * t220;
t59 = pkin(4) * t129 - pkin(9) * t128 + t187;
t326 = t112 * t327 + t260 * t59 + t263 * t47;
t321 = t174 * t327;
t320 = t349 / 0.2e1;
t319 = -t348 / 0.2e1;
t27 = -t83 * mrSges(7,1) + t82 * mrSges(7,2);
t316 = -t260 * t47 + t263 * t59;
t313 = t204 * mrSges(4,1) + t205 * mrSges(4,2);
t312 = t123 * mrSges(5,1) + t122 * mrSges(5,2);
t54 = t263 * t112 - t114 * t260;
t8 = -qJD(5) * t38 - t260 * t32 + t263 * t53;
t1 = pkin(5) * t123 - qJ(6) * t82 - qJD(6) * t152 + t8;
t309 = -t8 * mrSges(6,3) - t1 * mrSges(7,3);
t208 = -pkin(4) - t211;
t308 = -t37 * mrSges(6,3) - t22 * mrSges(7,3);
t307 = -t38 * mrSges(6,3) - t26 * mrSges(7,3);
t3 = qJ(6) * t83 + qJD(6) * t150 + t7;
t306 = -t1 * t263 - t260 * t3;
t305 = -t260 * t7 - t263 * t8;
t304 = mrSges(6,1) * t263 - mrSges(6,2) * t260;
t302 = mrSges(7,1) * t263 - mrSges(7,2) * t260;
t287 = t22 * t260 - t26 * t263;
t285 = -t263 * t38 + t351;
t282 = -qJ(6) * t128 - qJD(6) * t174;
t272 = t8 * mrSges(6,1) + t1 * mrSges(7,1) - t7 * mrSges(6,2) - t3 * mrSges(7,2);
t48 = qJD(4) * t114 - t264 * t137 + t138 * t261;
t271 = m(6) * (-qJD(5) * t286 - t8 * t260 + t369);
t12 = -pkin(5) * t83 + t33;
t269 = t3 * t263 * mrSges(7,3) - t32 * mrSges(5,2) + mrSges(6,3) * t369 + Ifges(5,5) * t122 - Ifges(5,6) * t123 - t12 * t302 + t409 * t392 + t410 * t391 + t411 * t388 + t415 * t377 + t416 * t375 + (-t304 - mrSges(5,1)) * t33 - t414 * t328 / 0.2e1 + t406 * t327 / 0.2e1 + t412 * qJD(5) + (t279 + t278 + t277 + t276 + t275 + t274) * qJD(5) / 0.2e1;
t252 = Ifges(3,4) * t330;
t250 = -pkin(4) - t373;
t241 = t256 + t370;
t239 = t366 * t260;
t238 = mrSges(3,3) * t330 - t346;
t237 = -mrSges(3,3) * t331 + t347;
t227 = Ifges(3,1) * t331 + t252 + t349;
t226 = t348 + (t265 * Ifges(3,2) + t365) * qJD(1);
t210 = Ifges(4,4) * t218;
t191 = qJD(2) * mrSges(4,1) + t219 * mrSges(4,3);
t190 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t218;
t188 = t208 - t373;
t182 = t209 * t263 + t256;
t181 = t336 * t260;
t170 = -mrSges(4,1) * t218 - mrSges(4,2) * t219;
t165 = -t219 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t210;
t164 = t218 * Ifges(4,2) + Ifges(4,6) * qJD(2) - t354;
t155 = -mrSges(5,2) * t257 + mrSges(5,3) * t311;
t126 = -mrSges(5,1) * t311 + mrSges(5,2) * t284;
t119 = Ifges(6,3) * t123;
t118 = Ifges(7,3) * t123;
t108 = mrSges(6,1) * t163 - mrSges(6,3) * t152;
t107 = mrSges(7,1) * t163 - mrSges(7,3) * t152;
t106 = -mrSges(6,2) * t163 + mrSges(6,3) * t150;
t105 = -mrSges(7,2) * t163 + mrSges(7,3) * t150;
t100 = -mrSges(7,1) * t150 + mrSges(7,2) * t152;
t85 = pkin(5) * t341 - t397;
t81 = Ifges(6,5) * t82;
t80 = Ifges(7,5) * t82;
t79 = Ifges(6,6) * t83;
t78 = Ifges(7,6) * t83;
t61 = t90 + t419;
t43 = -qJ(6) * t341 + t55;
t42 = -mrSges(6,2) * t123 + mrSges(6,3) * t83;
t41 = -mrSges(7,2) * t123 + mrSges(7,3) * t83;
t40 = mrSges(6,1) * t123 - mrSges(6,3) * t82;
t39 = mrSges(7,1) * t123 - mrSges(7,3) * t82;
t35 = -pkin(5) * t283 - t174 * t256 + t54;
t28 = -mrSges(6,1) * t83 + mrSges(6,2) * t82;
t21 = (t128 * t260 + t321) * pkin(5) + t48;
t10 = -qJD(5) * t55 + t316;
t9 = -t114 * t328 + t326;
t5 = -qJ(6) * t321 + (-qJD(5) * t114 + t282) * t260 + t326;
t4 = pkin(5) * t129 + t282 * t263 + (-t111 + (qJ(6) * t174 - t112) * t260) * qJD(5) + t316;
t2 = [(-t114 * t123 - t122 * t397) * mrSges(5,3) - t397 * t28 + (-t266 + t395) * t128 + (Ifges(5,1) * t122 - Ifges(5,4) * t123 + t178 * mrSges(5,2) + t12 * t301 + (mrSges(5,3) + t303) * t33 + t306 * mrSges(7,3) + t305 * mrSges(6,3) + (mrSges(6,3) * t285 + mrSges(7,3) * t287 + t302 * t60 + t304 * t86 + t414 * t376 + t411 * t383 + t409 * t386 + t410 * t387) * qJD(5) + (t298 + t300) * t392 + (t296 + t294) * t391 + (t290 + t292) * t388 + t415 * t375 - (qJD(5) * t406 + t416) * t260 / 0.2e1) * t174 + (-t356 / 0.2e1 + t393) * t129 - (t80 / 0.2e1 + t78 / 0.2e1 + t118 / 0.2e1 + t81 / 0.2e1 + t79 / 0.2e1 + t119 / 0.2e1 - Ifges(5,4) * t122 + t178 * mrSges(5,1) - t32 * mrSges(5,3) + t323 * t83 + t324 * t82 + (Ifges(5,2) - t322) * t123 + t272) * t283 + (-t149 * t230 + t151 * t229 - t171 * t221 - t172 * t220 - t179 * t205 - t180 * t204) * mrSges(4,3) + (-t230 * t204 + t229 * t205 + t218 * t379 - t220 * t381) * Ifges(4,4) + (-t229 * t204 + t218 * t380) * Ifges(4,2) + m(6) * (t10 * t37 + t38 * t9 + t48 * t86 + t54 * t8 + t55 * t7 - t357) + m(5) * (t114 * t32 + t177 * t187 + t178 * t194 + t47 * t90 - t48 * t89 - t357) + t236 * (mrSges(4,1) * t220 + mrSges(4,2) * t221) + m(7) * (t1 * t35 + t12 * t85 + t21 * t60 + t22 * t4 + t26 * t5 + t3 * t43) + t54 * t40 + t55 * t42 + m(4) * (t149 * t179 + t151 * t180 + t158 * t171 + t159 * t172) + t43 * t41 + t35 * t39 + t165 * t379 + t164 * t380 + t194 * t312 + t251 * t313 + t85 * t27 + t21 * t100 + t5 * t105 + t9 * t106 + t4 * t107 + t10 * t108 + t335 * t48 + t47 * t155 + t187 * t126 + t159 * t190 + t158 * t191 + (t230 * t205 + t221 * t381) * Ifges(4,1) + (Ifges(4,5) * t379 + Ifges(4,6) * t380 + (t227 / 0.2e1 - pkin(7) * t237 + t320 + (-0.2e1 * t389 + 0.3e1 / 0.2e1 * Ifges(3,4) * t265) * qJD(1)) * t265) * qJD(2) + (-t226 / 0.2e1 - pkin(7) * t238 + t319 + (-0.2e1 * t390 - 0.3e1 / 0.2e1 * t365 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t265) * qJD(1) + (m(4) * (t236 + t332) + qJD(1) * (-mrSges(4,1) * t229 + mrSges(4,2) * t230) + t170) * pkin(2)) * t329; (t407 - t393) * t284 + t399 * t155 + (-t177 * t186 - t211 * t33 + t212 * t32 + t399 * t90 + t400 * t89) * m(5) - t335 * t400 + t401 * t100 + t403 * t105 + (t1 * t181 + t12 * t188 + t182 * t3 + t22 * t405 + t26 * t403 + t401 * t60) * m(7) + t405 * t107 + (-t355 / 0.2e1 + t396) * t311 + (-t199 * t351 + t208 * t33 - t37 * t45 + t38 * t413 - t400 * t86) * m(6) + (t171 * t218 - t172 * t219 + (-t204 * t258 - t205 * t259) * pkin(2)) * mrSges(4,3) + (-t122 * t211 - t123 * t212) * mrSges(5,3) + t269 + (-t199 * t108 - t209 * t40 + (-t106 * t209 + t307) * qJD(5) + t309) * t260 - (Ifges(4,2) * t219 + t165 + t210) * t218 / 0.2e1 + (t199 * t106 + t209 * t42 + (-t108 * t209 + t308) * qJD(5)) * t263 - m(4) * (t171 * t175 + t172 * t176) + t164 * t381 - t46 * t106 - t45 * t108 + m(4) * (t149 * t259 + t151 * t258) * pkin(2) + t149 * mrSges(4,1) + t209 * t271 - t151 * mrSges(4,2) + t219 * (Ifges(4,1) * t218 + t354) / 0.2e1 + t181 * t39 + t182 * t41 - t186 * t126 + t188 * t27 - t176 * t190 - t175 * t191 - Ifges(4,6) * t204 + Ifges(4,5) * t205 + t208 * t28 + ((-t252 / 0.2e1 - t227 / 0.2e1 + t320 + qJD(1) * t389 + (t237 - t347) * pkin(7)) * t265 + (t226 / 0.2e1 + t319 + (t390 + t365 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t265) * qJD(1) + (t238 + t346) * pkin(7) + (-m(4) * t236 - t170) * pkin(2)) * t262) * qJD(1) - qJD(2) * (Ifges(4,5) * t218 + Ifges(4,6) * t219) / 0.2e1 - t236 * (-mrSges(4,1) * t219 + mrSges(4,2) * t218); -t311 * t155 - t218 * t190 - t219 * t191 + (-t100 - t335) * t284 + (t39 + t40 + t163 * (t105 + t106)) * t263 + (t41 + t42 - t163 * (t107 + t108)) * t260 + t312 + t313 + (-t163 * t287 - t284 * t60 - t306) * m(7) + (-t163 * t285 - t284 * t86 - t305) * m(6) + (t284 * t89 - t311 * t90 + t178) * m(5) + (-t171 * t219 - t172 * t218 + t248) * m(4); pkin(9) * t271 + (-pkin(9) * t40 + t309) * t260 + ((Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t284 + t396) * t311 - t335 * t90 + t402 * t105 + t269 + t404 * t107 - pkin(4) * t28 + ((-pkin(9) * t108 + t308) * t263 + (pkin(5) * t100 - pkin(9) * t106 + t307) * t260) * qJD(5) + t42 * t370 - t61 * t100 - t52 * t106 - t51 * t108 - t89 * t155 - t393 * t284 + t239 * t39 + t241 * t41 + t250 * t27 + (t1 * t239 + t12 * t250 + t241 * t3 + (t325 - t61) * t60 + t402 * t26 + t404 * t22) * m(7) + (-pkin(4) * t33 - t37 * t51 - t38 * t52 - t86 * t90) * m(6); (-t100 * t152 + t39) * pkin(5) + (t150 * t22 + t152 * t26) * mrSges(7,3) + (t150 * t37 + t152 * t38) * mrSges(6,3) + t81 + t80 + t79 + t78 + t272 + t118 - t25 * t105 - t37 * t106 + t26 * t107 + t38 * t108 + t119 - t60 * (mrSges(7,1) * t152 + mrSges(7,2) * t150) - t86 * (mrSges(6,1) * t152 + mrSges(6,2) * t150) + (-(-t22 + t25) * t26 + (-t152 * t60 + t1) * pkin(5)) * m(7) + (t424 * t150 - t420) * t386 + t414 * t385 + (t150 * t418 - t152 * t417) * t383 + (-t422 * t152 + t406 + t421) * t387; -t150 * t105 + t152 * t107 + 0.2e1 * (t12 / 0.2e1 + t26 * t387 + t22 * t385) * m(7) + t27;];
tauc  = t2(:);
