% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:51:14
% EndTime: 2019-03-09 03:51:37
% DurationCPUTime: 12.37s
% Computational Cost: add. (15367->605), mult. (40234->832), div. (0->0), fcn. (31863->10), ass. (0->256)
t255 = cos(pkin(10));
t335 = cos(qJ(3));
t279 = t335 * t255;
t246 = qJD(1) * t279;
t253 = sin(pkin(10));
t258 = sin(qJ(3));
t295 = t258 * t253;
t222 = qJD(1) * t295 - t246;
t252 = sin(pkin(11));
t257 = sin(qJ(5));
t254 = cos(pkin(11));
t260 = cos(qJ(5));
t297 = t254 * t260;
t267 = t252 * t257 - t297;
t169 = t267 * t222;
t224 = t267 * qJD(5);
t394 = t224 + t169;
t235 = t252 * t260 + t254 * t257;
t168 = t235 * t222;
t225 = t235 * qJD(5);
t289 = -t225 - t168;
t236 = t253 * t335 + t258 * t255;
t223 = t236 * qJD(1);
t186 = pkin(3) * t223 + qJ(4) * t222;
t327 = pkin(7) + qJ(2);
t243 = t327 * t255;
t238 = qJD(1) * t243;
t220 = t258 * t238;
t241 = t327 * t253;
t237 = qJD(1) * t241;
t280 = t335 * t237;
t190 = -t220 - t280;
t124 = t252 * t186 + t254 * t190;
t300 = t222 * t252;
t106 = pkin(8) * t300 + t124;
t326 = pkin(8) + qJ(4);
t240 = t326 * t252;
t242 = t326 * t254;
t286 = qJD(5) * t260;
t123 = t254 * t186 - t190 * t252;
t332 = pkin(8) * t254;
t92 = pkin(4) * t223 + t222 * t332 + t123;
t371 = qJD(4) * t297 - t260 * t106 - t240 * t286 + (-qJD(4) * t252 - qJD(5) * t242 - t92) * t257;
t195 = -t257 * t240 + t260 * t242;
t370 = -t235 * qJD(4) - qJD(5) * t195 + t106 * t257 - t260 * t92;
t199 = qJD(3) * t252 + t223 * t254;
t275 = t254 * qJD(3) - t223 * t252;
t146 = t199 * t260 + t257 * t275;
t256 = sin(qJ(6));
t259 = cos(qJ(6));
t384 = -t199 * t257 + t260 * t275;
t393 = -t146 * t256 + t259 * t384;
t86 = t146 * t259 + t256 * t384;
t392 = pkin(9) * t289 + t371;
t391 = -pkin(5) * t223 + pkin(9) * t394 + t370;
t281 = -pkin(2) * t255 - pkin(1);
t239 = qJD(1) * t281 + qJD(2);
t160 = pkin(3) * t222 - qJ(4) * t223 + t239;
t191 = -t258 * t237 + t335 * t238;
t185 = qJD(3) * qJ(4) + t191;
t112 = t254 * t160 - t185 * t252;
t73 = pkin(4) * t222 - pkin(8) * t199 + t112;
t113 = t252 * t160 + t254 * t185;
t91 = pkin(8) * t275 + t113;
t43 = t257 * t73 + t260 * t91;
t34 = pkin(9) * t384 + t43;
t309 = t256 * t34;
t218 = qJD(5) + t222;
t42 = -t257 * t91 + t260 * t73;
t33 = -pkin(9) * t146 + t42;
t32 = pkin(5) * t218 + t33;
t11 = t259 * t32 - t309;
t308 = t259 * t34;
t12 = t256 * t32 + t308;
t227 = t236 * qJD(3);
t214 = qJD(1) * t227;
t265 = t279 - t295;
t226 = t265 * qJD(3);
t213 = qJD(1) * t226;
t96 = qJD(5) * t384 - t213 * t267;
t97 = -qJD(5) * t146 - t213 * t235;
t30 = qJD(6) * t393 + t256 * t97 + t259 * t96;
t31 = -qJD(6) * t86 - t256 * t96 + t259 * t97;
t284 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t214;
t334 = Ifges(7,4) * t86;
t212 = qJD(6) + t218;
t343 = -t212 / 0.2e1;
t354 = -t86 / 0.2e1;
t356 = -t393 / 0.2e1;
t287 = qJD(5) * t257;
t302 = t213 * t254;
t142 = pkin(3) * t214 - qJ(4) * t213 - qJD(4) * t223;
t285 = qJD(1) * qJD(2);
t278 = t253 * t285;
t291 = qJD(2) * t246 - qJD(3) * t280;
t148 = -t258 * t278 + (qJD(4) - t220) * qJD(3) + t291;
t78 = t254 * t142 - t148 * t252;
t59 = pkin(4) * t214 - pkin(8) * t302 + t78;
t303 = t213 * t252;
t79 = t252 * t142 + t254 * t148;
t62 = -pkin(8) * t303 + t79;
t15 = t257 * t59 + t260 * t62 + t73 * t286 - t287 * t91;
t10 = pkin(9) * t97 + t15;
t16 = -qJD(5) * t43 - t257 * t62 + t260 * t59;
t8 = pkin(5) * t214 - pkin(9) * t96 + t16;
t2 = qJD(6) * t11 + t10 * t259 + t256 * t8;
t3 = -qJD(6) * t12 - t10 * t256 + t259 * t8;
t378 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t80 = Ifges(7,4) * t393;
t38 = Ifges(7,1) * t86 + Ifges(7,5) * t212 + t80;
t183 = -qJD(3) * pkin(3) + qJD(4) - t190;
t135 = -pkin(4) * t275 + t183;
t81 = -pkin(5) * t384 + t135;
t388 = t284 + t378 + (Ifges(7,5) * t393 - Ifges(7,6) * t86) * t343 + (t11 * t393 + t12 * t86) * mrSges(7,3) + (-Ifges(7,2) * t86 + t38 + t80) * t356 - t81 * (mrSges(7,1) * t86 + mrSges(7,2) * t393) + (Ifges(7,1) * t393 - t334) * t354;
t387 = Ifges(4,1) / 0.2e1;
t386 = -t222 / 0.2e1;
t385 = Ifges(5,3) + Ifges(4,2);
t360 = t30 / 0.2e1;
t359 = t31 / 0.2e1;
t37 = Ifges(7,2) * t393 + Ifges(7,6) * t212 + t334;
t382 = t37 / 0.2e1;
t352 = t96 / 0.2e1;
t351 = t97 / 0.2e1;
t341 = t214 / 0.2e1;
t381 = Ifges(4,4) * t386;
t380 = t222 / 0.2e1;
t379 = t223 * t387;
t193 = -t260 * t240 - t242 * t257;
t172 = -pkin(9) * t235 + t193;
t173 = -pkin(9) * t267 + t195;
t110 = t172 * t259 - t173 * t256;
t377 = qJD(6) * t110 + t391 * t256 + t259 * t392;
t111 = t172 * t256 + t173 * t259;
t376 = -qJD(6) * t111 - t256 * t392 + t391 * t259;
t154 = -pkin(4) * t300 + t191;
t367 = -pkin(5) * t289 - t154;
t108 = t168 * t256 + t169 * t259;
t188 = -t235 * t256 - t259 * t267;
t125 = qJD(6) * t188 - t224 * t259 - t225 * t256;
t294 = t125 - t108;
t107 = t168 * t259 - t169 * t256;
t189 = t235 * t259 - t256 * t267;
t126 = -qJD(6) * t189 + t224 * t256 - t225 * t259;
t293 = t126 - t107;
t366 = -t335 * t241 - t258 * t243;
t365 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t364 = (m(3) * qJ(2) + mrSges(3,3)) * (t253 ^ 2 + t255 ^ 2);
t320 = Ifges(5,2) * t252;
t323 = Ifges(5,4) * t254;
t271 = -t320 + t323;
t324 = Ifges(5,4) * t252;
t272 = t254 * Ifges(5,1) - t324;
t273 = mrSges(5,1) * t252 + mrSges(5,2) * t254;
t336 = t254 / 0.2e1;
t337 = -t252 / 0.2e1;
t363 = t183 * t273 + t271 * t275 / 0.2e1 + t272 * t199 / 0.2e1 + t239 * mrSges(4,2) + (-t112 * t254 - t113 * t252) * mrSges(5,3) + t381 + Ifges(4,5) * qJD(3) + t379 + (t199 * Ifges(5,4) + Ifges(5,2) * t275 + Ifges(5,6) * t222) * t337 + (t199 * Ifges(5,1) + Ifges(5,4) * t275 + Ifges(5,5) * t222) * t336;
t362 = Ifges(7,4) * t360 + Ifges(7,2) * t359 + Ifges(7,6) * t341;
t361 = Ifges(7,1) * t360 + Ifges(7,4) * t359 + Ifges(7,5) * t341;
t358 = Ifges(6,4) * t352 + Ifges(6,2) * t351 + Ifges(6,6) * t341;
t357 = Ifges(6,1) * t352 + Ifges(6,4) * t351 + Ifges(6,5) * t341;
t355 = t393 / 0.2e1;
t353 = t86 / 0.2e1;
t349 = -t384 / 0.2e1;
t348 = t384 / 0.2e1;
t347 = -t146 / 0.2e1;
t346 = t146 / 0.2e1;
t342 = t212 / 0.2e1;
t340 = -t218 / 0.2e1;
t339 = t218 / 0.2e1;
t325 = Ifges(4,4) * t223;
t322 = Ifges(6,4) * t146;
t321 = Ifges(5,5) * t254;
t319 = Ifges(5,6) * t252;
t313 = t214 * Ifges(5,5);
t312 = t214 * Ifges(5,6);
t187 = -pkin(3) * t265 - qJ(4) * t236 + t281;
t196 = -t258 * t241 + t243 * t335;
t130 = t254 * t187 - t196 * t252;
t103 = -pkin(4) * t265 - t236 * t332 + t130;
t131 = t252 * t187 + t254 * t196;
t298 = t236 * t252;
t114 = -pkin(8) * t298 + t131;
t55 = t257 * t103 + t260 * t114;
t264 = t236 * qJD(2);
t153 = qJD(1) * t264 + qJD(3) * t191;
t304 = t153 * t366;
t299 = t226 * t252;
t155 = pkin(3) * t227 - qJ(4) * t226 - qJD(4) * t236;
t163 = t265 * qJD(2) + qJD(3) * t366;
t102 = t252 * t155 + t254 * t163;
t159 = mrSges(5,1) * t303 + mrSges(5,2) * t302;
t292 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t275 - mrSges(5,2) * t199 - mrSges(4,3) * t223;
t283 = Ifges(6,5) * t96 + Ifges(6,6) * t97 + Ifges(6,3) * t214;
t282 = Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t248 = -pkin(4) * t254 - pkin(3);
t51 = -t97 * mrSges(6,1) + t96 * mrSges(6,2);
t9 = -t31 * mrSges(7,1) + t30 * mrSges(7,2);
t277 = t214 * mrSges(4,1) + t213 * mrSges(4,2);
t54 = t260 * t103 - t114 * t257;
t101 = t254 * t155 - t163 * t252;
t167 = pkin(4) * t298 - t366;
t270 = -t319 + t321;
t269 = t252 * t79 + t254 * t78;
t178 = t267 * t236;
t44 = -pkin(5) * t265 + pkin(9) * t178 + t54;
t177 = t235 * t236;
t48 = -pkin(9) * t177 + t55;
t21 = -t256 * t48 + t259 * t44;
t22 = t256 * t44 + t259 * t48;
t268 = t112 * t252 - t113 * t254;
t117 = -t177 * t259 + t178 * t256;
t118 = -t177 * t256 - t178 * t259;
t67 = pkin(4) * t227 - t226 * t332 + t101;
t77 = -pkin(8) * t299 + t102;
t23 = t103 * t286 - t114 * t287 + t257 * t67 + t260 * t77;
t24 = -qJD(5) * t55 - t257 * t77 + t260 * t67;
t164 = qJD(3) * t196 + t264;
t134 = pkin(4) * t299 + t164;
t122 = pkin(4) * t303 + t153;
t262 = t11 * mrSges(7,1) + t112 * mrSges(5,1) + t239 * mrSges(4,1) + t42 * mrSges(6,1) + t199 * Ifges(5,5) + t275 * Ifges(5,6) - Ifges(4,6) * qJD(3) - t325 / 0.2e1 + t212 * Ifges(7,3) + t86 * Ifges(7,5) + t393 * Ifges(7,6) + t218 * Ifges(6,3) + t146 * Ifges(6,5) + t384 * Ifges(6,6) - t113 * mrSges(5,2) - t12 * mrSges(7,2) - t43 * mrSges(6,2) + t385 * t380;
t206 = pkin(5) * t267 + t248;
t204 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t222;
t171 = mrSges(5,1) * t214 - mrSges(5,3) * t302;
t170 = -mrSges(5,2) * t214 - mrSges(5,3) * t303;
t166 = mrSges(5,1) * t222 - mrSges(5,3) * t199;
t165 = -mrSges(5,2) * t222 + mrSges(5,3) * t275;
t152 = (-qJD(3) * t238 - t278) * t258 + t291;
t143 = Ifges(6,4) * t384;
t133 = t213 * t272 + t313;
t132 = t213 * t271 + t312;
t121 = mrSges(6,1) * t218 - mrSges(6,3) * t146;
t120 = -mrSges(6,2) * t218 + mrSges(6,3) * t384;
t119 = pkin(5) * t177 + t167;
t116 = t224 * t236 - t226 * t235;
t115 = -t225 * t236 - t226 * t267;
t88 = -mrSges(6,1) * t384 + mrSges(6,2) * t146;
t75 = -mrSges(6,2) * t214 + mrSges(6,3) * t97;
t74 = mrSges(6,1) * t214 - mrSges(6,3) * t96;
t71 = t146 * Ifges(6,1) + t218 * Ifges(6,5) + t143;
t70 = Ifges(6,2) * t384 + t218 * Ifges(6,6) + t322;
t64 = mrSges(7,1) * t212 - mrSges(7,3) * t86;
t63 = -mrSges(7,2) * t212 + mrSges(7,3) * t393;
t60 = -t116 * pkin(5) + t134;
t56 = -t97 * pkin(5) + t122;
t45 = -mrSges(7,1) * t393 + mrSges(7,2) * t86;
t40 = -qJD(6) * t118 - t115 * t256 + t116 * t259;
t39 = qJD(6) * t117 + t115 * t259 + t116 * t256;
t26 = -mrSges(7,2) * t214 + mrSges(7,3) * t31;
t25 = mrSges(7,1) * t214 - mrSges(7,3) * t30;
t18 = pkin(9) * t116 + t23;
t17 = pkin(5) * t227 - pkin(9) * t115 + t24;
t14 = t259 * t33 - t309;
t13 = -t256 * t33 - t308;
t5 = -qJD(6) * t22 + t17 * t259 - t18 * t256;
t4 = qJD(6) * t21 + t17 * t256 + t18 * t259;
t1 = [(Ifges(7,1) * t118 + Ifges(7,4) * t117) * t360 + m(5) * (t101 * t112 + t102 * t113 + t130 * t78 + t131 * t79 + t164 * t183 - t304) + m(4) * (t152 * t196 + t163 * t191 - t164 * t190 - t304) + (-t11 * t39 + t117 * t2 - t118 * t3 + t12 * t40) * mrSges(7,3) + 0.2e1 * t364 * t285 + (Ifges(7,4) * t118 + Ifges(7,2) * t117) * t359 + (-Ifges(6,4) * t178 - Ifges(6,2) * t177) * t351 + (-Ifges(6,1) * t178 - Ifges(6,4) * t177) * t352 + (-t115 * t42 + t116 * t43 - t15 * t177 + t16 * t178) * mrSges(6,3) + (-Ifges(6,5) * t178 + Ifges(7,5) * t118 - Ifges(6,6) * t177 + Ifges(7,6) * t117) * t341 + t122 * (mrSges(6,1) * t177 - mrSges(6,2) * t178) + (Ifges(7,1) * t39 + Ifges(7,4) * t40) * t353 + (Ifges(7,4) * t39 + Ifges(7,2) * t40) * t355 - t178 * t357 - t177 * t358 + t118 * t361 + t117 * t362 + (-t190 * t226 - t191 * t227 - t196 * t214 - t213 * t366) * mrSges(4,3) - t366 * t159 + (t270 * t380 + t363 + t379) * t226 + m(7) * (t11 * t5 + t119 * t56 + t12 * t4 + t2 * t22 + t21 * t3 + t60 * t81) + m(6) * (t122 * t167 + t134 * t135 + t15 * t55 + t16 * t54 + t23 * t43 + t24 * t42) + t163 * t204 + t131 * t170 + t130 * t171 + t102 * t165 + t101 * t166 + t167 * t51 + t134 * t88 + t135 * (-mrSges(6,1) * t116 + mrSges(6,2) * t115) + (-t214 * Ifges(4,4) + (Ifges(4,1) + Ifges(5,1) * t254 ^ 2 / 0.2e1 + (-t323 + t320 / 0.2e1) * t252) * t213 - mrSges(5,3) * t269 + t132 * t337 + t133 * t336 + t270 * t341 + (mrSges(4,3) + t273) * t153) * t236 + t56 * (-mrSges(7,1) * t117 + mrSges(7,2) * t118) + t119 * t9 + t23 * t120 + t24 * t121 + t115 * t71 / 0.2e1 + t116 * t70 / 0.2e1 + t81 * (-mrSges(7,1) * t40 + mrSges(7,2) * t39) + t54 * t74 + t55 * t75 + t60 * t45 + t4 * t63 + t5 * t64 + t39 * t38 / 0.2e1 + t21 * t25 + t22 * t26 + (-t385 * t214 - t284 / 0.2e1 - t283 / 0.2e1 - mrSges(5,1) * t78 + mrSges(5,2) * t79 + mrSges(4,3) * t152 - Ifges(6,5) * t352 - Ifges(7,5) * t360 - Ifges(6,6) * t351 - Ifges(7,6) * t359 - (Ifges(6,3) + Ifges(7,3)) * t341 - (-Ifges(4,4) + t270) * t213 - t365 - t378) * t265 + t281 * t277 + (t222 * t282 + t262) * t227 + (-t223 * t227 / 0.2e1 + t226 * t386) * Ifges(4,4) - t292 * t164 + t40 * t382 + (Ifges(6,5) * t115 + Ifges(6,6) * t116) * t339 + (Ifges(7,5) * t39 + Ifges(7,6) * t40) * t342 + (Ifges(6,1) * t115 + Ifges(6,4) * t116) * t346 + (Ifges(6,4) * t115 + Ifges(6,2) * t116) * t348; t289 * t121 - t394 * t120 + (-t45 - t88 + t292) * t223 + (t165 * t254 - t166 * t252 + t204) * t222 - m(4) * (-t190 * t223 - t191 * t222) + t293 * t64 + t294 * t63 + t254 * t171 + t252 * t170 + t235 * t75 - t267 * t74 + t189 * t26 + t188 * t25 + t277 - t364 * qJD(1) ^ 2 + (t11 * t293 + t12 * t294 + t188 * t3 + t189 * t2 - t223 * t81) * m(7) + (-t135 * t223 + t15 * t235 - t16 * t267 + t289 * t42 - t394 * t43) * m(6) + (-t183 * t223 - t222 * t268 + t269) * m(5); (t125 / 0.2e1 - t108 / 0.2e1) * t38 + (t126 / 0.2e1 - t107 / 0.2e1) * t37 + t370 * t121 + t371 * t120 + (t122 * t248 - t135 * t154 + t15 * t195 + t16 * t193 + t370 * t42 + t371 * t43) * m(6) + t367 * t45 + (-t112 * t123 - t113 * t124 - t183 * t191 - pkin(3) * t153 - t268 * qJD(4) + (-t252 * t78 + t254 * t79) * qJ(4)) * m(5) + (-t15 * t267 - t16 * t235 + t289 * t43 + t394 * t42) * mrSges(6,3) + (-mrSges(6,1) * t289 - mrSges(6,2) * t394) * t135 + (Ifges(7,1) * t125 + Ifges(7,4) * t126) * t353 + (Ifges(7,1) * t108 + Ifges(7,4) * t107) * t354 + (Ifges(7,4) * t125 + Ifges(7,2) * t126) * t355 + (Ifges(7,4) * t108 + Ifges(7,2) * t107) * t356 + t235 * t357 + (Ifges(7,4) * t189 + Ifges(7,2) * t188) * t359 + (Ifges(7,1) * t189 + Ifges(7,4) * t188) * t360 + t189 * t361 + t188 * t362 + (-t224 / 0.2e1 - t169 / 0.2e1) * t71 + (-Ifges(6,5) * t224 - Ifges(6,6) * t225) * t339 + (-Ifges(6,1) * t224 - Ifges(6,4) * t225) * t346 + (-Ifges(6,4) * t224 - Ifges(6,2) * t225) * t348 + (-t225 / 0.2e1 - t168 / 0.2e1) * t70 - t267 * t358 + (Ifges(6,5) * t235 + Ifges(7,5) * t189 - Ifges(6,6) * t267 + Ifges(7,6) * t188) * t341 + (Ifges(6,4) * t235 - Ifges(6,2) * t267) * t351 + (Ifges(6,1) * t235 - Ifges(6,4) * t267) * t352 + t122 * (mrSges(6,1) * t267 + mrSges(6,2) * t235) + t376 * t64 + t377 * t63 + (t11 * t376 + t110 * t3 + t111 * t2 + t12 * t377 + t206 * t56 + t367 * t81) * m(7) + t248 * t51 - Ifges(4,6) * t214 - t190 * t204 + t206 * t9 + t193 * t74 + t195 * t75 + t56 * (-mrSges(7,1) * t188 + mrSges(7,2) * t189) - t124 * t165 - t123 * t166 - pkin(3) * t159 - t152 * mrSges(4,2) - t153 * mrSges(4,1) - t154 * t88 + t110 * t25 + t111 * t26 + (t381 - t190 * mrSges(4,3) + (t321 / 0.2e1 - t319 / 0.2e1) * t222 + (t387 - t282) * t223 + t363) * t222 + t292 * t191 + (-mrSges(7,1) * t293 + mrSges(7,2) * t294) * t81 + (-t11 * t294 + t12 * t293 + t188 * t2 - t189 * t3) * mrSges(7,3) + (t312 / 0.2e1 - t153 * mrSges(5,1) + t132 / 0.2e1 + t79 * mrSges(5,3) + qJD(4) * t165 + qJ(4) * t170) * t254 + (t313 / 0.2e1 + t153 * mrSges(5,2) + t133 / 0.2e1 - t78 * mrSges(5,3) - qJD(4) * t166 - qJ(4) * t171) * t252 + (t325 / 0.2e1 + t191 * mrSges(4,3) - t262) * t223 + (Ifges(4,5) + (Ifges(5,1) * t252 + t323) * t336 + (Ifges(5,2) * t254 + t324) * t337) * t213 + (Ifges(6,5) * t169 + Ifges(6,6) * t168) * t340 + (Ifges(7,5) * t125 + Ifges(7,6) * t126) * t342 + (Ifges(7,5) * t108 + Ifges(7,6) * t107) * t343 + (Ifges(6,1) * t169 + Ifges(6,4) * t168) * t347 + (Ifges(6,4) * t169 + Ifges(6,2) * t168) * t349; -t384 * t120 + t146 * t121 - t275 * t165 + t199 * t166 - t393 * t63 + t86 * t64 + t159 + t51 + t9 + (t11 * t86 - t12 * t393 + t56) * m(7) + (t146 * t42 - t384 * t43 + t122) * m(6) + (t112 * t199 - t113 * t275 + t153) * m(5); -m(7) * (t11 * t13 + t12 * t14) + (-Ifges(6,2) * t146 + t143 + t71) * t349 + t283 + t86 * t382 + t365 - t135 * (mrSges(6,1) * t146 + mrSges(6,2) * t384) - t42 * t120 + t43 * t121 - t14 * t63 - t13 * t64 + (t146 * t43 + t384 * t42) * mrSges(6,3) + (-t146 * t45 + t259 * t25 + t256 * t26 + (-t256 * t64 + t259 * t63) * qJD(6) + (-t146 * t81 + t2 * t256 + t259 * t3 + (-t11 * t256 + t12 * t259) * qJD(6)) * m(7)) * pkin(5) + (Ifges(6,5) * t384 - Ifges(6,6) * t146) * t340 + t70 * t346 + (Ifges(6,1) * t384 - t322) * t347 + t388; -t11 * t63 + t12 * t64 + t37 * t353 + t388;];
tauc  = t1(:);
