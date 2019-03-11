% Calculate time derivative of joint inertia matrix for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR12_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR12_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:35:19
% EndTime: 2019-03-09 14:35:34
% DurationCPUTime: 6.64s
% Computational Cost: add. (10263->605), mult. (24920->868), div. (0->0), fcn. (23716->10), ass. (0->288)
t228 = sin(qJ(6));
t232 = cos(qJ(6));
t191 = -mrSges(7,1) * t232 + mrSges(7,2) * t228;
t405 = mrSges(6,1) - t191;
t229 = sin(qJ(5));
t230 = sin(qJ(4));
t363 = cos(qJ(5));
t364 = cos(qJ(4));
t267 = t363 * t364;
t316 = qJD(5) * t229;
t317 = qJD(4) * t230;
t141 = -t229 * t317 - t230 * t316 + (qJD(4) + qJD(5)) * t267;
t224 = t228 ^ 2;
t225 = t232 ^ 2;
t399 = t224 + t225;
t273 = t399 * t141;
t181 = t229 * t230 - t267;
t332 = t181 * t228;
t283 = t332 / 0.2e1;
t277 = qJD(5) * t363;
t278 = qJD(4) * t364;
t292 = t363 * t230;
t296 = t229 * t364;
t140 = -qJD(4) * t292 - qJD(5) * t296 - t229 * t278 - t230 * t277;
t322 = t232 * t140;
t404 = qJD(6) * t283 + t322 / 0.2e1;
t314 = qJD(6) * t232;
t275 = t314 / 0.2e1;
t324 = t228 * t140;
t403 = t181 * t275 - t324 / 0.2e1;
t315 = qJD(6) * t228;
t247 = t181 * t315 + t322;
t366 = t228 / 0.2e1;
t365 = t232 / 0.2e1;
t276 = -t315 / 0.2e1;
t226 = sin(pkin(6));
t231 = sin(qJ(2));
t326 = t226 * t231;
t227 = cos(pkin(6));
t233 = cos(qJ(2));
t325 = t226 * t233;
t244 = -t227 * t364 + t230 * t325;
t211 = pkin(8) * t326;
t362 = pkin(1) * t233;
t298 = -pkin(2) - t362;
t136 = pkin(3) * t326 + t211 + (-pkin(9) + t298) * t227;
t234 = -pkin(2) - pkin(9);
t274 = -qJ(3) * t231 - pkin(1);
t156 = (t233 * t234 + t274) * t226;
t86 = t364 * t136 - t156 * t230;
t69 = pkin(4) * t326 + pkin(10) * t244 + t86;
t171 = -t227 * t230 - t325 * t364;
t87 = t230 * t136 + t364 * t156;
t75 = pkin(10) * t171 + t87;
t355 = t229 * t69 + t363 * t75;
t29 = pkin(11) * t326 + t355;
t214 = t227 * t231 * pkin(1);
t176 = pkin(8) * t325 + t214;
t164 = -t227 * qJ(3) - t176;
t155 = pkin(3) * t325 - t164;
t105 = -pkin(4) * t171 + t155;
t116 = t229 * t171 - t244 * t363;
t249 = t171 * t363 + t229 * t244;
t54 = -pkin(5) * t249 - pkin(11) * t116 + t105;
t16 = -t228 * t29 + t232 * t54;
t17 = t228 * t54 + t232 * t29;
t291 = qJD(2) * t326;
t149 = qJD(4) * t171 + t230 * t291;
t150 = qJD(4) * t244 + t291 * t364;
t62 = qJD(5) * t249 + t149 * t363 + t229 * t150;
t63 = qJD(5) * t116 + t229 * t149 - t150 * t363;
t310 = t227 * t362;
t203 = qJD(2) * t310;
t221 = t227 * qJD(3);
t370 = pkin(3) + pkin(8);
t135 = -t291 * t370 + t203 + t221;
t94 = -pkin(4) * t150 + t135;
t20 = pkin(5) * t63 - pkin(11) * t62 + t94;
t319 = qJD(2) * t233;
t290 = t226 * t319;
t202 = pkin(2) * t291;
t318 = qJD(3) * t231;
t123 = t202 + (-t318 + (pkin(9) * t231 - qJ(3) * t233) * qJD(2)) * t226;
t158 = (t325 * t370 + t214) * qJD(2);
t53 = -qJD(4) * t87 - t230 * t123 + t364 * t158;
t37 = pkin(4) * t290 - t149 * pkin(10) + t53;
t52 = t364 * t123 + t136 * t278 - t156 * t317 + t230 * t158;
t41 = pkin(10) * t150 + t52;
t8 = t229 * t37 + t69 * t277 - t316 * t75 + t363 * t41;
t5 = pkin(11) * t290 + t8;
t3 = -qJD(6) * t17 + t20 * t232 - t228 * t5;
t358 = t228 * t3;
t241 = -t358 + (-t16 * t232 - t17 * t228) * qJD(6);
t2 = qJD(6) * t16 + t20 * t228 + t232 * t5;
t359 = t2 * t232;
t402 = m(7) * (t241 + t359);
t401 = mrSges(4,2) - mrSges(3,1);
t104 = mrSges(6,1) * t326 - mrSges(6,3) * t116;
t98 = -t116 * t228 + t232 * t326;
t99 = t116 * t232 + t228 * t326;
t58 = -mrSges(7,1) * t98 + mrSges(7,2) * t99;
t400 = -t104 + t58;
t182 = t296 + t292;
t129 = -mrSges(7,2) * t182 + mrSges(7,3) * t332;
t330 = t181 * t232;
t130 = mrSges(7,1) * t182 + mrSges(7,3) * t330;
t398 = -t129 * t315 - t130 * t314;
t336 = t141 * t228;
t397 = -t182 * t314 - t336;
t334 = t141 * t232;
t396 = t182 * t315 - t334;
t65 = mrSges(7,2) * t249 + mrSges(7,3) * t98;
t66 = -mrSges(7,1) * t249 - mrSges(7,3) * t99;
t395 = -t66 * t314 - t65 * t315;
t394 = t230 * t52 + t364 * t53;
t393 = -t363 * t140 - t141 * t229;
t38 = qJD(6) * t98 + t228 * t290 + t232 * t62;
t18 = mrSges(7,1) * t63 - mrSges(7,3) * t38;
t39 = -qJD(6) * t99 - t228 * t62 + t232 * t290;
t19 = -mrSges(7,2) * t63 + mrSges(7,3) * t39;
t392 = -t228 * t18 + t232 * t19;
t32 = -t229 * t75 + t363 * t69;
t391 = -t140 * t32 - t355 * t141;
t15 = -mrSges(7,1) * t39 + mrSges(7,2) * t38;
t9 = -qJD(5) * t355 - t229 * t41 + t363 * t37;
t6 = -pkin(5) * t290 - t9;
t390 = -m(7) * t6 - t15;
t389 = Ifges(5,5) * t149 + Ifges(5,6) * t150 + Ifges(5,3) * t290;
t388 = -t16 * t314 - t17 * t315 - t358;
t204 = pkin(4) * t278 + qJD(3);
t78 = pkin(5) * t141 - pkin(11) * t140 + t204;
t215 = t230 * pkin(4) + qJ(3);
t124 = pkin(5) * t182 + pkin(11) * t181 + t215;
t356 = pkin(10) - t234;
t190 = t356 * t230;
t294 = t364 * t234;
t245 = -pkin(10) * t364 + t294;
t152 = -t190 * t363 + t229 * t245;
t82 = t124 * t228 + t152 * t232;
t151 = -t229 * t190 - t245 * t363;
t179 = t245 * qJD(4);
t262 = t356 * t317;
t88 = -qJD(5) * t151 + t179 * t363 + t229 * t262;
t26 = -qJD(6) * t82 - t228 * t88 + t232 * t78;
t342 = t228 * t26;
t81 = t124 * t232 - t152 * t228;
t387 = -t81 * t314 - t82 * t315 - t342;
t28 = -pkin(5) * t326 - t32;
t386 = -m(7) * t28 - t400;
t113 = -mrSges(5,2) * t290 + mrSges(5,3) * t150;
t160 = -mrSges(5,2) * t326 + mrSges(5,3) * t171;
t385 = m(5) * ((-t230 * t86 + t364 * t87) * qJD(4) + t394) + t160 * t278 + t230 * t113;
t56 = mrSges(6,1) * t290 - mrSges(6,3) * t62;
t384 = -m(6) * t9 - t390 - t56;
t383 = 0.2e1 * m(6);
t382 = 0.2e1 * m(7);
t381 = -0.2e1 * pkin(1);
t380 = 2 * mrSges(4,1);
t379 = -2 * mrSges(3,3);
t378 = -0.2e1 * mrSges(6,3);
t377 = 0.2e1 * t151;
t165 = (-pkin(2) * t233 + t274) * t226;
t376 = -0.2e1 * t165;
t375 = m(6) * pkin(4);
t374 = t38 / 0.2e1;
t373 = t39 / 0.2e1;
t372 = t98 / 0.2e1;
t371 = t99 / 0.2e1;
t219 = Ifges(7,5) * t314;
t368 = Ifges(7,6) * t276 + t219 / 0.2e1;
t367 = Ifges(7,5) * t366 + Ifges(7,6) * t365;
t361 = pkin(4) * t229;
t259 = mrSges(7,1) * t228 + mrSges(7,2) * t232;
t183 = t259 * qJD(6);
t360 = pkin(5) * t183;
t357 = Ifges(3,4) + Ifges(4,6);
t353 = Ifges(5,4) * t230;
t352 = Ifges(7,4) * t228;
t351 = Ifges(7,4) * t232;
t350 = Ifges(7,6) * t228;
t349 = pkin(4) * qJD(5);
t89 = qJD(5) * t152 + t229 * t179 - t262 * t363;
t347 = t151 * t89;
t168 = -pkin(8) * t291 + t203;
t346 = t168 * mrSges(3,2);
t345 = t181 * mrSges(6,3);
t344 = t182 * mrSges(6,3);
t25 = qJD(6) * t81 + t228 * t78 + t232 * t88;
t339 = t232 * t25;
t337 = t140 * t181;
t333 = t151 * t229;
t331 = t181 * t229;
t217 = pkin(11) + t361;
t329 = t217 * t228;
t328 = t217 * t232;
t308 = t363 * pkin(4);
t218 = -t308 - pkin(5);
t327 = t218 * t183;
t321 = Ifges(6,5) * t140 - Ifges(6,6) * t141;
t320 = Ifges(3,5) * t290 + Ifges(4,5) * t291;
t169 = t176 * qJD(2);
t313 = 0.2e1 * t169;
t312 = 0.2e1 * t233;
t311 = m(7) * t349;
t12 = Ifges(7,5) * t38 + Ifges(7,6) * t39 + Ifges(7,3) * t63;
t309 = Ifges(6,5) * t62 - Ifges(6,6) * t63 + Ifges(6,3) * t290;
t307 = pkin(4) * t316;
t306 = Ifges(5,4) * t364;
t297 = t228 * t363;
t295 = t232 * t363;
t289 = t234 * t317;
t282 = -t330 / 0.2e1;
t272 = mrSges(6,1) * t307;
t271 = t191 * t307;
t270 = pkin(4) * t277;
t269 = mrSges(4,1) * t290;
t268 = mrSges(5,2) * t278;
t265 = t290 / 0.2e1;
t261 = mrSges(7,3) * t270;
t260 = mrSges(6,2) * t270;
t258 = Ifges(7,1) * t232 - t352;
t257 = -Ifges(7,2) * t228 + t351;
t256 = t228 * t270;
t255 = t232 * t270;
t254 = -t140 * t151 + t181 * t89;
t253 = t224 * t261;
t252 = t225 * t261;
t248 = t181 * t314 - t324;
t186 = t257 * qJD(6);
t188 = t258 * qJD(6);
t194 = Ifges(7,2) * t232 + t352;
t196 = Ifges(7,1) * t228 + t351;
t246 = t232 * t186 + t228 * t188 - t194 * t315 + t196 * t314;
t243 = t399 * t363;
t242 = -t141 * mrSges(6,2) + mrSges(7,3) * t273 + t140 * t405 + t181 * t183;
t240 = -t342 + (-t228 * t82 - t232 * t81) * qJD(6);
t49 = Ifges(7,5) * t247 + Ifges(7,6) * t248 + Ifges(7,3) * t141;
t13 = Ifges(7,4) * t38 + Ifges(7,2) * t39 + Ifges(7,6) * t63;
t14 = Ifges(7,1) * t38 + Ifges(7,4) * t39 + Ifges(7,5) * t63;
t45 = Ifges(7,4) * t99 + Ifges(7,2) * t98 - Ifges(7,6) * t249;
t46 = Ifges(7,1) * t99 + Ifges(7,4) * t98 - Ifges(7,5) * t249;
t237 = t9 * mrSges(6,1) - t8 * mrSges(6,2) + mrSges(7,3) * t359 + t13 * t365 + t14 * t366 + t28 * t183 + t186 * t372 + t188 * t371 + t6 * t191 + t194 * t373 + t196 * t374 - t249 * t368 + t46 * t275 + t45 * t276 + t63 * t367 + t309;
t73 = mrSges(7,1) * t141 - mrSges(7,3) * t247;
t74 = -mrSges(7,2) * t141 + mrSges(7,3) * t248;
t236 = m(7) * (t339 + t387) + t232 * t74 - t228 * t73;
t101 = Ifges(7,6) * t182 - t181 * t257;
t102 = Ifges(7,5) * t182 - t181 * t258;
t50 = Ifges(7,4) * t247 + Ifges(7,2) * t248 + Ifges(7,6) * t141;
t51 = Ifges(7,1) * t247 + Ifges(7,4) * t248 + Ifges(7,5) * t141;
t235 = -t88 * mrSges(6,2) + mrSges(7,3) * t339 + t101 * t276 + t102 * t275 + t141 * t367 + t151 * t183 + t182 * t368 + t186 * t283 + t188 * t282 + t403 * t194 + t404 * t196 + t50 * t365 + t51 * t366 - t405 * t89 + t321;
t197 = Ifges(5,1) * t364 - t353;
t195 = -Ifges(5,2) * t230 + t306;
t192 = t230 * mrSges(5,1) + mrSges(5,2) * t364;
t189 = (-Ifges(5,1) * t230 - t306) * qJD(4);
t187 = (-Ifges(5,2) * t364 - t353) * qJD(4);
t184 = (mrSges(5,1) * t364 - mrSges(5,2) * t230) * qJD(4);
t180 = -mrSges(4,1) * t325 - mrSges(4,3) * t227;
t175 = -t211 + t310;
t166 = t227 * t298 + t211;
t163 = -t168 - t221;
t161 = mrSges(5,1) * t326 + mrSges(5,3) * t244;
t159 = t202 + (-qJ(3) * t319 - t318) * t226;
t148 = -Ifges(6,1) * t181 - Ifges(6,4) * t182;
t147 = -Ifges(6,4) * t181 - Ifges(6,2) * t182;
t146 = mrSges(6,1) * t182 - mrSges(6,2) * t181;
t120 = t259 * t181;
t118 = -mrSges(5,1) * t171 - mrSges(5,2) * t244;
t112 = mrSges(5,1) * t290 - mrSges(5,3) * t149;
t109 = -Ifges(5,1) * t244 + Ifges(5,4) * t171 + Ifges(5,5) * t326;
t108 = -Ifges(5,4) * t244 + Ifges(5,2) * t171 + Ifges(5,6) * t326;
t103 = -mrSges(6,2) * t326 + mrSges(6,3) * t249;
t100 = Ifges(7,3) * t182 + (-Ifges(7,5) * t232 + t350) * t181;
t93 = -mrSges(5,1) * t150 + mrSges(5,2) * t149;
t92 = Ifges(6,1) * t140 - Ifges(6,4) * t141;
t91 = Ifges(6,4) * t140 - Ifges(6,2) * t141;
t90 = mrSges(6,1) * t141 + mrSges(6,2) * t140;
t80 = Ifges(5,1) * t149 + Ifges(5,4) * t150 + Ifges(5,5) * t290;
t79 = Ifges(5,4) * t149 + Ifges(5,2) * t150 + Ifges(5,6) * t290;
t76 = -mrSges(6,1) * t249 + mrSges(6,2) * t116;
t71 = Ifges(6,1) * t116 + Ifges(6,4) * t249 + Ifges(6,5) * t326;
t70 = Ifges(6,4) * t116 + Ifges(6,2) * t249 + Ifges(6,6) * t326;
t64 = -mrSges(7,1) * t248 + mrSges(7,2) * t247;
t57 = -mrSges(6,2) * t290 - mrSges(6,3) * t63;
t44 = Ifges(7,5) * t99 + Ifges(7,6) * t98 - Ifges(7,3) * t249;
t24 = mrSges(6,1) * t63 + mrSges(6,2) * t62;
t22 = Ifges(6,1) * t62 - Ifges(6,4) * t63 + Ifges(6,5) * t290;
t21 = Ifges(6,4) * t62 - Ifges(6,2) * t63 + Ifges(6,6) * t290;
t1 = [0.2e1 * t163 * t180 + t171 * t79 + 0.2e1 * t155 * t93 + 0.2e1 * t52 * t160 + 0.2e1 * t53 * t161 + t149 * t109 + t150 * t108 + 0.2e1 * t135 * t118 + 0.2e1 * t87 * t113 + t116 * t22 + 0.2e1 * t8 * t103 + 0.2e1 * t9 * t104 + 0.2e1 * t105 * t24 + 0.2e1 * t86 * t112 + 0.2e1 * t94 * t76 + t98 * t13 + t99 * t14 + t62 * t71 + 0.2e1 * t2 * t65 + 0.2e1 * t3 * t66 + 0.2e1 * t6 * t58 + 0.2e1 * t32 * t56 + t39 * t45 + t38 * t46 + 0.2e1 * t28 * t15 + 0.2e1 * t17 * t19 + 0.2e1 * t16 * t18 - (t12 - t21) * t249 - t244 * t80 + ((t159 * mrSges(4,2) + t168 * mrSges(3,3)) * t312 + (-0.2e1 * t159 * mrSges(4,3) + (mrSges(4,1) + mrSges(3,3)) * t313 + t309 + t389) * t231 + ((t164 * t380 + mrSges(4,2) * t376 + t176 * t379 + (Ifges(4,5) - (2 * Ifges(3,6))) * t227 + (mrSges(3,1) * t381 - 0.2e1 * t231 * t357) * t226) * t231 + (t166 * t380 + t175 * t379 + mrSges(4,3) * t376 - Ifges(5,5) * t244 + Ifges(6,5) * t116 + Ifges(5,6) * t171 + Ifges(6,6) * t249 + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t227 + (mrSges(3,2) * t381 + t312 * t357) * t226 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + Ifges(5,3) + Ifges(6,3)) * t326) * t233) * qJD(2)) * t226 + 0.2e1 * m(4) * (t159 * t165 + t163 * t164 + t166 * t169) + 0.2e1 * m(3) * (t168 * t176 - t169 * t175) + 0.2e1 * m(5) * (t135 * t155 + t52 * t87 + t53 * t86) + (-t70 + t44) * t63 + 0.2e1 * t355 * t57 + (t105 * t94 + t32 * t9 + t355 * t8) * t383 + (t16 * t3 + t17 * t2 + t28 * t6) * t382 + (t401 * t313 + t320 - 0.2e1 * t346) * t227; t320 + m(7) * (t16 * t26 + t17 * t25 + t2 * t82 + t3 * t81) + t112 * t294 - t230 * t79 / 0.2e1 + t204 * t76 + t215 * t24 + t135 * t192 + t150 * t195 / 0.2e1 + t149 * t197 / 0.2e1 + t155 * t184 + t171 * t187 / 0.2e1 + t152 * t57 - t163 * mrSges(4,3) + t94 * t146 + t62 * t148 / 0.2e1 + t140 * t71 / 0.2e1 - t6 * t120 + t2 * t129 + t3 * t130 + t116 * t92 / 0.2e1 + t88 * t103 + t105 * t90 + t81 * t18 + t82 * t19 + t16 * t73 + t17 * t74 + t28 * t64 + t25 * t65 + t26 * t66 - (t49 / 0.2e1 - t91 / 0.2e1) * t249 - t244 * t189 / 0.2e1 + t14 * t282 + t13 * t283 + t403 * t45 + t404 * t46 - t346 + t364 * t80 / 0.2e1 + (Ifges(5,5) * t364 - Ifges(5,6) * t230) * t265 - t8 * t344 - t109 * t317 / 0.2e1 - Ifges(3,6) * t291 - Ifges(4,4) * t290 - t161 * t289 - t108 * t278 / 0.2e1 + (t12 / 0.2e1 - t21 / 0.2e1 - Ifges(6,6) * t265) * t182 + (t321 + (-Ifges(5,5) * t230 - Ifges(5,6) * t364) * qJD(4)) * t326 / 0.2e1 + m(6) * (t105 * t204 + t152 * t8 + t215 * t94 + t355 * t88) + (-t22 / 0.2e1 - Ifges(6,5) * t265) * t181 - pkin(2) * t269 + t384 * t151 + t385 * t234 + (-m(6) * t32 - t386) * t89 + t391 * mrSges(6,3) + (-t278 * t87 + t317 * t86 - t394) * mrSges(5,3) + t9 * t345 + t51 * t371 + t50 * t372 + t101 * t373 + t102 * t374 + (-m(4) * pkin(2) + t401) * t169 + (-m(4) * t163 + m(5) * t135 - mrSges(4,1) * t291 + t93) * qJ(3) + (-m(4) * t164 + m(5) * t155 + t118 - t180) * qJD(3) + (t44 / 0.2e1 - t70 / 0.2e1) * t141 + (-t147 / 0.2e1 + t100 / 0.2e1) * t63; t364 * t189 + 0.2e1 * qJ(3) * t184 - 0.2e1 * t89 * t120 + 0.2e1 * t25 * t129 + 0.2e1 * t26 * t130 + 0.2e1 * t204 * t146 + t64 * t377 - t230 * t187 + 0.2e1 * t215 * t90 + 0.2e1 * t81 * t73 + 0.2e1 * t82 * t74 + (-t195 * t364 - t230 * t197) * qJD(4) + (t152 * t88 + t204 * t215 + t347) * t383 + (t25 * t82 + t26 * t81 + t347) * t382 + (t378 * t88 + t49 - t91) * t182 + (t152 * t378 + t100 - t147) * t141 + (mrSges(6,3) * t377 - t101 * t228 + t102 * t232 + t148) * t140 + 0.2e1 * (mrSges(4,3) + t192 + (m(4) + m(5)) * qJ(3)) * qJD(3) + (t89 * t378 + t228 * t50 - t232 * t51 - t92 + (t101 * t232 + t102 * t228) * qJD(6)) * t181; -m(6) * t391 - t161 * t317 + t364 * t112 + m(4) * t169 + t269 + t397 * t66 - t396 * t65 + (m(7) * (-t16 * t228 + t17 * t232) + t103) * t141 + t384 * t181 + t386 * t140 + (m(6) * t8 + t392 + t402 + t57) * t182 + t385; t181 * t64 + (t129 * t232 - t130 * t228) * t141 + (t120 + 0.2e1 * t345) * t140 + m(7) * (t334 * t82 - t336 * t81 + t254) + m(6) * (t141 * t152 + t254) + (t141 * t378 + (-t129 * t228 - t130 * t232) * qJD(6) + m(6) * t88 + t236) * t182; 0.2e1 * m(6) * (t141 * t182 - t337) + 0.2e1 * m(7) * (t182 * t273 - t337); t65 * t255 - t66 * t256 + t19 * t328 + t57 * t361 - t18 * t329 + t56 * t308 + t218 * t15 - t52 * mrSges(5,2) + t53 * mrSges(5,1) + t103 * t270 + m(7) * (t218 * t6 + (-t16 * t297 + t17 * t295 + t229 * t28) * t349) + t237 + (t363 * t9 + t229 * t8 + (-t229 * t32 + t355 * t363) * qJD(5)) * t375 + t400 * t307 + (t395 + t402) * t217 + t388 * mrSges(7,3) + t389; t129 * t255 - t130 * t256 - t270 * t344 + t74 * t328 - Ifges(5,5) * t317 - t73 * t329 - mrSges(5,1) * t289 - t234 * t268 + t218 * t64 + m(7) * (t218 * t89 + (t295 * t82 - t297 * t81 + t333) * t349) + t235 + (-t363 * t89 + t229 * t88 + (t152 * t363 + t333) * qJD(5)) * t375 - Ifges(5,6) * t278 + t393 * mrSges(6,3) * pkin(4) + (-t345 - t120) * t307 + (m(7) * (t240 + t339) + t398) * t217 + t387 * mrSges(7,3); -t268 - mrSges(5,1) * t317 + ((t182 * t363 + t331) * qJD(5) - t393) * t375 + m(7) * (-t218 * t140 + t217 * t273 + (t182 * t243 + t331) * t349) + t242; 0.2e1 * t271 + 0.2e1 * t327 + 0.2e1 * (t217 * t243 + t218 * t229) * t311 - 0.2e1 * t260 - 0.2e1 * t272 + 0.2e1 * t253 + 0.2e1 * t252 + t246; t390 * pkin(5) + (m(7) * (t359 + t388) + t392 + t395) * pkin(11) + t237 + t241 * mrSges(7,3); (-m(7) * t89 - t64) * pkin(5) + (t236 + t398) * pkin(11) + t240 * mrSges(7,3) + t235; m(7) * (pkin(5) * t140 + pkin(11) * t273) + t242; t271 + t327 + (-pkin(5) * t229 + pkin(11) * t243) * t311 - t360 + t253 + t252 - t260 - t272 + t246; t246 - 0.2e1 * t360; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t12; mrSges(7,1) * t26 - mrSges(7,2) * t25 + t49; mrSges(7,1) * t397 + mrSges(7,2) * t396; t219 + (-mrSges(7,1) * t297 - mrSges(7,2) * t295) * t349 + (t191 * t217 - t350) * qJD(6); t219 + (pkin(11) * t191 - t350) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
