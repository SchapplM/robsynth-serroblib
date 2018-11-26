% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:31
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:30:33
% EndTime: 2018-11-23 16:30:56
% DurationCPUTime: 22.54s
% Computational Cost: add. (24327->743), mult. (80611->1033), div. (0->0), fcn. (68496->12), ass. (0->311)
t452 = Ifges(7,4) + Ifges(6,4);
t267 = sin(pkin(12));
t269 = sin(pkin(6));
t270 = cos(pkin(12));
t278 = cos(qJ(3));
t271 = cos(pkin(7));
t275 = sin(qJ(3));
t356 = t271 * t275;
t286 = (-t267 * t356 + t270 * t278) * t269;
t228 = qJD(1) * t286;
t268 = sin(pkin(7));
t340 = qJD(3) * t278;
t473 = t268 * t340 - t228;
t272 = cos(pkin(6));
t360 = t268 * t275;
t217 = t272 * t360 + (t267 * t278 + t270 * t356) * t269;
t210 = t217 * qJD(1);
t358 = t269 * t270;
t239 = -t268 * t358 + t271 * t272;
t231 = qJD(1) * t239 + qJD(3);
t274 = sin(qJ(4));
t277 = cos(qJ(4));
t186 = t210 * t277 + t231 * t274;
t355 = t271 * t278;
t359 = t268 * t278;
t282 = t272 * t359 + t269 * (-t267 * t275 + t270 * t355);
t211 = t282 * qJD(3);
t200 = qJD(1) * t211;
t153 = qJD(4) * t186 + t200 * t274;
t416 = t153 / 0.2e1;
t185 = -t210 * t274 + t231 * t277;
t152 = qJD(4) * t185 + t200 * t277;
t209 = t282 * qJD(1);
t207 = qJD(4) - t209;
t273 = sin(qJ(5));
t276 = cos(qJ(5));
t155 = t186 * t276 + t207 * t273;
t212 = t217 * qJD(3);
t201 = qJD(1) * t212;
t83 = -qJD(5) * t155 - t152 * t273 + t201 * t276;
t424 = t83 / 0.2e1;
t154 = -t186 * t273 + t207 * t276;
t82 = qJD(5) * t154 + t152 * t276 + t201 * t273;
t425 = t82 / 0.2e1;
t451 = Ifges(7,5) + Ifges(6,5);
t453 = Ifges(7,1) + Ifges(6,1);
t464 = t451 * t416 + t452 * t424 + t453 * t425;
t450 = Ifges(7,2) + Ifges(6,2);
t449 = Ifges(7,6) + Ifges(6,6);
t448 = Ifges(7,3) + Ifges(6,3);
t184 = qJD(5) - t185;
t469 = t452 * t155;
t444 = t154 * t450 + t184 * t449 + t469;
t472 = -t444 / 0.2e1;
t471 = t452 * t154;
t242 = -t277 * t271 + t274 * t360;
t343 = qJD(1) * t269;
t326 = t267 * t343;
t318 = t268 * t326;
t436 = -qJD(4) * t242 - t274 * t318 + t277 * t473;
t287 = (t267 * t355 + t270 * t275) * t269;
t227 = qJD(1) * t287;
t341 = qJD(3) * t275;
t470 = t268 * t341 - t227;
t468 = t452 * t276;
t467 = t452 * t273;
t261 = qJ(2) * t358;
t394 = pkin(1) * t272;
t333 = qJD(1) * t394;
t236 = qJD(1) * t261 + t267 * t333;
t288 = (t268 * t272 + t271 * t358) * pkin(9);
t202 = qJD(1) * t288 + t236;
t259 = t270 * t333;
t362 = t267 * t269;
t285 = pkin(2) * t272 + (-pkin(9) * t271 - qJ(2)) * t362;
t208 = qJD(1) * t285 + t259;
t229 = (-pkin(9) * t267 * t268 - pkin(2) * t270 - pkin(1)) * t269;
t222 = qJD(1) * t229 + qJD(2);
t296 = t208 * t271 + t222 * t268;
t144 = -t275 * t202 + t296 * t278;
t123 = -t231 * pkin(3) - t144;
t183 = Ifges(5,4) * t185;
t373 = t207 * Ifges(5,5);
t181 = -t208 * t268 + t271 * t222;
t121 = -pkin(3) * t209 - pkin(10) * t210 + t181;
t145 = t202 * t278 + t275 * t296;
t124 = t231 * pkin(10) + t145;
t64 = t121 * t274 + t124 * t277;
t57 = pkin(11) * t207 + t64;
t76 = -t185 * pkin(4) - t186 * pkin(11) + t123;
t22 = -t273 * t57 + t276 * t76;
t12 = -qJ(6) * t155 + t22;
t10 = pkin(5) * t184 + t12;
t23 = t273 * t76 + t276 * t57;
t13 = qJ(6) * t154 + t23;
t297 = t22 * t276 + t23 * t273;
t311 = mrSges(7,1) * t273 + mrSges(7,2) * t276;
t313 = mrSges(6,1) * t273 + mrSges(6,2) * t276;
t396 = t276 / 0.2e1;
t399 = -t273 / 0.2e1;
t405 = t184 / 0.2e1;
t412 = t155 / 0.2e1;
t414 = t154 / 0.2e1;
t63 = t121 * t277 - t274 * t124;
t56 = -pkin(4) * t207 - t63;
t44 = -pkin(5) * t154 + qJD(6) + t56;
t443 = t155 * t453 + t451 * t184 + t471;
t457 = t276 * t453 - t467;
t458 = -t273 * t450 + t468;
t459 = -t273 * t449 + t276 * t451;
t427 = t297 * mrSges(6,3) + (t10 * t276 + t13 * t273) * mrSges(7,3) - t311 * t44 - t313 * t56 - t458 * t414 - t457 * t412 - t459 * t405 - t444 * t399 - t443 * t396;
t466 = t427 - t123 * mrSges(5,2) - t183 / 0.2e1 - t373 / 0.2e1 + t63 * mrSges(5,3);
t372 = t207 * Ifges(5,6);
t375 = t185 * Ifges(5,2);
t385 = Ifges(5,4) * t186;
t110 = t372 + t375 + t385;
t465 = -t22 * mrSges(6,1) - t10 * mrSges(7,1) + t23 * mrSges(6,2) + t13 * mrSges(7,2) + t110 / 0.2e1;
t447 = t153 * t448 + t449 * t83 + t451 * t82;
t446 = t153 * t449 + t450 * t83 + t452 * t82;
t68 = t155 * Ifges(7,5) + t154 * Ifges(7,6) + t184 * Ifges(7,3);
t69 = t155 * Ifges(6,5) + t154 * Ifges(6,6) + t184 * Ifges(6,3);
t463 = t69 + t68;
t351 = t276 * t277;
t172 = t209 * t351 + t210 * t273;
t253 = -pkin(4) * t277 - pkin(11) * t274 - pkin(3);
t264 = pkin(10) * t351;
t335 = qJD(6) * t276;
t315 = pkin(4) * t274 - pkin(11) * t277;
t251 = t315 * qJD(4);
t339 = qJD(4) * t274;
t392 = pkin(10) * t273;
t345 = t276 * t251 + t339 * t392;
t363 = t209 * t274;
t104 = t209 * t315 + t145;
t176 = pkin(3) * t210 - pkin(10) * t209;
t103 = t277 * t144 + t274 * t176;
t91 = pkin(11) * t210 + t103;
t42 = t276 * t104 - t273 * t91;
t462 = -pkin(5) * t363 + qJ(6) * t172 - t274 * t335 + (pkin(5) * t274 - qJ(6) * t351) * qJD(4) + (-t264 + (qJ(6) * t274 - t253) * t273) * qJD(5) + t345 - t42;
t353 = t273 * t277;
t171 = -t209 * t353 + t210 * t276;
t336 = qJD(5) * t276;
t346 = t273 * t251 + t253 * t336;
t352 = t274 * t276;
t43 = t273 * t104 + t276 * t91;
t461 = -qJ(6) * t171 + (-pkin(10) * qJD(4) - qJ(6) * qJD(5)) * t352 + (-qJD(6) * t274 + (-pkin(10) * qJD(5) - qJ(6) * qJD(4)) * t277) * t273 + t346 - t43;
t338 = qJD(4) * t277;
t102 = -t274 * t144 + t176 * t277;
t90 = -pkin(4) * t210 - t102;
t460 = pkin(10) * t338 - t90 + (t273 * t338 + t274 * t336 + t171) * pkin(5);
t344 = t267 * t394 + t261;
t213 = t288 + t344;
t263 = t270 * t394;
t218 = t263 + t285;
t295 = t218 * t271 + t229 * t268;
t158 = -t275 * t213 + t295 * t278;
t403 = t201 / 0.2e1;
t417 = -t153 / 0.2e1;
t418 = t152 / 0.2e1;
t423 = Ifges(5,1) * t418 + Ifges(5,4) * t417 + Ifges(5,5) * t403;
t283 = qJD(2) * t286;
t116 = qJD(1) * t283 + qJD(3) * t144;
t342 = qJD(2) * t269;
t325 = t267 * t342;
t316 = qJD(1) * t325;
t293 = t268 * t316;
t165 = pkin(3) * t201 - pkin(10) * t200 + t293;
t29 = t277 * t116 + t121 * t338 - t124 * t339 + t274 * t165;
t25 = pkin(11) * t201 + t29;
t284 = qJD(2) * t287;
t117 = qJD(1) * t284 + qJD(3) * t145;
t54 = t153 * pkin(4) - t152 * pkin(11) + t117;
t6 = -qJD(5) * t23 - t25 * t273 + t276 * t54;
t1 = pkin(5) * t153 - qJ(6) * t82 - qJD(6) * t155 + t6;
t337 = qJD(5) * t273;
t5 = t276 * t25 + t273 * t54 + t76 * t336 - t337 * t57;
t2 = qJ(6) * t83 + qJD(6) * t154 + t5;
t456 = t6 * mrSges(6,1) + t1 * mrSges(7,1) - t5 * mrSges(6,2) - t2 * mrSges(7,2);
t455 = t63 * mrSges(5,1);
t454 = t64 * mrSges(5,2);
t442 = Ifges(4,5) * t200;
t441 = Ifges(4,6) * t201;
t440 = t116 * mrSges(4,2);
t386 = mrSges(4,3) * t210;
t350 = -mrSges(4,1) * t231 - mrSges(5,1) * t185 + mrSges(5,2) * t186 + t386;
t243 = t271 * t274 + t277 * t360;
t292 = -t276 * t243 + t273 * t359;
t438 = qJD(5) * t292 - t273 * t436 + t276 * t470;
t223 = -t273 * t243 - t276 * t359;
t437 = qJD(5) * t223 + t273 * t470 + t276 * t436;
t435 = qJD(4) * t243 + t274 * t473 + t277 * t318;
t233 = t273 * t253 + t264;
t434 = t273 * t451 + t276 * t449;
t433 = t276 * t450 + t467;
t432 = t273 * t453 + t468;
t30 = -t274 * t116 - t121 * t339 - t124 * t338 + t165 * t277;
t431 = -t274 * t30 + t277 * t29;
t430 = -t273 * t6 + t276 * t5;
t429 = t30 * mrSges(5,1) - t29 * mrSges(5,2);
t426 = (-Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t154 + (-Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1) * t155 + (-Ifges(7,3) / 0.2e1 - Ifges(6,3) / 0.2e1) * t184 + t64 * mrSges(5,3) - t68 / 0.2e1 - t69 / 0.2e1 + t385 / 0.2e1 - t123 * mrSges(5,1) + t372 / 0.2e1 + t465;
t374 = t186 * Ifges(5,1);
t111 = t183 + t373 + t374;
t419 = t111 / 0.2e1;
t415 = -t154 / 0.2e1;
t413 = -t155 / 0.2e1;
t406 = -t184 / 0.2e1;
t400 = t210 / 0.2e1;
t393 = pkin(5) * t273;
t388 = -qJ(6) - pkin(11);
t187 = -t218 * t268 + t271 * t229;
t136 = -pkin(3) * t282 - pkin(10) * t217 + t187;
t205 = t278 * t213;
t159 = t218 * t356 + t229 * t360 + t205;
t143 = pkin(10) * t239 + t159;
t87 = t274 * t136 + t277 * t143;
t67 = -pkin(11) * t282 + t87;
t142 = -t239 * pkin(3) - t158;
t190 = t217 * t274 - t277 * t239;
t191 = t217 * t277 + t239 * t274;
t94 = t190 * pkin(4) - t191 * pkin(11) + t142;
t34 = t273 * t94 + t276 * t67;
t387 = mrSges(4,3) * t209;
t384 = Ifges(5,4) * t274;
t383 = Ifges(5,4) * t277;
t371 = t210 * Ifges(4,4);
t26 = -pkin(4) * t201 - t30;
t370 = t26 * t274;
t114 = mrSges(5,1) * t201 - mrSges(5,3) * t152;
t32 = -mrSges(6,1) * t83 + mrSges(6,2) * t82;
t367 = -t114 + t32;
t134 = pkin(4) * t186 - pkin(11) * t185;
t46 = t273 * t134 + t276 * t63;
t157 = mrSges(5,1) * t207 - mrSges(5,3) * t186;
t99 = -mrSges(6,1) * t154 + mrSges(6,2) * t155;
t366 = t99 - t157;
t365 = t117 * t278;
t364 = t185 * t273;
t357 = t270 * (-mrSges(3,2) * t272 + mrSges(3,3) * t358) * qJD(1);
t354 = t273 * t274;
t347 = -t441 + t442;
t327 = Ifges(5,5) * t152 - Ifges(5,6) * t153 + Ifges(5,3) * t201;
t31 = -t83 * mrSges(7,1) + t82 * mrSges(7,2);
t33 = -t273 * t67 + t276 * t94;
t319 = qJD(5) * t388;
t45 = t276 * t134 - t273 * t63;
t86 = t136 * t277 - t274 * t143;
t317 = t268 * t325;
t314 = mrSges(6,1) * t276 - mrSges(6,2) * t273;
t312 = mrSges(7,1) * t276 - mrSges(7,2) * t273;
t167 = t191 * t276 - t273 * t282;
t166 = -t191 * t273 - t276 * t282;
t294 = -(-qJ(2) * t326 + t259) * t267 + t236 * t270;
t129 = qJD(3) * t158 + t283;
t169 = pkin(3) * t212 - pkin(10) * t211 + t317;
t41 = -t274 * t129 - t136 * t339 - t143 * t338 + t169 * t277;
t66 = pkin(4) * t282 - t86;
t40 = t277 * t129 + t136 * t338 - t143 * t339 + t274 * t169;
t37 = pkin(11) * t212 + t40;
t130 = t284 + (t275 * t295 + t205) * qJD(3);
t163 = -qJD(4) * t190 + t211 * t277;
t164 = qJD(4) * t191 + t211 * t274;
t60 = t164 * pkin(4) - t163 * pkin(11) + t130;
t7 = t273 * t60 + t276 * t37 + t94 * t336 - t337 * t67;
t240 = (mrSges(3,1) * t272 - mrSges(3,3) * t362) * qJD(1);
t38 = -pkin(4) * t212 - t41;
t8 = -qJD(5) * t34 - t273 * t37 + t276 * t60;
t266 = -pkin(5) * t276 - pkin(4);
t255 = t388 * t276;
t254 = t388 * t273;
t252 = (pkin(10) + t393) * t274;
t249 = t276 * t253;
t238 = -qJD(6) * t273 + t276 * t319;
t237 = t273 * t319 + t335;
t232 = -pkin(10) * t353 + t249;
t225 = -qJ(6) * t354 + t233;
t215 = -qJ(6) * t352 + t249 + (-pkin(5) - t392) * t277;
t206 = Ifges(4,4) * t209;
t193 = -qJD(5) * t233 + t345;
t192 = (-t276 * t339 - t277 * t337) * pkin(10) + t346;
t188 = -mrSges(4,2) * t231 + t387;
t175 = -mrSges(4,1) * t209 + mrSges(4,2) * t210;
t170 = mrSges(4,1) * t201 + mrSges(4,2) * t200;
t162 = t210 * Ifges(4,1) + t231 * Ifges(4,5) + t206;
t161 = t209 * Ifges(4,2) + t231 * Ifges(4,6) + t371;
t156 = -mrSges(5,2) * t207 + mrSges(5,3) * t185;
t115 = -mrSges(5,2) * t201 - mrSges(5,3) * t153;
t109 = t186 * Ifges(5,5) + t185 * Ifges(5,6) + t207 * Ifges(5,3);
t108 = mrSges(6,1) * t184 - mrSges(6,3) * t155;
t107 = mrSges(7,1) * t184 - mrSges(7,3) * t155;
t106 = -mrSges(6,2) * t184 + mrSges(6,3) * t154;
t105 = -mrSges(7,2) * t184 + mrSges(7,3) * t154;
t98 = -mrSges(7,1) * t154 + mrSges(7,2) * t155;
t97 = mrSges(5,1) * t153 + mrSges(5,2) * t152;
t96 = -qJD(5) * t167 - t163 * t273 + t212 * t276;
t95 = qJD(5) * t166 + t163 * t276 + t212 * t273;
t84 = t152 * Ifges(5,4) - t153 * Ifges(5,2) + t201 * Ifges(5,6);
t55 = pkin(5) * t364 + t64;
t51 = -mrSges(6,2) * t153 + mrSges(6,3) * t83;
t50 = -mrSges(7,2) * t153 + mrSges(7,3) * t83;
t49 = mrSges(6,1) * t153 - mrSges(6,3) * t82;
t48 = mrSges(7,1) * t153 - mrSges(7,3) * t82;
t47 = -pkin(5) * t166 + t66;
t39 = -qJ(6) * t364 + t46;
t27 = -qJ(6) * t185 * t276 + pkin(5) * t186 + t45;
t21 = qJ(6) * t166 + t34;
t14 = pkin(5) * t190 - qJ(6) * t167 + t33;
t11 = -pkin(5) * t96 + t38;
t9 = -pkin(5) * t83 + t26;
t4 = qJ(6) * t96 + qJD(6) * t166 + t7;
t3 = pkin(5) * t164 - qJ(6) * t95 - qJD(6) * t167 + t8;
t15 = [(-Ifges(5,2) * t417 + t451 * t425 + t447 / 0.2e1 + t448 * t416 + t449 * t424 - t29 * mrSges(5,3) - Ifges(5,6) * t403 - Ifges(5,4) * t418 + mrSges(5,1) * t117 - t84 / 0.2e1 + t456) * t190 + t463 * t164 / 0.2e1 - (mrSges(4,1) * t293 - t116 * mrSges(4,3) - Ifges(4,4) * t200 + Ifges(4,2) * t201 + t327 / 0.2e1 + Ifges(5,3) * t403 + Ifges(5,6) * t417 + Ifges(5,5) * t418 + t429) * t282 + m(4) * (t116 * t159 - t117 * t158 + t129 * t145 + (qJD(1) * t187 + t181) * t317) + (mrSges(4,2) * t293 + mrSges(4,3) * t117 + Ifges(4,1) * t200 - Ifges(4,4) * t201) * t217 + m(5) * (t117 * t142 + t29 * t87 + t30 * t86 + t40 * t64 + t41 * t63) + (0.2e1 * t357 + m(3) * ((t270 * t344 + (qJ(2) * t362 - t263) * t267) * qJD(1) + t294)) * t342 + (t164 * t451 + t452 * t96 + t453 * t95) * t412 - t212 * t454 + (mrSges(5,2) * t117 - mrSges(5,3) * t30 + 0.2e1 * t423) * t191 + (-m(4) * t144 + m(5) * t123 + t350) * t130 + (-t144 * t211 - t145 * t212 - t158 * t200 - t159 * t201) * mrSges(4,3) + t231 * (Ifges(4,5) * t211 - Ifges(4,6) * t212) / 0.2e1 + t175 * t317 + (-t440 + t442 / 0.2e1 - t441 / 0.2e1 - t117 * mrSges(4,1) + t347 / 0.2e1) * t239 + t443 * t95 / 0.2e1 + t444 * t96 / 0.2e1 + (t164 * t448 + t449 * t96 + t451 * t95) * t405 + (t164 * t449 + t450 * t96 + t452 * t95) * t414 + (-t163 * t63 - t164 * t64) * mrSges(5,3) + m(7) * (t1 * t14 + t10 * t3 + t11 * t44 + t13 * t4 + t2 * t21 + t47 * t9) + m(6) * (t22 * t8 + t23 * t7 + t26 * t66 + t33 * t6 + t34 * t5 + t38 * t56) + (t452 * t425 + t446 / 0.2e1 + t449 * t416 + t450 * t424 + t2 * mrSges(7,3) + t5 * mrSges(6,3) - t26 * mrSges(6,1) - t9 * mrSges(7,1)) * t166 - t212 * t161 / 0.2e1 + t211 * t162 / 0.2e1 + t207 * (Ifges(5,5) * t163 - Ifges(5,6) * t164 + Ifges(5,3) * t212) / 0.2e1 + t186 * (Ifges(5,1) * t163 - Ifges(5,4) * t164 + Ifges(5,5) * t212) / 0.2e1 + t185 * (Ifges(5,4) * t163 - Ifges(5,2) * t164 + Ifges(5,6) * t212) / 0.2e1 + t212 * t109 / 0.2e1 + t209 * (Ifges(4,4) * t211 - Ifges(4,2) * t212) / 0.2e1 + t181 * (mrSges(4,1) * t212 + mrSges(4,2) * t211) + t187 * t170 + t129 * t188 + t23 * (-mrSges(6,2) * t164 + mrSges(6,3) * t96) + t13 * (-mrSges(7,2) * t164 + mrSges(7,3) * t96) + t22 * (mrSges(6,1) * t164 - mrSges(6,3) * t95) + t10 * (mrSges(7,1) * t164 - mrSges(7,3) * t95) + t123 * (mrSges(5,1) * t164 + mrSges(5,2) * t163) - t164 * t110 / 0.2e1 + t40 * t156 + t41 * t157 + t142 * t97 + t86 * t114 + t87 * t115 + t8 * t108 + t4 * t105 + t7 * t106 + t3 * t107 + t56 * (-mrSges(6,1) * t96 + mrSges(6,2) * t95) + t44 * (-mrSges(7,1) * t96 + mrSges(7,2) * t95) + t11 * t98 + t38 * t99 + t66 * t32 + t14 * t48 + t33 * t49 + t21 * t50 + t34 * t51 + t47 * t31 - 0.2e1 * t240 * t325 + (t26 * mrSges(6,2) + t9 * mrSges(7,2) - t6 * mrSges(6,3) - t1 * mrSges(7,3) + 0.2e1 * t464) * t167 + (Ifges(4,1) * t211 - Ifges(4,4) * t212) * t400 + t163 * t419 + t212 * t455; -m(4) * (-t144 * t227 + t145 * t228) + t243 * t115 + t436 * t156 - (t50 + t51) * t292 + (t48 + t49) * t223 - t228 * t188 + (t31 + t367) * t242 - t350 * t227 + t271 * t170 + (t108 + t107) * t438 + (t105 + t106) * t437 + (-m(3) * t294 + t267 * t240 - t357) * t343 + t435 * (t98 + t366) + (t1 * t223 + t10 * t438 + t13 * t437 - t2 * t292 + t242 * t9 + t435 * t44) * m(7) + (t22 * t438 + t223 * t6 + t23 * t437 + t242 * t26 - t292 * t5 + t435 * t56) * m(6) + (-t123 * t227 - t30 * t242 + t29 * t243 - t435 * t63 + t436 * t64) * m(5) + (m(5) * (t123 * t341 - t365) - t175 * t326 - t278 * t97 + (-t200 * t278 - t201 * t275) * mrSges(4,3) + (t188 * t278 + t275 * t350) * qJD(3) + (t116 * t275 - t144 * t341 + t145 * t340 - t181 * t326 + t271 * t316 - t365) * m(4)) * t268; -t446 * t354 / 0.2e1 - (t111 * t209 + t447) * t277 / 0.2e1 + (t274 * t459 - t277 * t448) * t416 + (t274 * t458 - t277 * t449) * t424 + (t56 * mrSges(6,1) + t44 * mrSges(7,1) - t23 * mrSges(6,3) - t13 * mrSges(7,3) + t406 * t449 + t413 * t452 + t415 * t450 + t472) * t171 + (-t443 / 0.2e1 - t56 * mrSges(6,2) - t44 * mrSges(7,2) + t22 * mrSges(6,3) + t10 * mrSges(7,3) + t452 * t415 + t453 * t413 + t451 * t406) * t172 + (-t463 / 0.2e1 + t449 * t415 + t451 * t413 + t448 * t406 + t465) * t363 + (((-m(5) * t64 - t156) * pkin(10) - t375 / 0.2e1 - t426) * t274 + (t419 + (-m(5) * t63 + m(6) * t56 + t366) * pkin(10) + t374 / 0.2e1 - t466) * t277) * qJD(4) - t210 * t455 + t462 * t107 + (t1 * t215 + t10 * t462 + t13 * t461 + t2 * t225 + t252 * t9 + t44 * t460) * m(7) + (t274 * t457 - t277 * t451) * t425 - t123 * (mrSges(5,1) * t274 + mrSges(5,2) * t277) * t209 + (t193 - t42) * t108 - m(6) * (t22 * t42 + t23 * t43 + t56 * t90) + t252 * t31 - t231 * (Ifges(4,5) * t209 - Ifges(4,6) * t210) / 0.2e1 + t232 * t49 + t233 * t51 + (t192 - t43) * t106 + t347 - t440 - (Ifges(4,1) * t209 + t109 - t371) * t210 / 0.2e1 - (-Ifges(4,2) * t210 + t162 + t206) * t209 / 0.2e1 + (t56 * t314 + t44 * t312 + (t10 * t273 - t13 * t276) * mrSges(7,3) + (t22 * t273 - t23 * t276) * mrSges(6,3) + t433 * t415 + t432 * t413 + t434 * t406 + t443 * t399 + t276 * t472) * t274 * qJD(5) + ((t274 * t64 + t277 * t63) * t209 + t431) * mrSges(5,3) + m(5) * (-pkin(3) * t117 + pkin(10) * t431) - m(5) * (t102 * t63 + t103 * t64 + t123 * t145) + t225 * t50 + t460 * t98 + t461 * t105 + t215 * t48 - t181 * (mrSges(4,1) * t210 + mrSges(4,2) * t209) - t103 * t156 - t102 * t157 - pkin(3) * t97 - t90 * t99 + t313 * t370 + (-t188 + t387) * t144 + (-t350 + t386) * t145 - t186 * (Ifges(5,5) * t210 + (Ifges(5,1) * t277 - t384) * t209) / 0.2e1 - t185 * (Ifges(5,6) * t210 + (-Ifges(5,2) * t274 + t383) * t209) / 0.2e1 + m(6) * (pkin(10) * t370 + t192 * t23 + t193 * t22 + t232 * t6 + t233 * t5) + (t277 * t115 + t274 * t367) * pkin(10) + t2 * (mrSges(7,2) * t277 - mrSges(7,3) * t354) + t5 * (mrSges(6,2) * t277 - mrSges(6,3) * t354) + t1 * (-mrSges(7,1) * t277 - mrSges(7,3) * t352) + t6 * (-mrSges(6,1) * t277 - mrSges(6,3) * t352) + (Ifges(5,5) * t274 + Ifges(5,6) * t277) * t403 + t9 * t311 * t274 + t277 * t84 / 0.2e1 + t161 * t400 + (-mrSges(5,1) * t277 + mrSges(5,2) * t274 - mrSges(4,1)) * t117 + (Ifges(5,2) * t277 + t384) * t417 + (Ifges(5,1) * t274 + t383) * t418 + t274 * t423 + t352 * t464 + t210 * t454 - t207 * (Ifges(5,3) * t210 + (Ifges(5,5) * t277 - Ifges(5,6) * t274) * t209) / 0.2e1; t429 + ((-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t186 - t111 / 0.2e1 + t466) * t185 + t426 * t186 + ((m(7) * t44 + t98) * t393 - t427) * qJD(5) + t266 * t31 + t327 + t254 * t48 - t255 * t50 + t446 * t396 + (-t273 * t49 + t276 * t51 + (-m(6) * t297 - t273 * t106 - t276 * t108) * qJD(5) + m(6) * t430) * pkin(11) + t430 * mrSges(6,3) + t432 * t425 + t433 * t424 + t434 * t416 + (t237 - t39) * t105 + (-pkin(4) * t26 - t22 * t45 - t23 * t46 - t56 * t64) * m(6) + (t238 - t27) * t107 - t63 * t156 - t45 * t108 - t46 * t106 - t55 * t98 + (-t1 * t273 + t2 * t276) * mrSges(7,3) - t366 * t64 - pkin(4) * t32 + t273 * t464 - t9 * t312 - t26 * t314 + m(7) * (t1 * t254 + t10 * t238 + t13 * t237 - t2 * t255 + t266 * t9) - m(7) * (t10 * t27 + t13 * t39 + t44 * t55); t447 - t56 * (mrSges(6,1) * t155 + mrSges(6,2) * t154) - t44 * (mrSges(7,1) * t155 + mrSges(7,2) * t154) + (-(-t10 + t12) * t13 + (-t155 * t44 + t1) * pkin(5)) * m(7) + t23 * t108 - t12 * t105 - t22 * t106 + t13 * t107 + (-t155 * t98 + t48) * pkin(5) + (t154 * t22 + t155 * t23) * mrSges(6,3) + (t10 * t154 + t13 * t155) * mrSges(7,3) + (t154 * t453 - t469) * t413 + t444 * t412 + (t154 * t451 - t155 * t449) * t406 + (-t155 * t450 + t443 + t471) * t415 + t456; -t154 * t105 + t155 * t107 + 0.2e1 * (t9 / 0.2e1 + t10 * t412 + t13 * t415) * m(7) + t31;];
tauc  = t15(:);
