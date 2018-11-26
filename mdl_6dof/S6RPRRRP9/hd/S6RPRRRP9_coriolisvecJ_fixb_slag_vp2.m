% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2018-11-23 16:29
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:28:54
% EndTime: 2018-11-23 16:29:06
% DurationCPUTime: 11.60s
% Computational Cost: add. (7412->548), mult. (16570->741), div. (0->0), fcn. (10490->6), ass. (0->247)
t211 = sin(qJ(4));
t214 = cos(qJ(4));
t271 = qJD(3) * t214;
t215 = cos(qJ(3));
t273 = qJD(1) * t215;
t178 = -t211 * t273 + t271;
t179 = qJD(3) * t211 + t214 * t273;
t210 = sin(qJ(5));
t213 = cos(qJ(5));
t120 = t178 * t210 + t179 * t213;
t216 = -pkin(1) - pkin(7);
t200 = qJD(1) * t216 + qJD(2);
t280 = t200 * t215;
t169 = -qJD(3) * pkin(3) - t280;
t129 = -pkin(4) * t178 + t169;
t212 = sin(qJ(3));
t209 = t212 * qJD(1);
t205 = t209 + qJD(4);
t191 = pkin(3) * t212 - pkin(8) * t215 + qJ(2);
t162 = t191 * qJD(1);
t190 = t212 * t200;
t168 = qJD(3) * pkin(8) + t190;
t108 = t214 * t162 - t168 * t211;
t87 = -pkin(9) * t179 + t108;
t76 = pkin(4) * t205 + t87;
t109 = t162 * t211 + t168 * t214;
t88 = pkin(9) * t178 + t109;
t81 = t210 * t88;
t23 = t213 * t76 - t81;
t375 = qJ(6) * t120;
t18 = t23 - t375;
t198 = qJD(5) + t205;
t17 = pkin(5) * t198 + t18;
t83 = t213 * t88;
t24 = t210 * t76 + t83;
t242 = t213 * t178 - t179 * t210;
t350 = qJ(6) * t242;
t19 = t24 + t350;
t315 = -t198 / 0.2e1;
t326 = t120 / 0.2e1;
t327 = -t120 / 0.2e1;
t270 = qJD(3) * t215;
t244 = qJD(1) * t270;
t263 = qJD(3) * qJD(4);
t266 = qJD(4) * t215;
t134 = t214 * t263 + (-t211 * t266 - t212 * t271) * qJD(1);
t272 = qJD(3) * t212;
t372 = -t211 * t272 + t214 * t266;
t135 = -qJD(1) * t372 - t211 * t263;
t46 = qJD(5) * t242 + t134 * t213 + t135 * t210;
t241 = pkin(3) * t215 + pkin(8) * t212;
t176 = qJD(3) * t241 + qJD(2);
t149 = t176 * qJD(1);
t222 = -qJD(4) * t109 + t214 * t149;
t31 = -pkin(9) * t134 + (pkin(4) * qJD(1) - t200 * t211) * t270 + t222;
t252 = t200 * t270;
t267 = qJD(4) * t214;
t268 = qJD(4) * t211;
t58 = t211 * t149 + t162 * t267 - t168 * t268 + t214 * t252;
t35 = pkin(9) * t135 + t58;
t6 = -qJD(5) * t24 - t210 * t35 + t213 * t31;
t2 = pkin(5) * t244 - qJ(6) * t46 - qJD(6) * t120 + t6;
t47 = -qJD(5) * t120 - t134 * t210 + t135 * t213;
t264 = qJD(5) * t213;
t265 = qJD(5) * t210;
t5 = t210 * t31 + t213 * t35 + t76 * t264 - t265 * t88;
t3 = qJ(6) * t47 + qJD(6) * t242 + t5;
t343 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2);
t363 = Ifges(6,6) + Ifges(7,6);
t365 = Ifges(6,5) + Ifges(7,5);
t390 = Ifges(7,3) + Ifges(6,3);
t347 = t390 * t244 + t363 * t47 + t365 * t46;
t364 = Ifges(6,2) + Ifges(7,2);
t367 = Ifges(6,1) + Ifges(7,1);
t368 = -t242 / 0.2e1;
t366 = Ifges(6,4) + Ifges(7,4);
t389 = t366 * t242;
t379 = t367 * t120 + t365 * t198 + t389;
t388 = t120 * t366;
t380 = t363 * t198 + t364 * t242 + t388;
t69 = -pkin(5) * t242 + qJD(6) + t129;
t396 = t343 + t347 + (-t120 * t363 + t365 * t242) * t315 + (t120 * t19 + t17 * t242) * mrSges(7,3) + (t120 * t24 + t23 * t242) * mrSges(6,3) - t129 * (mrSges(6,1) * t120 + mrSges(6,2) * t242) - t69 * (mrSges(7,1) * t120 + mrSges(7,2) * t242) + t380 * t326 + (-t120 * t364 + t379 + t389) * t368 + (t367 * t242 - t388) * t327;
t187 = t241 * qJD(1);
t279 = t211 * t215;
t125 = t214 * t187 - t200 * t279;
t334 = -pkin(9) - pkin(8);
t254 = qJD(4) * t334;
t262 = pkin(9) * t212 * t214;
t395 = t214 * t254 - (pkin(4) * t215 + t262) * qJD(1) - t125;
t277 = t214 * t215;
t126 = t211 * t187 + t200 * t277;
t253 = t211 * t209;
t394 = pkin(9) * t253 - t211 * t254 + t126;
t182 = t210 * t214 + t211 * t213;
t345 = qJD(4) + qJD(5);
t124 = t345 * t182;
t160 = t182 * qJD(1);
t143 = t212 * t160;
t393 = t124 + t143;
t227 = t210 * t211 - t213 * t214;
t349 = t227 * t212;
t144 = qJD(1) * t349;
t374 = t345 * t227;
t392 = t374 + t144;
t391 = -t364 * t47 / 0.2e1 - t366 * t46 / 0.2e1 - t363 * t244 / 0.2e1;
t195 = t334 * t211;
t196 = t334 * t214;
t378 = t195 * t264 + t196 * t265 + t395 * t210 - t394 * t213;
t133 = t210 * t195 - t213 * t196;
t377 = -qJD(5) * t133 + t394 * t210 + t395 * t213;
t298 = Ifges(5,4) * t179;
t104 = t178 * Ifges(5,2) + t205 * Ifges(5,6) + t298;
t171 = Ifges(5,4) * t178;
t105 = t179 * Ifges(5,1) + t205 * Ifges(5,5) + t171;
t231 = t108 * t214 + t109 * t211;
t296 = Ifges(5,4) * t214;
t236 = -Ifges(5,2) * t211 + t296;
t297 = Ifges(5,4) * t211;
t238 = Ifges(5,1) * t214 - t297;
t239 = mrSges(5,1) * t211 + mrSges(5,2) * t214;
t294 = Ifges(5,6) * t211;
t295 = Ifges(5,5) * t214;
t309 = t214 / 0.2e1;
t311 = -t211 / 0.2e1;
t318 = t179 / 0.2e1;
t387 = -t231 * mrSges(5,3) + t104 * t311 + t105 * t309 + t169 * t239 + (-t294 + t295) * t205 / 0.2e1 + t236 * t178 / 0.2e1 + t238 * t318;
t385 = qJD(1) / 0.2e1;
t383 = t365 * t244 + t366 * t47 + t367 * t46;
t382 = -t393 * qJ(6) - qJD(6) * t227 + t378;
t381 = -pkin(5) * t273 + t392 * qJ(6) - qJD(6) * t182 + t377;
t145 = -pkin(4) * t253 + t190;
t376 = pkin(4) * t268 + t393 * pkin(5) - t145;
t284 = Ifges(4,5) * qJD(3);
t300 = Ifges(4,4) * t212;
t373 = t284 / 0.2e1 + (Ifges(4,1) * t215 - t300) * t385 + t387;
t314 = t198 / 0.2e1;
t329 = t242 / 0.2e1;
t247 = -Ifges(4,6) * qJD(3) / 0.2e1;
t308 = m(6) * t129;
t67 = -mrSges(6,1) * t242 + mrSges(6,2) * t120;
t354 = t67 + t308;
t153 = t227 * t215;
t353 = -qJD(3) * t153 - t124 * t212 - t160;
t151 = t182 * t215;
t352 = t227 * qJD(1) - qJD(3) * t151 + t212 * t374;
t351 = qJ(2) * (m(3) + m(4));
t175 = t214 * t191;
t243 = -t211 * t216 + pkin(4);
t116 = -pkin(9) * t277 + t212 * t243 + t175;
t278 = t212 * t216;
t197 = t214 * t278;
t140 = t211 * t191 + t197;
t128 = -pkin(9) * t279 + t140;
t65 = t210 * t116 + t213 * t128;
t348 = Ifges(5,5) * t134 + Ifges(5,6) * t135;
t59 = -t211 * t252 + t222;
t232 = -t211 * t59 + t214 * t58;
t346 = t59 * mrSges(5,1) - t58 * mrSges(5,2);
t137 = -mrSges(5,2) * t205 + mrSges(5,3) * t178;
t138 = mrSges(5,1) * t205 - mrSges(5,3) * t179;
t228 = -t211 * t137 - t214 * t138;
t344 = -m(5) * t231 + t228;
t257 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t258 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t259 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t299 = Ifges(4,4) * t215;
t342 = -t259 * t120 - t257 * t198 - t258 * t242 - t108 * mrSges(5,1) - t17 * mrSges(7,1) - t23 * mrSges(6,1) - t205 * Ifges(5,3) - t179 * Ifges(5,5) - t178 * Ifges(5,6) - t247 + (-t212 * Ifges(4,2) + t299) * t385 + t109 * mrSges(5,2) + t19 * mrSges(7,2) + t24 * mrSges(6,2) - t363 * t329 - t365 * t326 - t390 * t314;
t340 = m(5) / 0.2e1;
t339 = t46 / 0.2e1;
t338 = t47 / 0.2e1;
t325 = t134 / 0.2e1;
t324 = t135 / 0.2e1;
t321 = -t178 / 0.2e1;
t319 = -t179 / 0.2e1;
t313 = -t205 / 0.2e1;
t38 = -mrSges(7,2) * t244 + mrSges(7,3) * t47;
t39 = -mrSges(6,2) * t244 + mrSges(6,3) * t47;
t304 = t39 + t38;
t30 = t213 * t87 - t81;
t92 = -mrSges(7,2) * t198 + mrSges(7,3) * t242;
t93 = -mrSges(6,2) * t198 + mrSges(6,3) * t242;
t303 = t92 + t93;
t94 = mrSges(7,1) * t198 - mrSges(7,3) * t120;
t95 = mrSges(6,1) * t198 - mrSges(6,3) * t120;
t302 = t94 + t95;
t301 = mrSges(4,1) * t212;
t293 = qJ(2) * mrSges(4,1);
t292 = qJ(2) * mrSges(4,2);
t282 = qJD(3) * mrSges(4,2);
t274 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t178 - mrSges(5,2) * t179 - mrSges(4,3) * t273;
t269 = qJD(3) * t216;
t256 = t211 * t278;
t255 = Ifges(5,3) * t244 + t348;
t208 = -pkin(4) * t214 - pkin(3);
t180 = t200 * t272;
t250 = t215 * t269;
t248 = -t284 / 0.2e1;
t15 = -t47 * mrSges(7,1) + t46 * mrSges(7,2);
t101 = -pkin(4) * t135 + t180;
t29 = -t210 * t87 - t83;
t64 = t213 * t116 - t128 * t210;
t132 = t213 * t195 + t196 * t210;
t177 = pkin(4) * t279 - t215 * t216;
t240 = mrSges(5,1) * t214 - mrSges(5,2) * t211;
t237 = Ifges(5,1) * t211 + t296;
t235 = Ifges(5,2) * t214 + t297;
t233 = Ifges(5,5) * t211 + Ifges(5,6) * t214;
t230 = t108 * t211 - t109 * t214;
t114 = mrSges(5,1) * t244 - mrSges(5,3) * t134;
t115 = -mrSges(5,2) * t244 + mrSges(5,3) * t135;
t229 = -t211 * t114 + t214 * t115;
t136 = pkin(4) * t372 + t212 * t269;
t155 = t214 * t176;
t62 = t155 + (-t197 + (pkin(9) * t215 - t191) * t211) * qJD(4) + (t215 * t243 + t262) * qJD(3);
t85 = -qJD(4) * t256 + t211 * t176 + t191 * t267 + t214 * t250;
t68 = -pkin(9) * t372 + t85;
t9 = t116 * t264 - t128 * t265 + t210 * t62 + t213 * t68;
t10 = -qJD(5) * t65 - t210 * t68 + t213 * t62;
t207 = pkin(4) * t213 + pkin(5);
t193 = -mrSges(4,3) * t209 - t282;
t186 = (mrSges(4,2) * t215 + t301) * qJD(1);
t150 = t182 * t212;
t146 = pkin(5) * t227 + t208;
t139 = t175 - t256;
t122 = pkin(5) * t151 + t177;
t100 = -qJ(6) * t227 + t133;
t99 = -qJ(6) * t182 + t132;
t89 = pkin(4) * t179 + pkin(5) * t120;
t86 = -qJD(4) * t140 - t211 * t250 + t155;
t84 = -mrSges(5,1) * t135 + mrSges(5,2) * t134;
t78 = t134 * Ifges(5,1) + t135 * Ifges(5,4) + Ifges(5,5) * t244;
t77 = t134 * Ifges(5,4) + t135 * Ifges(5,2) + Ifges(5,6) * t244;
t75 = t182 * t272 + t215 * t374;
t73 = qJD(3) * t349 - t124 * t215;
t66 = -mrSges(7,1) * t242 + mrSges(7,2) * t120;
t51 = -pkin(5) * t75 + t136;
t48 = -qJ(6) * t151 + t65;
t40 = pkin(5) * t212 + qJ(6) * t153 + t64;
t37 = mrSges(6,1) * t244 - mrSges(6,3) * t46;
t36 = mrSges(7,1) * t244 - mrSges(7,3) * t46;
t22 = -pkin(5) * t47 + t101;
t21 = t30 - t375;
t20 = t29 - t350;
t16 = -mrSges(6,1) * t47 + mrSges(6,2) * t46;
t8 = qJ(6) * t75 - qJD(6) * t151 + t9;
t7 = pkin(5) * t270 - qJ(6) * t73 + qJD(6) * t153 + t10;
t1 = [m(5) * (t108 * t86 + t109 * t85 + t139 * t59 + t140 * t58) + t379 * t73 / 0.2e1 + t380 * t75 / 0.2e1 - t383 * t153 / 0.2e1 + (t248 + (-0.2e1 * t292 + 0.3e1 / 0.2e1 * t300) * qJD(1) + (m(5) * t169 - t274) * t216 - t373) * t272 + (qJD(1) * qJD(2) * mrSges(4,2) + t77 * t311 + t78 * t309 - t216 * t84 + t238 * t325 + t236 * t324 + (-t211 * t58 - t214 * t59) * mrSges(5,3) + (t105 * t311 - t214 * t104 / 0.2e1 + t233 * t313 + t235 * t321 + t237 * t319 + t169 * t240 + t230 * mrSges(5,3)) * qJD(4) + (t216 * t193 + ((-0.3e1 / 0.2e1 * Ifges(4,4) + t295 / 0.2e1 - t294 / 0.2e1) * t215 + 0.2e1 * t293 - t259 * t153 - t258 * t151) * qJD(1) + (-m(5) * t216 + t239) * t190 + t247 - t342) * qJD(3)) * t215 + m(7) * (t122 * t22 + t17 * t7 + t19 * t8 + t2 * t40 + t3 * t48 + t51 * t69) + m(6) * (t10 * t23 + t101 * t177 + t129 * t136 + t24 * t9 + t5 * t65 + t6 * t64) + (t363 * t75 + t365 * t73) * t314 + ((-0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2) + Ifges(5,3) / 0.2e1 + t257) * t244 + t365 * t339 + t363 * t338 + t343 + t346) * t212 + (t364 * t75 + t366 * t73) * t329 + (-t151 * t364 - t153 * t366) * t338 + (t366 * t75 + t367 * t73) * t326 + (-t151 * t366 - t153 * t367) * t339 + (-t151 * t5 + t153 * t6 - t23 * t73 + t24 * t75) * mrSges(6,3) + t101 * (mrSges(6,1) * t151 - mrSges(6,2) * t153) + t22 * (mrSges(7,1) * t151 - mrSges(7,2) * t153) + (-t151 * t3 + t153 * t2 - t17 * t73 + t19 * t75) * mrSges(7,3) + (t186 + ((2 * mrSges(3,3)) + t301 + 0.2e1 * t351) * qJD(1)) * qJD(2) + (t255 + t347 + t348) * t212 / 0.2e1 + t177 * t16 + t136 * t67 + t85 * t137 + t86 * t138 + t139 * t114 + t140 * t115 + t129 * (-mrSges(6,1) * t75 + mrSges(6,2) * t73) + t122 * t15 + t8 * t92 + t9 * t93 + t7 * t94 + t10 * t95 + t69 * (-mrSges(7,1) * t75 + mrSges(7,2) * t73) + t65 * t39 + t51 * t66 + t64 * t37 + t48 * t38 + t40 * t36 + t151 * t391; -t304 * t349 - (t36 + t37) * t150 + (-t15 - t16 - t84 + (t137 * t214 - t138 * t211 + t193) * qJD(3)) * t215 + (t228 * qJD(4) + (t66 + t67 - t274) * qJD(3) + t229) * t212 - m(5) * t230 * t270 + 0.2e1 * ((-qJD(4) * t231 + t232) * t340 + (t308 / 0.2e1 + (t169 - t280) * t340) * qJD(3)) * t212 + t353 * t303 + t352 * t302 + (-t150 * t2 + t17 * t352 + t19 * t353 - t215 * t22 + t272 * t69 - t3 * t349) * m(7) + (-t101 * t215 - t150 * t6 + t23 * t352 + t24 * t353 - t349 * t5) * m(6) + (-t186 + (-mrSges(3,3) - t351) * qJD(1) + t344) * qJD(1); (-pkin(3) * t180 + pkin(8) * t232 - t108 * t125 - t109 * t126 - t169 * t190) * m(5) + (t143 * t363 + t144 * t365) * t315 + (t143 * t364 + t144 * t366) * t368 + (t182 * t366 - t227 * t364) * t338 + (t143 * t366 + t144 * t367) * t327 + (t182 * t367 - t227 * t366) * t339 + t377 * t95 + (t101 * t208 - t129 * t145 + t132 * t6 + t133 * t5 + t23 * t377 + t24 * t378) * m(6) + t378 * t93 + t380 * (-t143 / 0.2e1 - t124 / 0.2e1) + t381 * t94 + (t100 * t3 + t146 * t22 + t17 * t381 + t19 * t382 + t2 * t99 + t376 * t69) * m(7) + t382 * t92 + t383 * t182 / 0.2e1 + t376 * t66 + ((t248 + (t292 - t300 / 0.2e1) * qJD(1) + t373) * t212 + ((-t293 + t299 / 0.2e1 + (Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1) * t212) * qJD(1) + t247 + t342) * t215 + (t182 * t365 - t227 * t363 + t233) * t270 / 0.2e1) * qJD(1) + ((-t193 - t282) * t215 + ((-mrSges(4,1) - t240) * qJD(3) + t274) * t212) * t200 + t229 * pkin(8) + t232 * mrSges(5,3) + t101 * (mrSges(6,1) * t227 + mrSges(6,2) * t182) + t22 * (mrSges(7,1) * t227 + mrSges(7,2) * t182) + (pkin(4) * t211 * t354 + pkin(8) * t344 + t387) * qJD(4) + (t17 * t392 - t182 * t2 - t19 * t393 - t227 * t3) * mrSges(7,3) + (mrSges(6,1) * t393 - mrSges(6,2) * t392) * t129 + (mrSges(7,1) * t393 - mrSges(7,2) * t392) * t69 + (-t182 * t6 - t227 * t5 + t23 * t392 - t24 * t393) * mrSges(6,3) + (-t124 * t364 - t366 * t374) * t329 + (-t124 * t366 - t367 * t374) * t326 + t379 * (-t144 / 0.2e1 - t374 / 0.2e1) + (-t124 * t363 - t365 * t374) * t314 + t211 * t78 / 0.2e1 + t208 * t16 - t145 * t67 + t146 * t15 - t126 * t137 - t125 * t138 + t132 * t37 + t133 * t39 + t99 * t36 + t100 * t38 - pkin(3) * t84 + t77 * t309 + t235 * t324 + t237 * t325 + t227 * t391; t255 + (-Ifges(5,2) * t179 + t105 + t171) * t321 - m(6) * (t23 * t29 + t24 * t30) + t346 + (t108 * t178 + t109 * t179) * mrSges(5,3) + (t213 * t37 + t304 * t210 + (-t210 * t302 + t213 * t303) * qJD(5) + m(6) * (t210 * t5 + t213 * t6 - t23 * t265 + t24 * t264) - t354 * t179) * pkin(4) + (t2 * t207 - t17 * t20 - t19 * t21 - t69 * t89 + (-t17 * t265 + t19 * t264 + t210 * t3) * pkin(4)) * m(7) + t207 * t36 - t169 * (mrSges(5,1) * t179 + mrSges(5,2) * t178) - t108 * t137 + t109 * t138 - t21 * t92 - t30 * t93 - t20 * t94 - t29 * t95 - t89 * t66 + (Ifges(5,5) * t178 - Ifges(5,6) * t179) * t313 + t104 * t318 + (Ifges(5,1) * t178 - t298) * t319 + t396; (-(-t17 + t18) * t19 + (-t120 * t69 + t2) * pkin(5)) * m(7) + (-t120 * t66 + t36) * pkin(5) - t18 * t92 - t23 * t93 + t19 * t94 + t24 * t95 + t396; -t242 * t92 + t120 * t94 + 0.2e1 * (t22 / 0.2e1 + t19 * t368 + t17 * t326) * m(7) + t15;];
tauc  = t1(:);
