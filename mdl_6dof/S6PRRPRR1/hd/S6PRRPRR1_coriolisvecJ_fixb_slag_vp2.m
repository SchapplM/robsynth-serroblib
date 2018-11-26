% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2018-11-23 15:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:14:23
% EndTime: 2018-11-23 15:14:32
% DurationCPUTime: 8.55s
% Computational Cost: add. (9973->516), mult. (26066->717), div. (0->0), fcn. (20300->12), ass. (0->247)
t217 = sin(qJ(3));
t319 = -qJ(4) - pkin(8);
t266 = qJD(3) * t319;
t221 = cos(qJ(3));
t281 = qJD(4) * t221;
t173 = t217 * t266 + t281;
t282 = qJD(4) * t217;
t174 = t221 * t266 - t282;
t211 = sin(pkin(12));
t213 = cos(pkin(12));
t187 = t211 * t221 + t213 * t217;
t222 = cos(qJ(2));
t212 = sin(pkin(6));
t290 = qJD(1) * t212;
t272 = t222 * t290;
t348 = -t173 * t211 + t174 * t213 + t187 * t272;
t186 = -t211 * t217 + t213 * t221;
t347 = t173 * t213 + t174 * t211 - t186 * t272;
t179 = t187 * qJD(3);
t360 = pkin(9) * t179 - t347;
t180 = t186 * qJD(3);
t359 = -pkin(9) * t180 + t348;
t216 = sin(qJ(5));
t220 = cos(qJ(5));
t196 = t319 * t217;
t197 = t319 * t221;
t146 = t196 * t213 + t197 * t211;
t113 = -pkin(9) * t187 + t146;
t147 = t196 * t211 - t197 * t213;
t114 = pkin(9) * t186 + t147;
t248 = t113 * t220 - t114 * t216;
t352 = qJD(5) * t248 + t216 * t359 - t220 * t360;
t284 = qJD(3) * t217;
t209 = pkin(3) * t284;
t149 = pkin(4) * t179 + t209;
t218 = sin(qJ(2));
t273 = t218 * t290;
t244 = t186 * t220 - t187 * t216;
t85 = qJD(5) * t244 - t179 * t216 + t180 * t220;
t133 = t186 * t216 + t187 * t220;
t86 = qJD(5) * t133 + t179 * t220 + t180 * t216;
t358 = pkin(5) * t86 - pkin(10) * t85 + t149 - t273;
t206 = -pkin(3) * t221 - pkin(2);
t167 = qJD(2) * t206 + qJD(4) - t272;
t177 = t186 * qJD(2);
t129 = -pkin(4) * t177 + t167;
t210 = qJD(3) + qJD(5);
t215 = sin(qJ(6));
t219 = cos(qJ(6));
t285 = qJD(2) * t221;
t287 = qJD(2) * t217;
t178 = -t211 * t285 - t213 * t287;
t322 = pkin(9) * t178;
t192 = qJD(2) * pkin(8) + t273;
t264 = qJ(4) * qJD(2) + t192;
t214 = cos(pkin(6));
t289 = qJD(1) * t217;
t271 = t214 * t289;
t143 = t221 * t264 + t271;
t134 = t211 * t143;
t288 = qJD(1) * t221;
t202 = t214 * t288;
t142 = -t217 * t264 + t202;
t138 = qJD(3) * pkin(3) + t142;
t88 = t138 * t213 - t134;
t67 = qJD(3) * pkin(4) + t322 + t88;
t323 = pkin(9) * t177;
t291 = t213 * t143;
t89 = t138 * t211 + t291;
t79 = t89 + t323;
t36 = t216 * t67 + t220 * t79;
t34 = pkin(10) * t210 + t36;
t245 = t177 * t216 - t178 * t220;
t265 = t177 * t220 + t178 * t216;
t54 = -pkin(5) * t265 - pkin(10) * t245 + t129;
t16 = -t215 * t34 + t219 * t54;
t17 = t215 * t54 + t219 * t34;
t252 = t16 * t219 + t17 * t215;
t238 = t252 * mrSges(7,3);
t115 = Ifges(6,4) * t265;
t306 = t245 * Ifges(6,1);
t342 = t115 / 0.2e1 + t306 / 0.2e1;
t104 = t210 * t219 - t215 * t245;
t116 = qJD(6) - t265;
t255 = Ifges(7,5) * t219 - Ifges(7,6) * t215;
t313 = Ifges(7,4) * t219;
t257 = -Ifges(7,2) * t215 + t313;
t314 = Ifges(7,4) * t215;
t259 = Ifges(7,1) * t219 - t314;
t260 = mrSges(7,1) * t215 + mrSges(7,2) * t219;
t326 = t219 / 0.2e1;
t327 = -t215 / 0.2e1;
t35 = -t216 * t79 + t220 * t67;
t33 = -pkin(5) * t210 - t35;
t105 = t210 * t215 + t219 * t245;
t334 = t105 / 0.2e1;
t310 = t105 * Ifges(7,4);
t46 = Ifges(7,2) * t104 + Ifges(7,6) * t116 + t310;
t103 = Ifges(7,4) * t104;
t47 = Ifges(7,1) * t105 + Ifges(7,5) * t116 + t103;
t356 = t116 * t255 / 0.2e1 + t33 * t260 + t46 * t327 + t47 * t326 + t104 * t257 / 0.2e1 + t259 * t334;
t357 = -t129 * mrSges(6,2) - t210 * Ifges(6,5) + t238 - t342 - t356;
t268 = Ifges(4,5) * qJD(3) / 0.2e1;
t155 = -pkin(4) * t186 + t206;
t64 = -pkin(5) * t244 - pkin(10) * t133 + t155;
t66 = t113 * t216 + t114 * t220;
t32 = t215 * t64 + t219 * t66;
t355 = -qJD(6) * t32 - t215 * t352 + t219 * t358;
t31 = -t215 * t66 + t219 * t64;
t354 = qJD(6) * t31 + t215 * t358 + t219 * t352;
t353 = -qJD(5) * t66 + t216 * t360 + t220 * t359;
t299 = mrSges(6,1) * t210 + mrSges(7,1) * t104 - mrSges(7,2) * t105 - mrSges(6,3) * t245;
t165 = qJD(2) * t179;
t166 = qJD(2) * t180;
t112 = mrSges(5,1) * t165 + mrSges(5,2) * t166;
t77 = qJD(5) * t265 - t165 * t216 + t166 * t220;
t78 = qJD(5) * t245 + t165 * t220 + t166 * t216;
t38 = mrSges(6,1) * t78 + mrSges(6,2) * t77;
t351 = t112 + t38;
t205 = pkin(3) * t213 + pkin(4);
t324 = pkin(3) * t211;
t172 = t205 * t216 + t220 * t324;
t91 = -t142 * t211 - t291;
t243 = t91 - t323;
t93 = t142 * t213 - t134;
t80 = t93 + t322;
t350 = qJD(5) * t172 - t216 * t80 + t220 * t243;
t128 = -mrSges(5,1) * t177 - mrSges(5,2) * t178;
t349 = -m(5) * t167 - t128;
t346 = -t192 * t217 + t202;
t345 = -t115 / 0.2e1 + t357;
t344 = -t16 * t215 + t17 * t219;
t208 = pkin(3) * t287;
t286 = qJD(2) * t218;
t270 = t212 * t286;
t175 = qJD(1) * t270 + qJD(3) * t208;
t127 = pkin(4) * t165 + t175;
t29 = pkin(5) * t78 - pkin(10) * t77 + t127;
t151 = t192 * t221 + t271;
t283 = qJD(3) * t221;
t274 = qJ(4) * t283;
t275 = qJ(4) * t284;
t55 = (-t151 * t213 - t211 * t346) * qJD(3) + (-t211 * (t221 * t272 - t275 + t281) + t213 * (-t217 * t272 - t274 - t282)) * qJD(2);
t224 = -t166 * pkin(9) + t55;
t253 = qJD(4) + t272;
t56 = t213 * (t346 * qJD(3) + (t221 * t253 - t275) * qJD(2)) + t211 * (-t151 * qJD(3) + (-t217 * t253 - t274) * qJD(2));
t53 = -pkin(9) * t165 + t56;
t8 = qJD(5) * t35 + t216 * t224 + t220 * t53;
t2 = qJD(6) * t16 + t215 * t29 + t219 * t8;
t3 = -qJD(6) * t17 - t215 * t8 + t219 * t29;
t50 = qJD(6) * t104 + t219 * t77;
t51 = -qJD(6) * t105 - t215 * t77;
t343 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t50 + Ifges(7,6) * t51;
t82 = pkin(5) * t245 - pkin(10) * t265;
t341 = t50 / 0.2e1;
t340 = t51 / 0.2e1;
t339 = t78 / 0.2e1;
t9 = qJD(5) * t36 + t216 * t53 - t220 * t224;
t338 = t248 * t9;
t293 = t212 * t218;
t181 = t214 * t221 - t217 * t293;
t182 = t214 * t217 + t221 * t293;
t117 = t181 * t213 - t182 * t211;
t118 = t181 * t211 + t182 * t213;
t247 = t117 * t220 - t118 * t216;
t337 = t247 * t9;
t336 = -t104 / 0.2e1;
t335 = -t105 / 0.2e1;
t332 = -t116 / 0.2e1;
t330 = -t178 / 0.2e1;
t329 = -t179 / 0.2e1;
t328 = t180 / 0.2e1;
t321 = t2 * t219;
t320 = t3 * t215;
t317 = Ifges(4,4) * t217;
t316 = Ifges(5,4) * t178;
t312 = qJD(2) * pkin(2);
t307 = t265 * Ifges(6,2);
t296 = Ifges(4,6) * qJD(3);
t292 = t212 * t222;
t280 = qJD(6) * t215;
t279 = qJD(6) * t219;
t278 = qJD(2) * qJD(3);
t269 = qJD(2) * t292;
t267 = -t296 / 0.2e1;
t148 = -pkin(4) * t178 + t208;
t262 = -t2 * t215 - t219 * t3;
t261 = mrSges(7,1) * t219 - mrSges(7,2) * t215;
t258 = Ifges(7,1) * t215 + t313;
t256 = Ifges(7,2) * t219 + t314;
t254 = Ifges(7,5) * t215 + Ifges(7,6) * t219;
t23 = mrSges(7,1) * t78 - mrSges(7,3) * t50;
t24 = -mrSges(7,2) * t78 + mrSges(7,3) * t51;
t250 = -t215 * t23 + t219 * t24;
t62 = -mrSges(7,2) * t116 + mrSges(7,3) * t104;
t63 = mrSges(7,1) * t116 - mrSges(7,3) * t105;
t249 = -t215 * t63 + t219 * t62;
t75 = t117 * t216 + t118 * t220;
t234 = qJD(3) * t214 + t269;
t125 = -t192 * t284 + t234 * t288;
t126 = -t192 * t283 - t234 * t289;
t246 = t125 * t221 - t126 * t217;
t171 = t205 * t220 - t216 * t324;
t106 = -mrSges(6,2) * t210 + mrSges(6,3) * t265;
t241 = -t106 - t249;
t59 = -t215 * t75 - t219 * t292;
t240 = t215 * t292 - t219 * t75;
t193 = -t272 - t312;
t236 = t151 * mrSges(4,3) + t296 / 0.2e1 + (t221 * Ifges(4,2) + t317) * qJD(2) / 0.2e1 - t193 * mrSges(4,1);
t207 = Ifges(4,4) * t285;
t235 = t193 * mrSges(4,2) + Ifges(4,1) * t287 / 0.2e1 + t207 / 0.2e1 + t268 - t346 * mrSges(4,3);
t230 = -qJD(6) * t252 - t320;
t12 = Ifges(7,4) * t50 + Ifges(7,2) * t51 + Ifges(7,6) * t78;
t13 = Ifges(7,1) * t50 + Ifges(7,4) * t51 + Ifges(7,5) * t78;
t229 = -t8 * mrSges(6,2) + mrSges(7,3) * t321 + t215 * t13 / 0.2e1 + t12 * t326 + t258 * t341 + t256 * t340 + t254 * t339 - Ifges(6,6) * t78 + Ifges(6,5) * t77 + (-mrSges(6,1) - t261) * t9 + t356 * qJD(6);
t228 = t17 * mrSges(7,2) - t116 * Ifges(7,3) - t105 * Ifges(7,5) - t104 * Ifges(7,6) + t210 * Ifges(6,6) + t307 / 0.2e1 + Ifges(6,4) * t245 - t129 * mrSges(6,1) - t16 * mrSges(7,1);
t227 = -t307 / 0.2e1 - t228;
t223 = qJD(2) ^ 2;
t195 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t285;
t194 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t287;
t185 = (mrSges(4,1) * t217 + mrSges(4,2) * t221) * t278;
t170 = Ifges(5,4) * t177;
t169 = pkin(10) + t172;
t168 = -pkin(5) - t171;
t160 = t171 * qJD(5);
t154 = qJD(3) * mrSges(5,1) + mrSges(5,3) * t178;
t153 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t177;
t141 = -qJD(3) * t182 - t217 * t269;
t140 = qJD(3) * t181 + t221 * t269;
t120 = -t178 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t170;
t119 = t177 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t316;
t92 = t140 * t213 + t141 * t211;
t90 = -t140 * t211 + t141 * t213;
t81 = -mrSges(6,1) * t265 + mrSges(6,2) * t245;
t71 = Ifges(7,3) * t78;
t57 = t148 + t82;
t40 = t216 * t243 + t220 * t80;
t26 = qJD(5) * t75 + t216 * t92 - t220 * t90;
t25 = qJD(5) * t247 + t216 * t90 + t220 * t92;
t22 = t215 * t82 + t219 * t35;
t21 = -t215 * t35 + t219 * t82;
t20 = t215 * t57 + t219 * t40;
t19 = -t215 * t40 + t219 * t57;
t18 = -mrSges(7,1) * t51 + mrSges(7,2) * t50;
t15 = qJD(6) * t240 - t215 * t25 + t219 * t270;
t14 = qJD(6) * t59 + t215 * t270 + t219 * t25;
t1 = [t25 * t106 + t14 * t62 + t140 * t195 + t141 * t194 + t15 * t63 + t92 * t153 + t90 * t154 - t247 * t18 + t59 * t23 - t240 * t24 - t299 * t26 + (-t247 * t77 - t75 * t78) * mrSges(6,3) + (-t117 * t166 - t118 * t165) * mrSges(5,3) + (-t181 * t221 - t182 * t217) * mrSges(4,3) * t278 + ((-mrSges(3,2) * t223 - t185 - t351) * t222 + (-mrSges(3,1) * t223 + (t128 + qJD(2) * (-mrSges(4,1) * t221 + mrSges(4,2) * t217) + t81) * qJD(2)) * t218) * t212 + m(7) * (t14 * t17 + t15 * t16 - t2 * t240 + t26 * t33 + t3 * t59 - t337) + m(6) * (t25 * t36 - t26 * t35 - t337 + t75 * t8 + (-t127 * t222 + t129 * t286) * t212) + m(5) * (t117 * t55 + t118 * t56 + t88 * t90 + t89 * t92 + (t167 * t286 - t175 * t222) * t212) + m(4) * (t125 * t182 + t126 * t181 + t140 * t151 + t141 * t346 + (t193 - t272) * t270); (-t248 * t77 - t35 * t85 - t36 * t86 - t66 * t78) * mrSges(6,3) - t248 * t18 - (t127 * mrSges(6,1) - Ifges(6,4) * t77 + t71 / 0.2e1 - t8 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2)) * t78 + t343) * t244 + ((t217 * t194 - t221 * t195) * t222 + (-m(6) * t129 + t349 - t81) * t218 + ((-t151 * t221 + t217 * t346) * t222 + (-t312 - t193) * t218) * m(4)) * t290 + t347 * t153 + (t146 * t55 + t147 * t56 + t167 * t209 + t175 * t206 + t347 * t89 + t348 * t88) * m(5) + t348 * t154 + m(4) * ((-t151 * t217 - t221 * t346) * qJD(3) + t246) * pkin(8) + t227 * t86 + (t259 * t341 + t257 * t340 + t255 * t339 + t127 * mrSges(6,2) + Ifges(6,1) * t77 - Ifges(6,4) * t78 + t12 * t327 + t13 * t326 + (mrSges(6,3) + t260) * t9 + t262 * mrSges(7,3) + (-t219 * t46 / 0.2e1 + t47 * t327 + t254 * t332 + t258 * t335 + t256 * t336 + t33 * t261 - t344 * mrSges(7,3)) * qJD(6)) * t133 + t167 * (mrSges(5,1) * t179 + mrSges(5,2) * t180) + t246 * mrSges(4,3) + (0.3e1 / 0.2e1 * t221 ^ 2 - 0.3e1 / 0.2e1 * t217 ^ 2) * Ifges(4,4) * t278 + t206 * t112 - pkin(2) * t185 + t175 * (-mrSges(5,1) * t186 + mrSges(5,2) * t187) + t155 * t38 + t149 * t81 + t352 * t106 + t299 * t353 + (t127 * t155 + t129 * t149 + t35 * t353 + t352 * t36 + t66 * t8 - t338) * m(6) + t354 * t62 + t355 * t63 + (t16 * t355 + t17 * t354 + t2 * t32 + t3 * t31 - t33 * t353 - t338) * m(7) + t120 * t328 + t119 * t329 + (-t146 * t166 - t147 * t165 - t179 * t89 - t180 * t88 + t186 * t56 - t187 * t55) * mrSges(5,3) + (-t165 * t187 + t166 * t186 + t177 * t328 - t179 * t330) * Ifges(5,4) + (-t165 * t186 + t177 * t329) * Ifges(5,2) + (t342 - t357) * t85 + t31 * t23 + t32 * t24 + (Ifges(5,5) * t328 + Ifges(5,6) * t329 + (-pkin(8) * t194 + t235 + t268) * t221 + (-pkin(8) * t195 + pkin(3) * t128 + t267 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t285 - t236) * t217) * qJD(3) + (t166 * t187 + t180 * t330) * Ifges(5,1); -t346 * t195 + (-t306 / 0.2e1 + t345) * t265 + ((t268 - t207 / 0.2e1 - t235) * t221 + (t267 + (t317 / 0.2e1 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t221) * qJD(2) + t349 * pkin(3) + t236) * t217) * qJD(2) - t241 * t160 + (-t16 * t19 - t17 * t20 + (t230 + t321) * t169 + t168 * t9 + t350 * t33 + t344 * t160) * m(7) + (-t129 * t148 - t171 * t9 + t172 * t8 + (t160 - t40) * t36 - t350 * t35) * m(6) - t299 * t350 + t229 + ((t211 * t56 + t213 * t55) * pkin(3) - t88 * t91 - t89 * t93) * m(5) + ((-t215 * t62 - t219 * t63) * t169 - t238) * qJD(6) - (Ifges(5,2) * t178 + t120 + t170) * t177 / 0.2e1 - t227 * t245 + (-t171 * t77 - t172 * t78 + t245 * t36 + t265 * t35) * mrSges(6,3) + t250 * t169 + t151 * t194 - qJD(3) * (Ifges(5,5) * t177 + Ifges(5,6) * t178) / 0.2e1 - t167 * (-mrSges(5,1) * t178 + mrSges(5,2) * t177) + t168 * t18 - Ifges(5,6) * t165 + Ifges(5,5) * t166 - t93 * t153 - t91 * t154 - t148 * t81 + t126 * mrSges(4,1) - t125 * mrSges(4,2) - t40 * t106 + t119 * t330 + (t177 * t88 - t178 * t89 + (-t165 * t211 - t166 * t213) * pkin(3)) * mrSges(5,3) + t178 * (Ifges(5,1) * t177 + t316) / 0.2e1 - mrSges(7,3) * t320 + t55 * mrSges(5,1) - t56 * mrSges(5,2) - t20 * t62 - t19 * t63; t249 * qJD(6) + t299 * t245 + t241 * t265 - t177 * t153 - t178 * t154 + t215 * t24 + t219 * t23 + (t116 * t344 - t245 * t33 - t262) * m(7) + (t245 * t35 - t265 * t36 + t127) * m(6) + (-t177 * t89 - t178 * t88 + t175) * m(5) + t351; (mrSges(6,3) * t36 + t228) * t245 + (-t63 * t279 - t62 * t280 + t250) * pkin(10) + t229 + t230 * mrSges(7,3) + t299 * t36 + (t35 * mrSges(6,3) + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t245 + t345) * t265 - t35 * t106 - pkin(5) * t18 - t22 * t62 - t21 * t63 + (-t16 * t21 - t17 * t22 - t33 * t36 + (-t16 * t279 - t17 * t280 - t320 + t321) * pkin(10) - pkin(5) * t9) * m(7); t71 - t33 * (mrSges(7,1) * t105 + mrSges(7,2) * t104) + (Ifges(7,1) * t104 - t310) * t335 + t46 * t334 + (Ifges(7,5) * t104 - Ifges(7,6) * t105) * t332 - t16 * t62 + t17 * t63 + (t104 * t16 + t105 * t17) * mrSges(7,3) + (-Ifges(7,2) * t105 + t103 + t47) * t336 + t343;];
tauc  = t1(:);
