% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:30:16
% EndTime: 2019-03-08 22:30:34
% DurationCPUTime: 8.04s
% Computational Cost: add. (5848->560), mult. (13806->783), div. (0->0), fcn. (9034->10), ass. (0->281)
t203 = sin(qJ(3));
t204 = sin(qJ(2));
t206 = cos(qJ(5));
t199 = sin(pkin(6));
t278 = qJD(1) * t199;
t202 = sin(qJ(5));
t208 = cos(qJ(2));
t287 = t202 * t208;
t115 = (t203 * t287 + t204 * t206) * t278;
t273 = qJD(3) * t203;
t193 = pkin(3) * t273;
t207 = cos(qJ(3));
t232 = pkin(9) * t203 - qJ(4) * t207;
t270 = qJD(4) * t203;
t213 = qJD(3) * t232 - t270;
t121 = t193 + t213;
t209 = -pkin(3) - pkin(9);
t247 = -qJ(4) * t203 - pkin(2);
t143 = t207 * t209 + t247;
t271 = qJD(3) * t207;
t329 = pkin(4) + pkin(8);
t158 = t329 * t271;
t171 = t329 * t203;
t268 = qJD(5) * t206;
t269 = qJD(5) * t202;
t37 = t206 * t121 - t143 * t269 + t202 * t158 + t171 * t268;
t374 = -t115 + t37;
t284 = t206 * t208;
t114 = (-t202 * t204 + t203 * t284) * t278;
t152 = t202 * t171;
t223 = -pkin(10) * t202 * t203 + pkin(5) * t207;
t244 = -t121 * t202 + t206 * t158;
t246 = pkin(10) * t207 - t143;
t373 = -t114 + t223 * qJD(3) + (t206 * t246 - t152) * qJD(5) + t244;
t267 = qJD(5) * t207;
t250 = t202 * t267;
t272 = qJD(3) * t206;
t214 = t203 * t272 + t250;
t372 = -pkin(10) * t214 - t374;
t312 = pkin(10) - t209;
t255 = t204 * t278;
t160 = qJD(2) * pkin(8) + t255;
t200 = cos(pkin(6));
t277 = qJD(1) * t203;
t177 = t200 * t277;
t119 = t207 * t160 + t177;
t274 = qJD(2) * t207;
t102 = pkin(4) * t274 + t119;
t194 = t203 * qJD(2);
t191 = pkin(3) * t194;
t125 = qJD(2) * t232 + t191;
t57 = t206 * t102 - t125 * t202;
t371 = -qJD(2) * t223 + t312 * t269 - t57;
t163 = t312 * t206;
t251 = t206 * t194;
t58 = t202 * t102 + t206 * t125;
t370 = pkin(10) * t251 + qJD(5) * t163 + t58;
t254 = t208 * t278;
t105 = qJD(2) * t143 - t254;
t245 = pkin(4) * qJD(2) + t160;
t228 = t245 * t203;
t289 = t200 * t207;
t178 = qJD(1) * t289;
t362 = qJD(4) - t178;
t80 = qJD(3) * t209 + t228 + t362;
t42 = -t105 * t202 + t206 * t80;
t43 = t105 * t206 + t202 * t80;
t230 = t202 * t42 - t206 * t43;
t310 = Ifges(6,4) * t202;
t234 = Ifges(6,2) * t206 + t310;
t309 = Ifges(6,4) * t206;
t236 = Ifges(6,1) * t202 + t309;
t239 = mrSges(6,1) * t206 - mrSges(6,2) * t202;
t305 = Ifges(6,6) * t206;
t308 = Ifges(6,5) * t202;
t316 = -t206 / 0.2e1;
t317 = -t202 / 0.2e1;
t185 = t194 + qJD(5);
t318 = -t185 / 0.2e1;
t149 = -t202 * t274 + t272;
t323 = -t149 / 0.2e1;
t148 = -qJD(3) * t202 - t206 * t274;
t324 = -t148 / 0.2e1;
t311 = Ifges(6,4) * t149;
t73 = Ifges(6,2) * t148 + Ifges(6,6) * t185 + t311;
t144 = Ifges(6,4) * t148;
t74 = Ifges(6,1) * t149 + Ifges(6,5) * t185 + t144;
t198 = qJD(3) * qJ(4);
t88 = t198 + t102;
t369 = t230 * mrSges(6,3) + (t305 + t308) * t318 + t234 * t324 + t236 * t323 + t88 * t239 + t316 * t73 + t317 * t74;
t205 = cos(qJ(6));
t32 = -pkin(10) * t149 + t42;
t24 = pkin(5) * t185 + t32;
t201 = sin(qJ(6));
t33 = pkin(10) * t148 + t43;
t296 = t201 * t33;
t12 = t205 * t24 - t296;
t295 = t205 * t33;
t13 = t201 * t24 + t295;
t265 = qJD(2) * qJD(3);
t248 = t207 * t265;
t182 = Ifges(7,3) * t248;
t243 = t205 * t148 - t149 * t201;
t87 = t148 * t201 + t149 * t205;
t315 = Ifges(7,4) * t87;
t263 = qJD(5) + qJD(6);
t175 = t194 + t263;
t320 = -t175 / 0.2e1;
t331 = -t87 / 0.2e1;
t333 = -t243 / 0.2e1;
t81 = Ifges(7,4) * t243;
t36 = Ifges(7,1) * t87 + Ifges(7,5) * t175 + t81;
t63 = -pkin(5) * t148 + t88;
t368 = t182 + (Ifges(7,5) * t243 - Ifges(7,6) * t87) * t320 + (t12 * t243 + t13 * t87) * mrSges(7,3) + (-Ifges(7,2) * t87 + t36 + t81) * t333 - t63 * (mrSges(7,1) * t87 + mrSges(7,2) * t243) + (Ifges(7,1) * t243 - t315) * t331;
t153 = t206 * t171;
t67 = pkin(5) * t203 + t202 * t246 + t153;
t285 = t206 * t207;
t92 = t206 * t143 + t152;
t75 = -pkin(10) * t285 + t92;
t31 = t201 * t67 + t205 * t75;
t367 = -qJD(6) * t31 + t372 * t201 + t205 * t373;
t30 = -t201 * t75 + t205 * t67;
t366 = qJD(6) * t30 + t201 * t373 - t372 * t205;
t355 = mrSges(5,1) + mrSges(4,3);
t365 = mrSges(5,2) - mrSges(4,1);
t108 = -t198 - t119;
t164 = -pkin(3) * t207 + t247;
t120 = qJD(2) * t164 - t254;
t304 = qJD(2) * pkin(2);
t161 = -t254 - t304;
t356 = qJD(3) / 0.2e1;
t357 = -qJD(3) / 0.2e1;
t358 = -qJD(2) / 0.2e1;
t364 = t119 * mrSges(4,3) + t120 * mrSges(5,2) + Ifges(4,6) * t356 + (Ifges(4,4) * t203 + t207 * Ifges(4,2)) * qJD(2) / 0.2e1 + Ifges(5,5) * t357 + (-Ifges(5,6) * t203 - t207 * Ifges(5,3)) * t358 - t108 * mrSges(5,1) - t161 * mrSges(4,1) + t369;
t264 = qJD(3) * qJD(5);
t113 = qJD(2) * t214 - t206 * t264;
t276 = qJD(2) * t199;
t249 = qJD(1) * t276;
t242 = t208 * t249;
t66 = t203 * t242 + (t207 * t245 + t177) * qJD(3);
t279 = qJD(3) * t191 + t204 * t249;
t78 = qJD(2) * t213 + t279;
t17 = -t105 * t269 + t202 * t66 + t206 * t78 + t80 * t268;
t11 = pkin(10) * t113 + t17;
t112 = -t202 * t264 + (t202 * t273 - t206 * t267) * qJD(2);
t18 = -qJD(5) * t43 - t202 * t78 + t206 * t66;
t8 = pkin(5) * t248 - pkin(10) * t112 + t18;
t2 = qJD(6) * t12 + t11 * t205 + t201 * t8;
t28 = qJD(6) * t243 + t112 * t205 + t113 * t201;
t29 = -qJD(6) * t87 - t112 * t201 + t113 * t205;
t344 = qJD(6) * t13;
t3 = -t11 * t201 + t205 * t8 - t344;
t363 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t28 + Ifges(7,6) * t29;
t286 = t205 * t206;
t288 = t201 * t202;
t224 = -t286 + t288;
t225 = t201 * t206 + t205 * t202;
t216 = t203 * t225;
t124 = qJD(2) * t216;
t95 = t263 * t225;
t292 = -t95 - t124;
t123 = -t194 * t288 + t205 * t251;
t266 = qJD(6) * t201;
t94 = -t201 * t269 - t202 * t266 + t263 * t286;
t293 = t94 + t123;
t361 = -t12 * t292 - t13 * t293 - t2 * t225 + t224 * t3;
t35 = Ifges(7,2) * t243 + Ifges(7,6) * t175 + t315;
t359 = t35 / 0.2e1;
t162 = t312 * t202;
t103 = t162 * t201 - t163 * t205;
t354 = qJD(6) * t103 + t201 * t371 - t370 * t205;
t104 = -t162 * t205 - t163 * t201;
t353 = -qJD(6) * t104 + t370 * t201 + t205 * t371;
t256 = -pkin(5) * t206 - pkin(4);
t346 = pkin(5) * t268 - (qJD(2) * t256 - t160) * t203 + t362;
t44 = -mrSges(7,1) * t243 + mrSges(7,2) * t87;
t345 = -m(7) * t63 - t44;
t343 = t17 * t202 + t18 * t206;
t118 = t160 * t203 - t178;
t342 = -qJD(4) - t118;
t340 = t207 * t263;
t339 = t18 * mrSges(6,1) - t17 * mrSges(6,2) + Ifges(6,5) * t112 + Ifges(6,6) * t113 + t363;
t337 = -m(5) / 0.2e1;
t336 = t28 / 0.2e1;
t335 = t29 / 0.2e1;
t334 = t73 / 0.2e1;
t332 = t243 / 0.2e1;
t330 = t87 / 0.2e1;
t126 = t224 * t207;
t326 = t126 / 0.2e1;
t127 = t225 * t207;
t325 = -t127 / 0.2e1;
t322 = -t225 / 0.2e1;
t321 = -t224 / 0.2e1;
t319 = t175 / 0.2e1;
t307 = Ifges(6,5) * t206;
t306 = Ifges(6,6) * t202;
t290 = t199 * t204;
t257 = t203 * t290;
t132 = t257 - t289;
t252 = t208 * t276;
t218 = qJD(3) * t200 + t252;
t77 = t160 * t271 + t218 * t277;
t303 = t132 * t77;
t169 = -mrSges(5,1) * t274 - qJD(3) * mrSges(5,3);
t96 = -mrSges(6,1) * t148 + mrSges(6,2) * t149;
t294 = -t169 + t96;
t154 = (mrSges(5,2) * t207 - mrSges(5,3) * t203) * qJD(2);
t283 = t154 + (-mrSges(4,1) * t207 + mrSges(4,2) * t203) * qJD(2);
t282 = qJD(3) * t178 + t207 * t242;
t168 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t274;
t281 = t168 - t169;
t280 = qJD(3) * t365 + t194 * t355;
t172 = t329 * t207;
t275 = qJD(2) * t204;
t262 = pkin(8) * t77 / 0.2e1;
t261 = -t44 - t294;
t260 = -0.3e1 / 0.2e1 * Ifges(4,4) - 0.3e1 / 0.2e1 * Ifges(5,6);
t259 = Ifges(4,5) / 0.2e1 - Ifges(5,4) / 0.2e1;
t258 = -Ifges(4,6) / 0.2e1 + Ifges(5,5) / 0.2e1;
t253 = t199 * t275;
t241 = -t254 / 0.2e1;
t238 = mrSges(6,1) * t202 + mrSges(6,2) * t206;
t237 = Ifges(6,1) * t206 - t310;
t235 = -Ifges(6,2) * t202 + t309;
t82 = mrSges(6,1) * t248 - mrSges(6,3) * t112;
t83 = -mrSges(6,2) * t248 + mrSges(6,3) * t113;
t229 = t202 * t83 + t206 * t82;
t221 = -t132 * t202 + t199 * t284;
t99 = t132 * t206 + t199 * t287;
t49 = t201 * t99 - t205 * t221;
t48 = t201 * t221 + t205 * t99;
t116 = -mrSges(6,2) * t185 + mrSges(6,3) * t148;
t117 = mrSges(6,1) * t185 - mrSges(6,3) * t149;
t227 = t206 * t116 - t202 * t117;
t76 = -t160 * t273 + t282;
t133 = t200 * t203 + t207 * t290;
t219 = -qJ(4) * t271 - t270;
t217 = -m(4) * t119 + m(5) * t108 - t281;
t197 = qJD(3) * qJD(4);
t60 = -qJD(3) * t228 + t197 + t282;
t107 = -qJD(3) * pkin(3) - t342;
t190 = Ifges(4,4) * t274;
t212 = t107 * mrSges(5,1) + t118 * mrSges(4,3) + t12 * mrSges(7,1) + t161 * mrSges(4,2) + t42 * mrSges(6,1) + Ifges(4,1) * t194 / 0.2e1 + Ifges(4,5) * t356 + t190 / 0.2e1 + Ifges(5,4) * t357 + (-t203 * Ifges(5,2) - Ifges(5,6) * t207) * t358 + t175 * Ifges(7,3) + t87 * Ifges(7,5) + t243 * Ifges(7,6) + t185 * Ifges(6,3) + t149 * Ifges(6,5) + t148 * Ifges(6,6) - t120 * mrSges(5,3) - t13 * mrSges(7,2) - t43 * mrSges(6,2);
t210 = qJD(2) ^ 2;
t187 = pkin(5) * t202 + qJ(4);
t183 = Ifges(6,3) * t248;
t157 = t329 * t273;
t156 = -qJ(4) * t274 + t191;
t141 = (mrSges(4,1) * t203 + mrSges(4,2) * t207) * t265;
t140 = (-mrSges(5,2) * t203 - mrSges(5,3) * t207) * t265;
t131 = pkin(5) * t285 + t172;
t129 = t193 + t219;
t106 = -pkin(5) * t250 + (-pkin(8) + t256) * t273;
t101 = t178 - t228;
t98 = -qJD(3) * t257 + t207 * t218;
t97 = qJD(3) * t133 + t203 * t252;
t93 = qJD(2) * t219 + t279;
t91 = -t143 * t202 + t153;
t68 = -t197 - t76;
t65 = mrSges(7,1) * t175 - mrSges(7,3) * t87;
t64 = -mrSges(7,2) * t175 + mrSges(7,3) * t243;
t59 = -mrSges(6,1) * t113 + mrSges(6,2) * t112;
t54 = t112 * Ifges(6,1) + t113 * Ifges(6,4) + Ifges(6,5) * t248;
t53 = t112 * Ifges(6,4) + t113 * Ifges(6,2) + Ifges(6,6) * t248;
t52 = -t224 * t273 + t225 * t340;
t51 = qJD(3) * t216 + t224 * t340;
t41 = -pkin(5) * t113 + t60;
t40 = qJD(5) * t99 + t202 * t97 + t206 * t253;
t39 = qJD(5) * t221 - t202 * t253 + t206 * t97;
t38 = -qJD(5) * t92 + t244;
t22 = -mrSges(7,2) * t248 + mrSges(7,3) * t29;
t21 = mrSges(7,1) * t248 - mrSges(7,3) * t28;
t16 = t205 * t32 - t296;
t15 = -t201 * t32 - t295;
t14 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t10 = t28 * Ifges(7,1) + t29 * Ifges(7,4) + Ifges(7,5) * t248;
t9 = t28 * Ifges(7,4) + t29 * Ifges(7,2) + Ifges(7,6) * t248;
t7 = -qJD(6) * t49 - t201 * t40 + t205 * t39;
t6 = qJD(6) * t48 + t201 * t39 + t205 * t40;
t1 = [-t221 * t83 + t40 * t116 + t39 * t117 + t48 * t21 + t49 * t22 + t6 * t64 + t7 * t65 + t99 * t82 + t280 * t97 + (t59 + t14) * t133 + (t168 - t261) * t98 + (-t210 * t204 * mrSges(3,1) + (-mrSges(3,2) * t210 - t140 - t141) * t208) * t199 + (t283 * t290 + t355 * qJD(3) * (t132 * t207 - t133 * t203)) * qJD(2) + m(4) * (t118 * t97 + t119 * t98 + t303 + t133 * t76 + (t161 - t254) * t253) + m(5) * (t107 * t97 - t108 * t98 + t303 - t133 * t68 + (t120 * t275 - t208 * t93) * t199) + m(7) * (t12 * t7 + t13 * t6 + t133 * t41 + t2 * t49 + t3 * t48 + t63 * t98) + m(6) * (t133 * t60 - t17 * t221 + t18 * t99 + t39 * t42 + t40 * t43 + t88 * t98); (-t114 + t38) * t117 + (-t93 * mrSges(5,3) + t182 / 0.2e1 + t183 / 0.2e1 + t355 * t77 + (mrSges(4,2) * t275 - t208 * t280) * t278 + 0.2e1 * (t107 * t241 + t262) * m(5) + 0.2e1 * (t118 * t241 + t262) * m(4) + (pkin(8) * t217 + t258 * qJD(3) + t194 * t260 - t364) * qJD(3) + t339) * t203 + t366 * t64 + (t106 * t63 + t12 * t367 + t13 * t366 + t131 * t41 + t2 * t31 + t3 * t30) * m(7) + t367 * t65 + (-Ifges(7,4) * t127 + Ifges(7,2) * t126) * t335 + (-Ifges(7,1) * t127 + Ifges(7,4) * t126) * t336 + m(6) * (-t157 * t88 + t17 * t92 + t172 * t60 + t18 * t91 + t37 * t43 + t38 * t42) + (-t12 * t51 + t126 * t2 + t127 * t3 + t13 * t52) * mrSges(7,3) + t41 * (-mrSges(7,1) * t126 - mrSges(7,2) * t127) + t374 * t116 + m(5) * (t120 * t129 + t164 * t93) + t52 * t359 + 0.2e1 * (t120 * t337 + (-t304 / 0.2e1 - t161 / 0.2e1) * m(4)) * t255 - t283 * t255 + (-t112 * t236 / 0.2e1 - t113 * t234 / 0.2e1 - t68 * mrSges(5,1) + t76 * mrSges(4,3) + t54 * t317 + t53 * t316 + t60 * t239 + t93 * mrSges(5,2) + (-t17 * t206 + t18 * t202) * mrSges(6,3) + (-mrSges(4,1) * t275 + (-m(6) * t88 + t217 + t345 - t96) * t208) * t278 + (t185 * (t306 - t307) / 0.2e1 + t237 * t323 + t235 * t324 - t88 * t238 + t74 * t316 + t202 * t334 + (t43 * t202 + t42 * t206) * mrSges(6,3)) * qJD(5) + (((-t308 / 0.2e1 - t305 / 0.2e1 - t260) * t207 + Ifges(7,5) * t325 + Ifges(7,6) * t326 + (Ifges(6,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,2) + Ifges(7,3) / 0.2e1) * t203) * qJD(2) + t259 * qJD(3) + t212) * qJD(3) + (0.2e1 * t68 * t337 + m(4) * t76 + (m(4) * t118 + m(5) * t107 + t280) * qJD(3)) * pkin(8)) * t207 + t172 * t59 + t164 * t140 + t129 * t154 - t157 * t96 + t131 * t14 - pkin(2) * t141 + t106 * t44 + t91 * t82 + t92 * t83 + t63 * (-mrSges(7,1) * t52 + mrSges(7,2) * t51) + t51 * t36 / 0.2e1 - m(6) * (t114 * t42 + t115 * t43) + (Ifges(7,5) * t51 + Ifges(7,6) * t52) * t319 + t10 * t325 + t9 * t326 + (Ifges(7,1) * t51 + Ifges(7,4) * t52) * t330 + (Ifges(7,4) * t51 + Ifges(7,2) * t52) * t332 + t30 * t21 + t31 * t22; t369 * qJD(5) + t361 * mrSges(7,3) - t343 * mrSges(6,3) + (-Ifges(7,1) * t95 - Ifges(7,4) * t94) * t330 + (-Ifges(7,4) * t95 - Ifges(7,2) * t94) * t332 + t365 * t77 + (((Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1) * t194 + (-qJ(4) * mrSges(5,1) + t258) * qJD(3) + t364) * t203 + ((-pkin(3) * mrSges(5,1) + Ifges(7,5) * t321 + Ifges(7,6) * t322 + t307 / 0.2e1 - t306 / 0.2e1 + t259) * qJD(3) - Ifges(5,6) * t274 / 0.2e1 - t190 / 0.2e1 + (-Ifges(5,2) / 0.2e1 + Ifges(5,3) / 0.2e1 - Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t194 - t212) * t207) * qJD(2) + (-pkin(3) * t77 - qJ(4) * t68 - t107 * t119 + t108 * t342 - t120 * t156) * m(5) + ((-m(6) * t230 + t227) * qJD(5) + m(6) * t343 + t229) * t209 + t41 * (mrSges(7,1) * t225 - mrSges(7,2) * t224) + (-Ifges(7,4) * t224 - Ifges(7,2) * t225) * t335 + (-Ifges(7,1) * t224 - Ifges(7,4) * t225) * t336 + (-t124 / 0.2e1 - t95 / 0.2e1) * t36 + (-Ifges(7,5) * t95 - Ifges(7,6) * t94) * t319 + (qJ(4) * t60 - t42 * t57 - t43 * t58 + (qJD(4) - t101) * t88) * m(6) + (-t123 / 0.2e1 - t94 / 0.2e1) * t35 + (mrSges(7,1) * t293 + mrSges(7,2) * t292) * t63 + t294 * qJD(4) + t281 * t118 - t280 * t119 + t346 * t44 + t206 * t54 / 0.2e1 + t187 * t14 - t156 * t154 - t58 * t116 - t57 * t117 - t101 * t96 + t103 * t21 + t104 * t22 - t76 * mrSges(4,2) - t68 * mrSges(5,3) + qJ(4) * t59 + t353 * t65 + t354 * t64 + (t103 * t3 + t104 * t2 + t12 * t353 + t13 * t354 + t187 * t41 + t346 * t63) * m(7) + t53 * t317 + t113 * t235 / 0.2e1 + t112 * t237 / 0.2e1 + t60 * t238 + (Ifges(7,5) * t124 + Ifges(7,6) * t123) * t320 + t10 * t321 + t9 * t322 + (Ifges(7,1) * t124 + Ifges(7,4) * t123) * t331 + (Ifges(7,4) * t124 + Ifges(7,2) * t123) * t333; t225 * t22 - t224 * t21 + t292 * t65 + t293 * t64 + t227 * qJD(5) + t261 * qJD(3) + (mrSges(5,1) * t271 + (t154 + t227) * t203) * qJD(2) + t229 + (-qJD(3) * t63 - t361) * m(7) + (-qJD(3) * t88 - t185 * t230 + t343) * m(6) + (qJD(3) * t108 + t120 * t194 + t77) * m(5); (t148 * t42 + t149 * t43) * mrSges(6,3) + (-Ifges(6,2) * t149 + t144 + t74) * t324 - m(7) * (t12 * t15 + t13 * t16) + t183 + t87 * t359 + (-m(7) * t12 * t266 + (m(7) * t2 - qJD(6) * t65 + t22) * t201 + (t21 + t64 * qJD(6) + m(7) * (t3 + t344)) * t205 + t345 * t149) * pkin(5) - t88 * (mrSges(6,1) * t149 + mrSges(6,2) * t148) - t42 * t116 + t43 * t117 - t16 * t64 - t15 * t65 + t339 + (Ifges(6,5) * t148 - Ifges(6,6) * t149) * t318 + (Ifges(6,1) * t148 - t311) * t323 + t149 * t334 + t368; -t12 * t64 + t13 * t65 + t35 * t330 + t363 + t368;];
tauc  = t1(:);
