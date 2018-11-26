% Calculate vector of centrifugal and coriolis load on the joints for
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:19
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:18:43
% EndTime: 2018-11-23 15:18:51
% DurationCPUTime: 7.98s
% Computational Cost: add. (5848->561), mult. (13806->784), div. (0->0), fcn. (9034->10), ass. (0->280)
t204 = sin(qJ(3));
t205 = sin(qJ(2));
t207 = cos(qJ(5));
t200 = sin(pkin(6));
t278 = qJD(1) * t200;
t203 = sin(qJ(5));
t209 = cos(qJ(2));
t287 = t203 * t209;
t115 = (t204 * t287 + t205 * t207) * t278;
t273 = qJD(3) * t204;
t194 = pkin(3) * t273;
t208 = cos(qJ(3));
t232 = pkin(9) * t204 - qJ(4) * t208;
t270 = qJD(4) * t204;
t214 = qJD(3) * t232 - t270;
t121 = t194 + t214;
t210 = -pkin(3) - pkin(9);
t247 = -qJ(4) * t204 - pkin(2);
t143 = t208 * t210 + t247;
t271 = qJD(3) * t208;
t329 = pkin(4) + pkin(8);
t158 = t329 * t271;
t171 = t329 * t204;
t268 = qJD(5) * t207;
t269 = qJD(5) * t203;
t37 = t207 * t121 - t143 * t269 + t203 * t158 + t171 * t268;
t374 = -t115 + t37;
t284 = t207 * t209;
t114 = (-t203 * t205 + t204 * t284) * t278;
t152 = t203 * t171;
t223 = -pkin(10) * t203 * t204 + pkin(5) * t208;
t244 = -t121 * t203 + t207 * t158;
t246 = pkin(10) * t208 - t143;
t373 = -t114 + t223 * qJD(3) + (t207 * t246 - t152) * qJD(5) + t244;
t267 = qJD(5) * t208;
t250 = t203 * t267;
t272 = qJD(3) * t207;
t215 = t204 * t272 + t250;
t372 = -pkin(10) * t215 - t374;
t312 = pkin(10) - t210;
t255 = t205 * t278;
t160 = qJD(2) * pkin(8) + t255;
t201 = cos(pkin(6));
t277 = qJD(1) * t204;
t177 = t201 * t277;
t119 = t208 * t160 + t177;
t274 = qJD(2) * t208;
t102 = pkin(4) * t274 + t119;
t195 = t204 * qJD(2);
t192 = pkin(3) * t195;
t125 = qJD(2) * t232 + t192;
t57 = t207 * t102 - t125 * t203;
t371 = -qJD(2) * t223 + t312 * t269 - t57;
t163 = t312 * t207;
t251 = t207 * t195;
t58 = t203 * t102 + t207 * t125;
t370 = pkin(10) * t251 + qJD(5) * t163 + t58;
t254 = t209 * t278;
t105 = qJD(2) * t143 - t254;
t245 = pkin(4) * qJD(2) + t160;
t228 = t245 * t204;
t289 = t201 * t208;
t178 = qJD(1) * t289;
t362 = qJD(4) - t178;
t80 = qJD(3) * t210 + t228 + t362;
t42 = -t105 * t203 + t207 * t80;
t43 = t105 * t207 + t203 * t80;
t230 = t203 * t42 - t207 * t43;
t311 = Ifges(6,4) * t203;
t234 = Ifges(6,2) * t207 + t311;
t310 = Ifges(6,4) * t207;
t236 = Ifges(6,1) * t203 + t310;
t239 = mrSges(6,1) * t207 - mrSges(6,2) * t203;
t306 = Ifges(6,6) * t207;
t309 = Ifges(6,5) * t203;
t316 = -t207 / 0.2e1;
t317 = -t203 / 0.2e1;
t186 = t195 + qJD(5);
t318 = -t186 / 0.2e1;
t149 = -t203 * t274 + t272;
t323 = -t149 / 0.2e1;
t148 = -qJD(3) * t203 - t207 * t274;
t324 = -t148 / 0.2e1;
t302 = t149 * Ifges(6,4);
t73 = t148 * Ifges(6,2) + t186 * Ifges(6,6) + t302;
t144 = Ifges(6,4) * t148;
t74 = Ifges(6,1) * t149 + Ifges(6,5) * t186 + t144;
t199 = qJD(3) * qJ(4);
t88 = t199 + t102;
t369 = t230 * mrSges(6,3) + (t306 + t309) * t318 + t234 * t324 + t236 * t323 + t88 * t239 + t316 * t73 + t317 * t74;
t206 = cos(qJ(6));
t32 = -pkin(10) * t149 + t42;
t24 = pkin(5) * t186 + t32;
t202 = sin(qJ(6));
t33 = pkin(10) * t148 + t43;
t296 = t202 * t33;
t12 = t206 * t24 - t296;
t295 = t206 * t33;
t13 = t202 * t24 + t295;
t265 = qJD(2) * qJD(3);
t248 = t208 * t265;
t183 = Ifges(7,3) * t248;
t243 = t206 * t148 - t149 * t202;
t87 = t148 * t202 + t149 * t206;
t315 = Ifges(7,4) * t87;
t263 = qJD(5) + qJD(6);
t175 = t195 + t263;
t320 = -t175 / 0.2e1;
t331 = -t87 / 0.2e1;
t333 = -t243 / 0.2e1;
t81 = Ifges(7,4) * t243;
t36 = Ifges(7,1) * t87 + Ifges(7,5) * t175 + t81;
t63 = -pkin(5) * t148 + t88;
t368 = t183 + (Ifges(7,5) * t243 - Ifges(7,6) * t87) * t320 + (t12 * t243 + t13 * t87) * mrSges(7,3) + (-Ifges(7,2) * t87 + t36 + t81) * t333 - t63 * (mrSges(7,1) * t87 + mrSges(7,2) * t243) + (Ifges(7,1) * t243 - t315) * t331;
t153 = t207 * t171;
t67 = pkin(5) * t204 + t203 * t246 + t153;
t285 = t207 * t208;
t92 = t207 * t143 + t152;
t75 = -pkin(10) * t285 + t92;
t31 = t202 * t67 + t206 * t75;
t367 = -qJD(6) * t31 + t372 * t202 + t206 * t373;
t30 = -t202 * t75 + t206 * t67;
t366 = qJD(6) * t30 + t202 * t373 - t372 * t206;
t355 = mrSges(5,1) + mrSges(4,3);
t365 = mrSges(5,2) - mrSges(4,1);
t108 = -t199 - t119;
t164 = -pkin(3) * t208 + t247;
t120 = qJD(2) * t164 - t254;
t305 = qJD(2) * pkin(2);
t161 = -t254 - t305;
t356 = qJD(3) / 0.2e1;
t357 = -qJD(3) / 0.2e1;
t358 = -qJD(2) / 0.2e1;
t364 = t119 * mrSges(4,3) + t120 * mrSges(5,2) + Ifges(4,6) * t356 + (Ifges(4,4) * t204 + t208 * Ifges(4,2)) * qJD(2) / 0.2e1 + Ifges(5,5) * t357 + (-Ifges(5,6) * t204 - t208 * Ifges(5,3)) * t358 - t108 * mrSges(5,1) - t161 * mrSges(4,1) + t369;
t264 = qJD(3) * qJD(5);
t113 = qJD(2) * t215 - t207 * t264;
t276 = qJD(2) * t200;
t249 = qJD(1) * t276;
t242 = t209 * t249;
t66 = t204 * t242 + (t208 * t245 + t177) * qJD(3);
t279 = qJD(3) * t192 + t205 * t249;
t78 = qJD(2) * t214 + t279;
t17 = -t105 * t269 + t203 * t66 + t207 * t78 + t80 * t268;
t11 = pkin(10) * t113 + t17;
t112 = -t203 * t264 + (t203 * t273 - t207 * t267) * qJD(2);
t18 = -qJD(5) * t43 - t203 * t78 + t207 * t66;
t8 = pkin(5) * t248 - pkin(10) * t112 + t18;
t2 = qJD(6) * t12 + t11 * t206 + t202 * t8;
t28 = qJD(6) * t243 + t112 * t206 + t113 * t202;
t29 = -qJD(6) * t87 - t112 * t202 + t113 * t206;
t344 = qJD(6) * t13;
t3 = -t11 * t202 + t206 * t8 - t344;
t363 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t28 + Ifges(7,6) * t29;
t286 = t206 * t207;
t288 = t202 * t203;
t224 = -t286 + t288;
t225 = t202 * t207 + t206 * t203;
t217 = t225 * t204;
t124 = qJD(2) * t217;
t95 = t263 * t225;
t292 = -t95 - t124;
t123 = -t195 * t288 + t206 * t251;
t266 = qJD(6) * t202;
t94 = -t202 * t269 - t203 * t266 + t263 * t286;
t293 = t94 + t123;
t361 = -t12 * t292 - t13 * t293 - t2 * t225 + t224 * t3;
t35 = Ifges(7,2) * t243 + Ifges(7,6) * t175 + t315;
t359 = t35 / 0.2e1;
t162 = t312 * t203;
t103 = t162 * t202 - t163 * t206;
t354 = qJD(6) * t103 + t202 * t371 - t370 * t206;
t104 = -t162 * t206 - t163 * t202;
t353 = -qJD(6) * t104 + t370 * t202 + t206 * t371;
t256 = -pkin(5) * t207 - pkin(4);
t346 = pkin(5) * t268 - (qJD(2) * t256 - t160) * t204 + t362;
t44 = -mrSges(7,1) * t243 + mrSges(7,2) * t87;
t345 = -m(7) * t63 - t44;
t343 = t17 * t203 + t18 * t207;
t118 = t160 * t204 - t178;
t342 = -qJD(4) - t118;
t340 = t208 * t263;
t339 = t18 * mrSges(6,1) - t17 * mrSges(6,2) + Ifges(6,5) * t112 + Ifges(6,6) * t113 + t363;
t337 = -m(5) / 0.2e1;
t336 = t28 / 0.2e1;
t335 = t29 / 0.2e1;
t334 = t73 / 0.2e1;
t332 = t243 / 0.2e1;
t330 = t87 / 0.2e1;
t126 = t224 * t208;
t326 = t126 / 0.2e1;
t127 = t225 * t208;
t325 = -t127 / 0.2e1;
t322 = -t225 / 0.2e1;
t321 = -t224 / 0.2e1;
t319 = t175 / 0.2e1;
t308 = Ifges(6,5) * t207;
t307 = Ifges(6,6) * t203;
t290 = t200 * t205;
t257 = t204 * t290;
t132 = t257 - t289;
t252 = t209 * t276;
t77 = t160 * t271 + (qJD(3) * t201 + t252) * t277;
t304 = t132 * t77;
t169 = -mrSges(5,1) * t274 - qJD(3) * mrSges(5,3);
t96 = -mrSges(6,1) * t148 + mrSges(6,2) * t149;
t294 = -t169 + t96;
t154 = (mrSges(5,2) * t208 - mrSges(5,3) * t204) * qJD(2);
t283 = t154 + (-mrSges(4,1) * t208 + mrSges(4,2) * t204) * qJD(2);
t282 = qJD(3) * t178 + t208 * t242;
t168 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t274;
t281 = t169 - t168;
t280 = qJD(3) * t365 + t195 * t355;
t172 = t329 * t208;
t275 = qJD(2) * t205;
t262 = pkin(8) * t77 / 0.2e1;
t261 = -t44 - t294;
t260 = -0.3e1 / 0.2e1 * Ifges(4,4) - 0.3e1 / 0.2e1 * Ifges(5,6);
t259 = -Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t258 = Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1;
t253 = t200 * t275;
t241 = -t254 / 0.2e1;
t238 = mrSges(6,1) * t203 + mrSges(6,2) * t207;
t237 = Ifges(6,1) * t207 - t311;
t235 = -Ifges(6,2) * t203 + t310;
t82 = mrSges(6,1) * t248 - mrSges(6,3) * t112;
t83 = -mrSges(6,2) * t248 + mrSges(6,3) * t113;
t229 = t203 * t83 + t207 * t82;
t221 = -t132 * t203 + t200 * t284;
t99 = t132 * t207 + t200 * t287;
t49 = t202 * t99 - t206 * t221;
t48 = t202 * t221 + t206 * t99;
t116 = -mrSges(6,2) * t186 + mrSges(6,3) * t148;
t117 = mrSges(6,1) * t186 - mrSges(6,3) * t149;
t227 = t207 * t116 - t203 * t117;
t76 = -t160 * t273 + t282;
t133 = t201 * t204 + t208 * t290;
t219 = -qJ(4) * t271 - t270;
t218 = -m(4) * t119 + m(5) * t108 + t281;
t198 = qJD(3) * qJD(4);
t60 = -qJD(3) * t228 + t198 + t282;
t107 = -qJD(3) * pkin(3) - t342;
t191 = Ifges(4,4) * t274;
t213 = t107 * mrSges(5,1) + t118 * mrSges(4,3) + t12 * mrSges(7,1) + t161 * mrSges(4,2) + t42 * mrSges(6,1) + Ifges(4,1) * t195 / 0.2e1 + Ifges(4,5) * t356 + t191 / 0.2e1 + Ifges(5,4) * t357 + (-t204 * Ifges(5,2) - Ifges(5,6) * t208) * t358 + t175 * Ifges(7,3) + t87 * Ifges(7,5) + t243 * Ifges(7,6) + t186 * Ifges(6,3) + t149 * Ifges(6,5) + t148 * Ifges(6,6) - t120 * mrSges(5,3) - t13 * mrSges(7,2) - t43 * mrSges(6,2);
t211 = qJD(2) ^ 2;
t188 = pkin(5) * t203 + qJ(4);
t184 = Ifges(6,3) * t248;
t157 = t329 * t273;
t156 = -qJ(4) * t274 + t192;
t141 = (mrSges(4,1) * t204 + mrSges(4,2) * t208) * t265;
t140 = (-mrSges(5,2) * t204 - mrSges(5,3) * t208) * t265;
t131 = pkin(5) * t285 + t172;
t129 = t194 + t219;
t106 = -pkin(5) * t250 + (-pkin(8) + t256) * t273;
t101 = t178 - t228;
t98 = qJD(3) * t133 + t204 * t252;
t97 = qJD(3) * t257 - t201 * t271 - t208 * t252;
t93 = qJD(2) * t219 + t279;
t91 = -t143 * t203 + t153;
t68 = -t198 - t76;
t65 = mrSges(7,1) * t175 - mrSges(7,3) * t87;
t64 = -mrSges(7,2) * t175 + mrSges(7,3) * t243;
t59 = -mrSges(6,1) * t113 + mrSges(6,2) * t112;
t54 = t112 * Ifges(6,1) + t113 * Ifges(6,4) + Ifges(6,5) * t248;
t53 = t112 * Ifges(6,4) + t113 * Ifges(6,2) + Ifges(6,6) * t248;
t52 = -t224 * t273 + t225 * t340;
t51 = qJD(3) * t217 + t224 * t340;
t41 = -pkin(5) * t113 + t60;
t40 = qJD(5) * t221 - t203 * t253 + t207 * t98;
t39 = qJD(5) * t99 + t203 * t98 + t207 * t253;
t38 = -qJD(5) * t92 + t244;
t22 = -mrSges(7,2) * t248 + mrSges(7,3) * t29;
t21 = mrSges(7,1) * t248 - mrSges(7,3) * t28;
t16 = t206 * t32 - t296;
t15 = -t202 * t32 - t295;
t14 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t10 = t28 * Ifges(7,1) + t29 * Ifges(7,4) + Ifges(7,5) * t248;
t9 = t28 * Ifges(7,4) + t29 * Ifges(7,2) + Ifges(7,6) * t248;
t7 = -qJD(6) * t49 - t202 * t39 + t206 * t40;
t6 = qJD(6) * t48 + t202 * t40 + t206 * t39;
t1 = [-t221 * t83 + t39 * t116 + t40 * t117 + t48 * t21 + t49 * t22 + t6 * t64 + t7 * t65 + t99 * t82 + t280 * t98 + (t59 + t14) * t133 + (-t168 + t261) * t97 + (-t211 * t205 * mrSges(3,1) + (-mrSges(3,2) * t211 - t140 - t141) * t209) * t200 + (t283 * t290 + t355 * qJD(3) * (t132 * t208 - t133 * t204)) * qJD(2) + m(5) * (t107 * t98 + t108 * t97 + t304 - t133 * t68 + (t120 * t275 - t209 * t93) * t200) + m(4) * (t118 * t98 - t119 * t97 + t304 + t133 * t76 + (t161 - t254) * t253) + m(7) * (t12 * t7 + t13 * t6 + t133 * t41 + t2 * t49 + t3 * t48 - t63 * t97) + m(6) * (t133 * t60 - t17 * t221 + t18 * t99 + t39 * t43 + t40 * t42 - t88 * t97); -m(6) * (t114 * t42 + t115 * t43) + t52 * t359 + 0.2e1 * (t120 * t337 + (-t305 / 0.2e1 - t161 / 0.2e1) * m(4)) * t255 - t283 * t255 + t172 * t59 - t157 * t96 + t164 * t140 + t129 * t154 + t131 * t14 - pkin(2) * t141 + t106 * t44 + t91 * t82 + t92 * t83 + t63 * (-mrSges(7,1) * t52 + mrSges(7,2) * t51) + t51 * t36 / 0.2e1 + t30 * t21 + t31 * t22 + t374 * t116 + (-t114 + t38) * t117 + (Ifges(7,5) * t51 + Ifges(7,6) * t52) * t319 + t10 * t325 + t9 * t326 + (Ifges(7,1) * t51 + Ifges(7,4) * t52) * t330 + (Ifges(7,4) * t51 + Ifges(7,2) * t52) * t332 + (-t93 * mrSges(5,3) + t184 / 0.2e1 + t183 / 0.2e1 + t355 * t77 + (mrSges(4,2) * t275 - t209 * t280) * t278 + 0.2e1 * (t107 * t241 + t262) * m(5) + 0.2e1 * (t118 * t241 + t262) * m(4) + (pkin(8) * t218 + t258 * qJD(3) + t195 * t260 - t364) * qJD(3) + t339) * t204 + (-t112 * t236 / 0.2e1 - t113 * t234 / 0.2e1 + t60 * t239 + t93 * mrSges(5,2) + t54 * t317 - t68 * mrSges(5,1) + t76 * mrSges(4,3) + t53 * t316 + (-t17 * t207 + t18 * t203) * mrSges(6,3) + (-mrSges(4,1) * t275 + (-m(6) * t88 + t218 + t345 - t96) * t209) * t278 + (t203 * t334 + t74 * t316 + t186 * (t307 - t308) / 0.2e1 + t235 * t324 + t237 * t323 - t88 * t238 + (t43 * t203 + t42 * t207) * mrSges(6,3)) * qJD(5) + (t259 * qJD(3) + (Ifges(7,5) * t325 + Ifges(7,6) * t326 + (-t309 / 0.2e1 - t306 / 0.2e1 - t260) * t208 + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t204) * qJD(2) + t213) * qJD(3) + (0.2e1 * t68 * t337 + m(4) * t76 + (m(4) * t118 + m(5) * t107 + t280) * qJD(3)) * pkin(8)) * t208 + m(5) * (t120 * t129 + t164 * t93) + (-t12 * t51 + t126 * t2 + t127 * t3 + t13 * t52) * mrSges(7,3) + t41 * (-mrSges(7,1) * t126 - mrSges(7,2) * t127) + (-Ifges(7,4) * t127 + Ifges(7,2) * t126) * t335 + (-Ifges(7,1) * t127 + Ifges(7,4) * t126) * t336 + t366 * t64 + (t106 * t63 + t12 * t367 + t13 * t366 + t131 * t41 + t2 * t31 + t3 * t30) * m(7) + t367 * t65 + m(6) * (-t157 * t88 + t17 * t92 + t172 * t60 + t18 * t91 + t37 * t43 + t38 * t42); t369 * qJD(5) - t343 * mrSges(6,3) + (-t123 / 0.2e1 - t94 / 0.2e1) * t35 + t41 * (mrSges(7,1) * t225 - mrSges(7,2) * t224) + (-Ifges(7,4) * t224 - Ifges(7,2) * t225) * t335 + (-Ifges(7,1) * t224 - Ifges(7,4) * t225) * t336 + (-t124 / 0.2e1 - t95 / 0.2e1) * t36 + (-Ifges(7,5) * t95 - Ifges(7,6) * t94) * t319 + (-Ifges(7,1) * t95 - Ifges(7,4) * t94) * t330 + (-Ifges(7,4) * t95 - Ifges(7,2) * t94) * t332 + (mrSges(7,1) * t293 + mrSges(7,2) * t292) * t63 + t294 * qJD(4) - t280 * t119 - t281 * t118 + t207 * t54 / 0.2e1 + t188 * t14 - t156 * t154 - t58 * t116 - t57 * t117 - t101 * t96 + t103 * t21 + t104 * t22 - t76 * mrSges(4,2) - t68 * mrSges(5,3) + qJ(4) * t59 + (-pkin(3) * t77 - qJ(4) * t68 - t107 * t119 + t108 * t342 - t120 * t156) * m(5) + t361 * mrSges(7,3) + ((-m(6) * t230 + t227) * qJD(5) + m(6) * t343 + t229) * t210 + t53 * t317 + (Ifges(7,5) * t124 + Ifges(7,6) * t123) * t320 + t10 * t321 + t9 * t322 + (Ifges(7,1) * t124 + Ifges(7,4) * t123) * t331 + (Ifges(7,4) * t124 + Ifges(7,2) * t123) * t333 + (((Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1) * t195 + (-qJ(4) * mrSges(5,1) + t258) * qJD(3) + t364) * t204 + ((Ifges(7,5) * t321 + Ifges(7,6) * t322 + t308 / 0.2e1 - t307 / 0.2e1 - pkin(3) * mrSges(5,1) + t259) * qJD(3) - t191 / 0.2e1 + (-Ifges(5,2) / 0.2e1 + Ifges(5,3) / 0.2e1 - Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t195 - Ifges(5,6) * t274 / 0.2e1 - t213) * t208) * qJD(2) + t346 * t44 + (qJ(4) * t60 - t42 * t57 - t43 * t58 + (qJD(4) - t101) * t88) * m(6) + t353 * t65 + t354 * t64 + (t103 * t3 + t104 * t2 + t12 * t353 + t13 * t354 + t188 * t41 + t346 * t63) * m(7) + t365 * t77 + t113 * t235 / 0.2e1 + t112 * t237 / 0.2e1 + t60 * t238; t225 * t22 - t224 * t21 + t292 * t65 + t293 * t64 + t227 * qJD(5) + t261 * qJD(3) + (mrSges(5,1) * t271 + (t154 + t227) * t204) * qJD(2) + t229 + (-qJD(3) * t63 - t361) * m(7) + (-qJD(3) * t88 - t186 * t230 + t343) * m(6) + (qJD(3) * t108 + t120 * t195 + t77) * m(5); (t148 * t42 + t149 * t43) * mrSges(6,3) + t184 - m(7) * (t12 * t15 + t13 * t16) + (-Ifges(6,2) * t149 + t144 + t74) * t324 + t339 + t87 * t359 - t88 * (mrSges(6,1) * t149 + mrSges(6,2) * t148) - t42 * t116 + t43 * t117 - t15 * t65 - t16 * t64 + (Ifges(6,5) * t148 - Ifges(6,6) * t149) * t318 + (Ifges(6,1) * t148 - t302) * t323 + t149 * t334 + (-m(7) * t12 * t266 + (m(7) * t2 - qJD(6) * t65 + t22) * t202 + (t21 + t64 * qJD(6) + m(7) * (t3 + t344)) * t206 + t345 * t149) * pkin(5) + t368; -t12 * t64 + t13 * t65 + t35 * t330 + t363 + t368;];
tauc  = t1(:);
