% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2018-11-23 15:11
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:11:11
% EndTime: 2018-11-23 15:11:21
% DurationCPUTime: 9.97s
% Computational Cost: add. (5359->516), mult. (14184->686), div. (0->0), fcn. (10313->10), ass. (0->254)
t360 = Ifges(6,4) + Ifges(7,4);
t180 = sin(qJ(2));
t176 = sin(pkin(6));
t262 = qJD(1) * t176;
t240 = t180 * t262;
t179 = sin(qJ(3));
t257 = qJD(3) * t179;
t365 = pkin(3) * t257 - t240;
t361 = Ifges(6,1) + Ifges(7,1);
t347 = Ifges(7,5) + Ifges(6,5);
t359 = Ifges(6,2) + Ifges(7,2);
t358 = Ifges(7,6) + Ifges(6,6);
t175 = sin(pkin(11));
t182 = cos(qJ(3));
t273 = cos(pkin(11));
t225 = t273 * t179;
t153 = t175 * t182 + t225;
t143 = t153 * qJD(3);
t194 = -t175 * t179 + t182 * t273;
t144 = t194 * qJD(3);
t364 = pkin(4) * t143 - pkin(9) * t144 + t365;
t295 = -qJ(4) - pkin(8);
t227 = qJD(3) * t295;
t255 = qJD(4) * t182;
t138 = t179 * t227 + t255;
t139 = -qJD(4) * t179 + t182 * t227;
t183 = cos(qJ(2));
t239 = t183 * t262;
t344 = t138 * t273 + t175 * t139 - t194 * t239;
t181 = cos(qJ(5));
t363 = t360 * t181;
t178 = sin(qJ(5));
t362 = t360 * t178;
t141 = t194 * qJD(2);
t305 = -t141 / 0.2e1;
t357 = Ifges(5,2) * t305;
t356 = -t344 * t178 + t364 * t181;
t173 = -pkin(3) * t182 - pkin(2);
t103 = -pkin(4) * t194 - pkin(9) * t153 + t173;
t253 = qJD(5) * t181;
t355 = t103 * t253 + t364 * t178 + t344 * t181;
t258 = qJD(2) * t182;
t142 = -qJD(2) * t225 - t175 * t258;
t120 = qJD(3) * t181 + t142 * t178;
t354 = t360 * t120;
t327 = -t358 * t178 + t347 * t181;
t326 = -t359 * t178 + t363;
t325 = t361 * t181 - t362;
t254 = qJD(5) * t178;
t269 = t141 * t178;
t353 = t254 - t269;
t121 = qJD(3) * t178 - t142 * t181;
t352 = t360 * t121;
t162 = t295 * t179;
t163 = t295 * t182;
t117 = t175 * t162 - t163 * t273;
t109 = t181 * t117;
t199 = -qJ(6) * t144 - qJD(6) * t153;
t270 = qJ(6) * t153;
t351 = pkin(5) * t143 + t199 * t181 + (-t109 + (-t103 + t270) * t178) * qJD(5) + t356;
t235 = t153 * t253;
t350 = -qJ(6) * t235 + (-qJD(5) * t117 + t199) * t178 + t355;
t349 = -t117 * t254 + t355;
t48 = t178 * t103 + t109;
t348 = -qJD(5) * t48 + t356;
t340 = qJD(5) - t141;
t336 = t120 * t359 + t340 * t358 + t352;
t335 = t121 * t361 + t347 * t340 + t354;
t114 = t153 * t239;
t85 = t138 * t175 - t273 * t139;
t345 = t114 - t85;
t343 = t178 * t347 + t181 * t358;
t342 = t181 * t359 + t362;
t341 = t178 * t361 + t363;
t339 = -Ifges(5,6) / 0.2e1;
t133 = qJD(2) * t143;
t134 = qJD(2) * t144;
t73 = qJD(5) * t120 + t134 * t181;
t74 = -qJD(5) * t121 - t134 * t178;
t338 = t358 * t133 + t359 * t74 + t360 * t73;
t337 = t347 * t133 + t360 * t74 + t361 * t73;
t135 = qJD(2) * t173 + qJD(4) - t239;
t334 = t135 * mrSges(5,2);
t299 = pkin(3) * t175;
t171 = pkin(9) + t299;
t263 = qJ(6) + t171;
t224 = qJD(5) * t263;
t157 = qJD(2) * pkin(8) + t240;
t223 = qJ(4) * qJD(2) + t157;
t177 = cos(pkin(6));
t261 = qJD(1) * t179;
t238 = t177 * t261;
t113 = t182 * t223 + t238;
t104 = t175 * t113;
t264 = t177 * t182;
t168 = qJD(1) * t264;
t112 = -t179 * t223 + t168;
t55 = t112 * t273 - t104;
t260 = qJD(2) * t179;
t249 = pkin(3) * t260;
t88 = -pkin(4) * t142 - pkin(9) * t141 + t249;
t24 = t178 * t88 + t181 * t55;
t333 = qJ(6) * t269 + qJD(6) * t181 - t178 * t224 - t24;
t23 = -t178 * t55 + t181 * t88;
t268 = t141 * t181;
t332 = pkin(5) * t142 + qJ(6) * t268 - qJD(6) * t178 - t181 * t224 - t23;
t289 = mrSges(5,3) * t142;
t274 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t120 + mrSges(6,2) * t121 - t289;
t226 = t273 * t113;
t53 = t112 * t175 + t226;
t331 = t353 * pkin(5) - t53;
t330 = Ifges(5,5) * qJD(3);
t252 = qJD(2) * qJD(3);
t232 = t179 * t252;
t123 = t157 * t182 + t238;
t256 = qJD(3) * t182;
t329 = -t123 * qJD(3) + (-qJ(4) * t256 + (-qJD(4) - t239) * t179) * qJD(2);
t215 = mrSges(7,1) * t178 + mrSges(7,2) * t181;
t217 = mrSges(6,1) * t178 + mrSges(6,2) * t181;
t108 = qJD(3) * pkin(3) + t112;
t50 = t108 * t273 - t104;
t45 = -qJD(3) * pkin(4) - t50;
t31 = -t120 * pkin(5) + qJD(6) + t45;
t328 = -t31 * t215 - t45 * t217;
t324 = -t253 + t268;
t265 = t176 * t183;
t236 = qJD(2) * t265;
t222 = t182 * t236;
t99 = qJD(1) * t222 + qJD(3) * t168 - t157 * t257;
t72 = (-qJ(4) * t257 + t255) * qJD(2) + t99;
t27 = t175 * t329 + t273 * t72;
t51 = t175 * t108 + t226;
t46 = qJD(3) * pkin(9) + t51;
t259 = qJD(2) * t180;
t237 = t176 * t259;
t140 = pkin(3) * t232 + qJD(1) * t237;
t59 = pkin(4) * t133 - pkin(9) * t134 + t140;
t64 = -pkin(4) * t141 + pkin(9) * t142 + t135;
t3 = t178 * t59 + t181 * t27 + t64 * t253 - t254 * t46;
t18 = t178 * t64 + t181 * t46;
t4 = -qJD(5) * t18 - t178 * t27 + t181 * t59;
t322 = -t178 * t4 + t181 * t3;
t320 = qJD(3) * t339;
t287 = Ifges(5,4) * t142;
t319 = t357 + t287 / 0.2e1;
t309 = t121 / 0.2e1;
t318 = -t328 + t326 * t120 / 0.2e1 + t325 * t309 + t327 * t340 / 0.2e1;
t12 = qJ(6) * t120 + t18;
t17 = -t178 * t46 + t181 * t64;
t11 = -qJ(6) * t121 + t17;
t7 = pkin(5) * t340 + t11;
t317 = t135 * mrSges(5,1) + t17 * mrSges(6,1) + t7 * mrSges(7,1) - t18 * mrSges(6,2) - t12 * mrSges(7,2) + t319 + t320;
t184 = qJD(2) ^ 2;
t315 = t73 / 0.2e1;
t314 = t74 / 0.2e1;
t312 = -t120 / 0.2e1;
t310 = -t121 / 0.2e1;
t308 = t133 / 0.2e1;
t307 = -t340 / 0.2e1;
t304 = t142 / 0.2e1;
t303 = -t178 / 0.2e1;
t300 = t181 / 0.2e1;
t26 = t175 * t72 - t273 * t329;
t266 = t176 * t180;
t145 = -t179 * t266 + t264;
t146 = t177 * t179 + t182 * t266;
t93 = -t145 * t273 + t146 * t175;
t296 = t26 * t93;
t41 = mrSges(7,1) * t133 - mrSges(7,3) * t73;
t42 = mrSges(6,1) * t133 - mrSges(6,3) * t73;
t294 = t41 + t42;
t43 = -mrSges(7,2) * t133 + mrSges(7,3) * t74;
t44 = -mrSges(6,2) * t133 + mrSges(6,3) * t74;
t293 = t43 + t44;
t75 = -mrSges(7,2) * t340 + mrSges(7,3) * t120;
t76 = -mrSges(6,2) * t340 + mrSges(6,3) * t120;
t292 = t75 + t76;
t77 = mrSges(7,1) * t340 - mrSges(7,3) * t121;
t78 = mrSges(6,1) * t340 - mrSges(6,3) * t121;
t291 = t77 + t78;
t290 = mrSges(5,3) * t133;
t288 = Ifges(4,4) * t179;
t136 = Ifges(5,4) * t141;
t116 = -t273 * t162 - t163 * t175;
t281 = t116 * t26;
t278 = t134 * mrSges(5,3);
t277 = t141 * mrSges(5,3);
t275 = t142 * Ifges(5,1);
t272 = Ifges(4,5) * qJD(3);
t271 = Ifges(4,6) * qJD(3);
t267 = t153 * t178;
t62 = -mrSges(7,1) * t120 + mrSges(7,2) * t121;
t250 = -t62 - t274;
t246 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t245 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t244 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t243 = mrSges(4,3) * t260;
t242 = mrSges(4,3) * t258;
t241 = t273 * pkin(3);
t28 = -t74 * mrSges(7,1) + t73 * mrSges(7,2);
t87 = t133 * mrSges(5,1) + t134 * mrSges(5,2);
t47 = t181 * t103 - t117 * t178;
t172 = -t241 - pkin(4);
t1 = pkin(5) * t133 - qJ(6) * t73 - qJD(6) * t121 + t4;
t2 = qJ(6) * t74 + qJD(6) * t120 + t3;
t221 = -t1 * t181 - t2 * t178;
t220 = -t3 * t178 - t4 * t181;
t219 = t12 * t181 - t7 * t178;
t218 = mrSges(6,1) * t181 - mrSges(6,2) * t178;
t216 = mrSges(7,1) * t181 - mrSges(7,2) * t178;
t202 = -t17 * t181 - t18 * t178;
t201 = t17 * t178 - t18 * t181;
t100 = -t157 * t256 + (-qJD(3) * t177 - t236) * t261;
t200 = -t100 * t179 + t99 * t182;
t94 = t175 * t145 + t146 * t273;
t65 = -t178 * t94 - t181 * t265;
t198 = t178 * t265 - t181 * t94;
t195 = (mrSges(4,1) * t179 + mrSges(4,2) * t182) * qJD(2);
t186 = t4 * mrSges(6,1) + t1 * mrSges(7,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t174 = Ifges(4,4) * t258;
t161 = -qJD(3) * mrSges(4,2) + t242;
t160 = qJD(3) * mrSges(4,1) - t243;
t159 = -t181 * pkin(5) + t172;
t158 = -qJD(2) * pkin(2) - t239;
t151 = t263 * t181;
t150 = t263 * t178;
t149 = qJD(3) * t195;
t148 = Ifges(4,1) * t260 + t174 + t272;
t147 = t271 + (t182 * Ifges(4,2) + t288) * qJD(2);
t131 = Ifges(6,3) * t133;
t130 = Ifges(7,3) * t133;
t126 = -qJD(3) * mrSges(5,2) + t277;
t122 = -t157 * t179 + t168;
t111 = -qJD(3) * t146 - t179 * t236;
t110 = qJD(3) * t145 + t222;
t101 = -mrSges(5,1) * t141 - mrSges(5,2) * t142;
t96 = t136 - t275 + t330;
t84 = pkin(5) * t267 + t116;
t71 = Ifges(6,5) * t73;
t70 = Ifges(7,5) * t73;
t69 = Ifges(6,6) * t74;
t68 = Ifges(7,6) * t74;
t54 = t110 * t273 + t175 * t111;
t52 = t110 * t175 - t111 * t273;
t36 = t121 * Ifges(6,5) + t120 * Ifges(6,6) + Ifges(6,3) * t340;
t35 = t121 * Ifges(7,5) + t120 * Ifges(7,6) + Ifges(7,3) * t340;
t34 = (t178 * t144 + t235) * pkin(5) + t85;
t32 = -qJ(6) * t267 + t48;
t30 = -pkin(5) * t194 - t181 * t270 + t47;
t29 = -mrSges(6,1) * t74 + mrSges(6,2) * t73;
t16 = qJD(5) * t198 - t178 * t54 + t181 * t237;
t15 = qJD(5) * t65 + t178 * t237 + t181 * t54;
t10 = -pkin(5) * t74 + t26;
t5 = [-t94 * t290 + t110 * t161 + t111 * t160 + t54 * t126 - t293 * t198 + t294 * t65 + t291 * t16 + t292 * t15 + (-t145 * t182 - t146 * t179) * mrSges(4,3) * t252 + (t28 + t29 + t278) * t93 - t250 * t52 + ((-mrSges(3,2) * t184 - t149 - t87) * t183 + (-mrSges(3,1) * t184 + (t101 + qJD(2) * (-mrSges(4,1) * t182 + mrSges(4,2) * t179)) * qJD(2)) * t180) * t176 + m(5) * (t296 + t27 * t94 - t50 * t52 + t51 * t54 + (t135 * t259 - t140 * t183) * t176) + m(4) * (t100 * t145 + t110 * t123 + t111 * t122 + t146 * t99 + (t158 - t239) * t237) + m(7) * (t1 * t65 + t10 * t93 + t12 * t15 + t16 * t7 - t198 * t2 + t31 * t52) + m(6) * (t15 * t18 + t16 * t17 - t198 * t3 + t4 * t65 + t45 * t52 + t296); (((-t122 * t182 - t123 * t179) * qJD(3) + t200) * pkin(8) - (pkin(2) * t259 + t158 * t180 + (-t122 * t179 + t123 * t182) * t183) * t262) * m(4) + ((t179 * t160 - t182 * t161) * t183 - t180 * t101) * t262 + (-0.3e1 / 0.2e1 * t179 ^ 2 + 0.3e1 / 0.2e1 * t182 ^ 2) * Ifges(4,4) * t252 + t274 * t85 + (t35 / 0.2e1 + t36 / 0.2e1 - t51 * mrSges(5,3) + t244 * t340 + t246 * t121 + t245 * t120 + t317 + t319) * t143 + t250 * t114 + t344 * t126 + (t116 * t134 - t117 * t133) * mrSges(5,3) + t200 * mrSges(4,3) - (t140 * mrSges(5,1) - Ifges(5,4) * t134 + t70 / 0.2e1 + t68 / 0.2e1 + t130 / 0.2e1 + t71 / 0.2e1 + t69 / 0.2e1 + t131 / 0.2e1 - t27 * mrSges(5,3) + t245 * t74 + t246 * t73 + (Ifges(5,2) + t244) * t133 + t186) * t194 + t173 * t87 - pkin(2) * t149 + (t117 * t27 + t365 * t135 + t140 * t173 + t344 * t51 + t345 * t50 + t281) * m(5) + (Ifges(5,1) * t134 - Ifges(5,4) * t133 + t140 * mrSges(5,2) + t10 * t215 + (mrSges(5,3) + t217) * t26 + t221 * mrSges(7,3) + t220 * mrSges(6,3) + (mrSges(6,3) * t201 - mrSges(7,3) * t219 + t216 * t31 + t218 * t45 + t342 * t312 + t341 * t310 + t343 * t307 - t336 * t181 / 0.2e1) * qJD(5) + t325 * t315 + t326 * t314 + t327 * t308 + t337 * t300 + (qJD(5) * t335 + t338) * t303) * t153 + t348 * t78 + t349 * t76 + (t17 * t348 + t18 * t349 + t3 * t48 - t345 * t45 + t4 * t47 + t281) * m(6) + t350 * t75 + t351 * t77 + (t1 * t30 + t10 * t84 + t2 * t32 + t351 * t7 + (-t114 + t34) * t31 + t350 * t12) * m(7) + (t335 * t300 + t336 * t303 + t334 - t50 * mrSges(5,3) - t275 / 0.2e1 + t136 / 0.2e1 + t96 / 0.2e1 + t318 + (-t12 * t178 - t7 * t181) * mrSges(7,3) + t202 * mrSges(6,3)) * t144 + (Ifges(5,5) * t144 / 0.2e1 + t143 * t339 + (t148 / 0.2e1 + t158 * mrSges(4,2) - t122 * mrSges(4,3) - pkin(8) * t160 + t272 / 0.2e1) * t182 + (-t147 / 0.2e1 + t158 * mrSges(4,1) + pkin(3) * t101 - t123 * mrSges(4,3) - pkin(8) * t161 - t271 / 0.2e1 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t258) * t179) * qJD(3) + t30 * t41 + t32 * t43 + t47 * t42 + t48 * t44 + t34 * t62 + t84 * t28 + t116 * t29; (m(6) * t172 - mrSges(5,1) - t218) * t26 + t50 * t277 + (t253 / 0.2e1 - t268 / 0.2e1) * t335 - t290 * t299 - t51 * t289 - t241 * t278 + (t17 * t324 - t18 * t353 + t322) * mrSges(6,3) + (-t1 * t178 - t12 * t353 + t181 * t2 + t324 * t7) * mrSges(7,3) + (t357 + t320 - t358 * t312 - t347 * t310 + (-Ifges(7,3) - Ifges(6,3)) * t307 + t317) * t142 + (-t135 * t249 + t50 * t53 - t51 * t55 + (t175 * t27 - t26 * t273) * pkin(3)) * m(5) + t252 * Ifges(4,5) * t182 / 0.2e1 + t147 * t260 / 0.2e1 + (t96 + t136) * t305 + (-t254 / 0.2e1 + t269 / 0.2e1) * t336 - t101 * t249 - Ifges(4,6) * t232 / 0.2e1 + (t36 + t35 + t287) * t304 + (m(6) * t171 * t202 + t318) * qJD(5) - m(6) * (t17 * t23 + t18 * t24) - t10 * t216 - (-Ifges(4,2) * t260 + t148 + t174) * t258 / 0.2e1 - t158 * t195 + t341 * t315 + t342 * t314 + t343 * t308 + t172 * t29 + t159 * t28 - t150 * t41 + t151 * t43 + Ifges(5,5) * t134 + (t243 + t160) * t123 - Ifges(5,6) * t133 - t55 * t126 + (m(6) * t322 - t178 * t42 + t181 * t44 - t78 * t253 - t76 * t254) * t171 - t179 * t184 * (Ifges(4,1) * t182 - t288) / 0.2e1 + t331 * t62 + (-m(6) * t45 - t274) * t53 + t332 * t77 + t333 * t75 + (-t1 * t150 + t10 * t159 + t12 * t333 + t151 * t2 + t31 * t331 + t332 * t7) * m(7) + (Ifges(5,1) * t304 - t330 / 0.2e1 - t334 + t326 * t312 + t325 * t310 + t327 * t307 + t328) * t141 + t337 * t178 / 0.2e1 + t338 * t300 - t27 * mrSges(5,2) + (t242 - t161) * t122 - t24 * t76 - t23 * t78 - t99 * mrSges(4,2) + t100 * mrSges(4,1); -t141 * t126 - t250 * t142 + (t292 * t340 + t294) * t181 + (-t291 * t340 + t293) * t178 + t87 + (t142 * t31 + t219 * t340 - t221) * m(7) + (t142 * t45 - t201 * t340 - t220) * m(6) + (-t141 * t51 - t142 * t50 + t140) * m(5); (-(t11 - t7) * t12 + (-t121 * t31 + t1) * pkin(5)) * m(7) + (t12 * t121 + t120 * t7) * mrSges(7,3) + (t120 * t17 + t121 * t18) * mrSges(6,3) + t130 + t131 + t71 + t70 + t69 + t68 + (-t121 * t62 + t41) * pkin(5) + t186 - t11 * t75 - t17 * t76 + t12 * t77 + t18 * t78 - t31 * (mrSges(7,1) * t121 + mrSges(7,2) * t120) - t45 * (mrSges(6,1) * t121 + mrSges(6,2) * t120) + (t361 * t120 - t352) * t310 + t336 * t309 + (t120 * t347 - t121 * t358) * t307 + (-t359 * t121 + t335 + t354) * t312; -t120 * t75 + t121 * t77 + 0.2e1 * (t10 / 0.2e1 + t12 * t312 + t7 * t309) * m(7) + t28;];
tauc  = t5(:);
