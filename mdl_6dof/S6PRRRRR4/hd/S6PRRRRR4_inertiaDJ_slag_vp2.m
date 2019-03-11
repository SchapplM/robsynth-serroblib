% Calculate time derivative of joint inertia matrix for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:54:56
% EndTime: 2019-03-09 00:55:11
% DurationCPUTime: 6.87s
% Computational Cost: add. (10732->623), mult. (29527->928), div. (0->0), fcn. (29888->14), ass. (0->262)
t232 = sin(qJ(6));
t237 = cos(qJ(6));
t238 = cos(qJ(5));
t356 = (t232 ^ 2 + t237 ^ 2) * t238;
t205 = -mrSges(7,1) * t237 + mrSges(7,2) * t232;
t360 = t205 - mrSges(6,1);
t233 = sin(qJ(5));
t234 = sin(qJ(4));
t239 = cos(qJ(4));
t193 = t233 * t234 - t238 * t239;
t355 = qJD(4) + qJD(5);
t151 = t355 * t193;
t194 = t233 * t239 + t234 * t238;
t281 = qJD(6) * t237;
t251 = -t232 * t151 + t194 * t281;
t282 = qJD(6) * t232;
t293 = t237 * t151;
t250 = t194 * t282 + t293;
t343 = -pkin(11) - pkin(10);
t213 = t343 * t239;
t277 = t343 * t234;
t165 = -t213 * t233 - t238 * t277;
t274 = qJD(4) * t343;
t203 = t234 * t274;
t265 = t239 * t274;
t106 = -qJD(5) * t165 + t203 * t238 + t233 * t265;
t152 = t355 * t194;
t285 = qJD(4) * t234;
t278 = pkin(4) * t285;
t93 = pkin(5) * t152 + pkin(12) * t151 + t278;
t223 = -pkin(4) * t239 - pkin(3);
t140 = pkin(5) * t193 - pkin(12) * t194 + t223;
t166 = -t238 * t213 + t233 * t277;
t96 = t140 * t237 - t166 * t232;
t32 = qJD(6) * t96 + t106 * t237 + t232 * t93;
t97 = t140 * t232 + t166 * t237;
t33 = -qJD(6) * t97 - t106 * t232 + t237 * t93;
t359 = -t33 * t232 + t237 * t32;
t206 = -mrSges(5,1) * t239 + mrSges(5,2) * t234;
t358 = -m(5) * pkin(3) - mrSges(4,1) + t206;
t335 = t232 / 0.2e1;
t333 = t237 / 0.2e1;
t267 = -t282 / 0.2e1;
t230 = cos(pkin(7));
t235 = sin(qJ(3));
t228 = sin(pkin(7));
t240 = cos(qJ(3));
t300 = t228 * t240;
t189 = pkin(2) * t230 * t235 + pkin(9) * t300;
t175 = pkin(10) * t230 + t189;
t176 = (-pkin(3) * t240 - pkin(10) * t235 - pkin(2)) * t228;
t125 = t175 * t239 + t176 * t234;
t301 = t228 * t235;
t217 = pkin(9) * t301;
t330 = pkin(2) * t240;
t188 = t230 * t330 - t217;
t241 = cos(qJ(2));
t291 = t240 * t241;
t236 = sin(qJ(2));
t296 = t235 * t236;
t357 = t230 * t291 - t296;
t183 = t230 * t239 - t234 * t301;
t287 = qJD(3) * t228;
t271 = t240 * t287;
t158 = qJD(4) * t183 + t239 * t271;
t184 = t230 * t234 + t239 * t301;
t159 = -qJD(4) * t184 - t234 * t271;
t286 = qJD(3) * t235;
t272 = t228 * t286;
t354 = -Ifges(5,5) * t158 - Ifges(5,6) * t159 - Ifges(5,3) * t272;
t103 = pkin(11) * t183 + t125;
t124 = -t175 * t234 + t176 * t239;
t90 = -pkin(4) * t300 - pkin(11) * t184 + t124;
t324 = t103 * t238 + t233 * t90;
t179 = (pkin(3) * t235 - pkin(10) * t240) * t287;
t180 = t188 * qJD(3);
t81 = -qJD(4) * t125 + t179 * t239 - t180 * t234;
t61 = pkin(4) * t272 - pkin(11) * t158 + t81;
t284 = qJD(4) * t239;
t80 = -t175 * t285 + t176 * t284 + t179 * t234 + t180 * t239;
t64 = pkin(11) * t159 + t80;
t17 = -qJD(5) * t324 - t233 * t64 + t238 * t61;
t353 = 2 * m(6);
t352 = 2 * m(7);
t351 = 0.2e1 * pkin(4);
t350 = -2 * mrSges(4,3);
t349 = -2 * mrSges(6,3);
t107 = qJD(5) * t166 + t233 * t203 - t238 * t265;
t348 = 0.2e1 * t107;
t347 = 0.2e1 * t165;
t346 = m(6) / 0.2e1;
t134 = t183 * t233 + t184 * t238;
t253 = -t134 * t237 + t232 * t300;
t254 = t183 * t238 - t184 * t233;
t75 = qJD(5) * t254 + t158 * t238 + t159 * t233;
t39 = qJD(6) * t253 - t232 * t75 + t237 * t272;
t344 = t39 / 0.2e1;
t116 = -t134 * t232 - t237 * t300;
t342 = t116 / 0.2e1;
t224 = Ifges(7,5) * t281;
t341 = Ifges(7,6) * t267 + t224 / 0.2e1;
t320 = Ifges(7,4) * t232;
t261 = Ifges(7,1) * t237 - t320;
t201 = t261 * qJD(6);
t340 = t201 / 0.2e1;
t339 = Ifges(7,5) * t335 + Ifges(7,6) * t333;
t319 = Ifges(7,4) * t237;
t210 = Ifges(7,1) * t232 + t319;
t337 = t210 / 0.2e1;
t336 = -t232 / 0.2e1;
t334 = t234 / 0.2e1;
t332 = t239 / 0.2e1;
t229 = sin(pkin(6));
t231 = cos(pkin(6));
t109 = t231 * t271 + (t357 * qJD(3) + (-t230 * t296 + t291) * qJD(2)) * t229;
t294 = t236 * t240;
t295 = t235 * t241;
t252 = t230 * t295 + t294;
t148 = t229 * t252 + t231 * t301;
t182 = -t228 * t229 * t241 + t230 * t231;
t114 = t148 * t239 + t182 * t234;
t288 = qJD(2) * t229;
t273 = t236 * t288;
t264 = t228 * t273;
t59 = -qJD(4) * t114 - t109 * t234 + t239 * t264;
t113 = -t148 * t234 + t182 * t239;
t60 = qJD(4) * t113 + t109 * t239 + t234 * t264;
t70 = t113 * t233 + t114 * t238;
t21 = qJD(5) * t70 + t233 * t60 - t238 * t59;
t255 = t113 * t238 - t114 * t233;
t329 = t21 * t255;
t108 = t231 * t272 + (t252 * qJD(3) + (t230 * t294 + t295) * qJD(2)) * t229;
t20 = qJD(5) * t255 + t233 * t59 + t238 * t60;
t147 = -t229 * t357 - t231 * t300;
t43 = t147 * t237 - t232 * t70;
t5 = qJD(6) * t43 + t108 * t232 + t20 * t237;
t328 = t237 * t5;
t283 = qJD(5) * t238;
t309 = t103 * t233;
t16 = -qJD(5) * t309 + t233 * t61 + t238 * t64 + t283 * t90;
t13 = pkin(12) * t272 + t16;
t46 = -pkin(12) * t300 + t324;
t174 = t217 + (-pkin(3) - t330) * t230;
t136 = -pkin(4) * t183 + t174;
t65 = -pkin(5) * t254 - pkin(12) * t134 + t136;
t26 = t232 * t65 + t237 * t46;
t181 = t189 * qJD(3);
t123 = -pkin(4) * t159 + t181;
t76 = qJD(5) * t134 + t158 * t233 - t159 * t238;
t27 = pkin(5) * t76 - pkin(12) * t75 + t123;
t3 = -qJD(6) * t26 - t13 * t232 + t237 * t27;
t327 = t3 * t232;
t44 = t147 * t232 + t237 * t70;
t6 = -qJD(6) * t44 + t108 * t237 - t20 * t232;
t326 = t6 * t232;
t38 = qJD(6) * t116 + t232 * t272 + t237 * t75;
t22 = -mrSges(7,1) * t39 + mrSges(7,2) * t38;
t67 = mrSges(6,1) * t272 - mrSges(6,3) * t75;
t325 = t22 - t67;
t323 = mrSges(7,3) * t237;
t322 = Ifges(5,4) * t234;
t321 = Ifges(5,4) * t239;
t318 = Ifges(7,6) * t232;
t317 = pkin(4) * qJD(5);
t316 = t233 * mrSges(6,1);
t315 = t233 * t255;
t313 = t238 * mrSges(6,2);
t122 = -mrSges(6,1) * t300 - mrSges(6,3) * t134;
t71 = -mrSges(7,1) * t116 - mrSges(7,2) * t253;
t311 = -t122 + t71;
t308 = t107 * t165;
t82 = t147 * t108;
t307 = t147 * t181;
t306 = t165 * t233;
t305 = t194 * t232;
t304 = t194 * t237;
t298 = t232 * t238;
t297 = t233 * t205;
t292 = t237 * t238;
t290 = -Ifges(6,5) * t151 - Ifges(6,6) * t152;
t289 = -mrSges(4,1) * t230 - mrSges(5,1) * t183 + mrSges(5,2) * t184 + mrSges(4,3) * t301;
t9 = Ifges(7,5) * t38 + Ifges(7,6) * t39 + Ifges(7,3) * t76;
t279 = -Ifges(6,5) * t75 + Ifges(6,6) * t76 - Ifges(6,3) * t272;
t14 = -pkin(5) * t272 - t17;
t275 = m(7) * t14 + t22;
t77 = mrSges(7,1) * t251 - mrSges(7,2) * t250;
t268 = m(7) * t107 + t77;
t266 = t281 / 0.2e1;
t263 = mrSges(7,3) * t356;
t262 = mrSges(7,1) * t232 + mrSges(7,2) * t237;
t260 = -Ifges(7,2) * t232 + t319;
t259 = -t107 * t255 + t165 * t21;
t25 = -t232 * t46 + t237 * t65;
t257 = -t234 * t81 + t239 * t80;
t52 = t238 * t90 - t309;
t199 = t260 * qJD(6);
t208 = Ifges(7,2) * t237 + t320;
t249 = t199 * t237 + t201 * t232 - t208 * t282 + t210 * t281;
t248 = -t326 + (-t232 * t44 - t237 * t43) * qJD(6);
t54 = -Ifges(7,5) * t250 - Ifges(7,6) * t251 + Ifges(7,3) * t152;
t196 = t262 * qJD(6);
t246 = -t20 * mrSges(6,2) + mrSges(7,3) * t248 - t196 * t255 + t21 * t360 + t5 * t323;
t2 = qJD(6) * t25 + t13 * t237 + t232 * t27;
t23 = mrSges(7,1) * t76 - mrSges(7,3) * t38;
t24 = -mrSges(7,2) * t76 + mrSges(7,3) * t39;
t78 = mrSges(7,2) * t254 + mrSges(7,3) * t116;
t79 = -mrSges(7,1) * t254 + mrSges(7,3) * t253;
t245 = m(7) * (t2 * t237 - t25 * t281 - t26 * t282 - t327) + t237 * t24 - t232 * t23 - t79 * t281 - t78 * t282;
t142 = -mrSges(7,2) * t193 - mrSges(7,3) * t305;
t143 = mrSges(7,1) * t193 - mrSges(7,3) * t304;
t85 = mrSges(7,1) * t152 + mrSges(7,3) * t250;
t86 = -mrSges(7,2) * t152 - mrSges(7,3) * t251;
t244 = m(7) * (-t281 * t96 - t282 * t97 + t359) + t237 * t86 - t232 * t85 - t143 * t281 - t142 * t282;
t10 = Ifges(7,4) * t38 + Ifges(7,2) * t39 + Ifges(7,6) * t76;
t11 = Ifges(7,1) * t38 + Ifges(7,4) * t39 + Ifges(7,5) * t76;
t45 = pkin(5) * t300 - t52;
t48 = -Ifges(7,4) * t253 + Ifges(7,2) * t116 - Ifges(7,6) * t254;
t49 = -Ifges(7,1) * t253 + Ifges(7,4) * t116 - Ifges(7,5) * t254;
t243 = t2 * t323 - t254 * t341 + t14 * t205 + t17 * mrSges(6,1) + t38 * t337 + t208 * t344 + t45 * t196 + t48 * t267 + t49 * t266 + t76 * t339 + t11 * t335 + t10 * t333 + t199 * t342 - t253 * t340 + (-t327 + (-t232 * t26 - t237 * t25) * qJD(6)) * mrSges(7,3) - t16 * mrSges(6,2) - t279;
t119 = Ifges(7,6) * t193 + t194 * t260;
t120 = Ifges(7,5) * t193 + t194 * t261;
t55 = -Ifges(7,4) * t250 - Ifges(7,2) * t251 + Ifges(7,6) * t152;
t56 = -Ifges(7,1) * t250 - Ifges(7,4) * t251 + Ifges(7,5) * t152;
t242 = t120 * t266 - t293 * t337 + t152 * t339 + t165 * t196 - t199 * t305 / 0.2e1 + t304 * t340 + t193 * t341 + t56 * t335 + t55 * t333 - t106 * mrSges(6,2) + t290 - t251 * t208 / 0.2e1 + (t194 * t210 + t119) * t267 + t360 * t107 + ((-t232 * t97 - t237 * t96) * qJD(6) + t359) * mrSges(7,3);
t225 = Ifges(5,5) * t284;
t222 = -pkin(4) * t238 - pkin(5);
t221 = pkin(4) * t233 + pkin(12);
t216 = Ifges(4,5) * t271;
t211 = Ifges(5,1) * t234 + t321;
t209 = Ifges(5,2) * t239 + t322;
t202 = (Ifges(5,1) * t239 - t322) * qJD(4);
t200 = (-Ifges(5,2) * t234 + t321) * qJD(4);
t197 = (mrSges(5,1) * t234 + mrSges(5,2) * t239) * qJD(4);
t192 = -mrSges(4,2) * t230 + mrSges(4,3) * t300;
t178 = (mrSges(4,1) * t235 + mrSges(4,2) * t240) * t287;
t164 = -mrSges(5,1) * t300 - mrSges(5,3) * t184;
t163 = mrSges(5,2) * t300 + mrSges(5,3) * t183;
t157 = Ifges(6,1) * t194 - Ifges(6,4) * t193;
t156 = Ifges(6,4) * t194 - Ifges(6,2) * t193;
t155 = mrSges(6,1) * t193 + mrSges(6,2) * t194;
t138 = t262 * t194;
t132 = -mrSges(5,2) * t272 + mrSges(5,3) * t159;
t131 = mrSges(5,1) * t272 - mrSges(5,3) * t158;
t129 = Ifges(5,1) * t184 + Ifges(5,4) * t183 - Ifges(5,5) * t300;
t128 = Ifges(5,4) * t184 + Ifges(5,2) * t183 - Ifges(5,6) * t300;
t121 = mrSges(6,2) * t300 + mrSges(6,3) * t254;
t118 = Ifges(7,3) * t193 + (Ifges(7,5) * t237 - t318) * t194;
t104 = -mrSges(5,1) * t159 + mrSges(5,2) * t158;
t102 = -Ifges(6,1) * t151 - Ifges(6,4) * t152;
t101 = -Ifges(6,4) * t151 - Ifges(6,2) * t152;
t100 = mrSges(6,1) * t152 - mrSges(6,2) * t151;
t92 = Ifges(5,1) * t158 + Ifges(5,4) * t159 + Ifges(5,5) * t272;
t91 = Ifges(5,4) * t158 + Ifges(5,2) * t159 + Ifges(5,6) * t272;
t87 = -mrSges(6,1) * t254 + mrSges(6,2) * t134;
t84 = Ifges(6,1) * t134 + Ifges(6,4) * t254 - Ifges(6,5) * t300;
t83 = Ifges(6,4) * t134 + Ifges(6,2) * t254 - Ifges(6,6) * t300;
t68 = -mrSges(6,2) * t272 - mrSges(6,3) * t76;
t47 = -Ifges(7,5) * t253 + Ifges(7,6) * t116 - Ifges(7,3) * t254;
t30 = mrSges(6,1) * t76 + mrSges(6,2) * t75;
t29 = Ifges(6,1) * t75 - Ifges(6,4) * t76 + Ifges(6,5) * t272;
t28 = Ifges(6,4) * t75 - Ifges(6,2) * t76 + Ifges(6,6) * t272;
t1 = [0.2e1 * m(7) * (t43 * t6 + t44 * t5 - t329) + 0.2e1 * m(6) * (t20 * t70 - t329 + t82) + 0.2e1 * m(5) * (t113 * t59 + t114 * t60 + t82) + 0.2e1 * m(4) * (t109 * t148 + t182 * t264 + t82); t109 * t192 + t113 * t131 + t114 * t132 + t20 * t121 + t60 * t163 + t59 * t164 + t182 * t178 + t43 * t23 + t44 * t24 + t5 * t78 + t6 * t79 + t70 * t68 - t325 * t255 + t311 * t21 + (t104 + t30) * t147 + (-mrSges(3,1) * t236 - mrSges(3,2) * t241) * t288 + (t87 + t289) * t108 + ((-mrSges(4,1) * t240 + mrSges(4,2) * t235) * t264 + (t147 * t240 - t148 * t235) * qJD(3) * mrSges(4,3)) * t228 + m(4) * (-pkin(2) * t228 ^ 2 * t273 - t108 * t188 + t109 * t189 + t148 * t180 + t307) + m(5) * (t108 * t174 + t113 * t81 + t114 * t80 + t124 * t59 + t125 * t60 + t307) + m(6) * (t108 * t136 + t123 * t147 + t16 * t70 + t17 * t255 + t20 * t324 - t21 * t52) + m(7) * (-t14 * t255 + t2 * t44 + t21 * t45 + t25 * t6 + t26 * t5 + t3 * t43); 0.2e1 * t289 * t181 + 0.2e1 * t324 * t68 + (t123 * t136 + t16 * t324 + t17 * t52) * t353 + t183 * t91 + t184 * t92 + 0.2e1 * t180 * t192 + 0.2e1 * t174 * t104 + t159 * t128 + 0.2e1 * t80 * t163 + 0.2e1 * t81 * t164 + t158 * t129 + t134 * t29 + 0.2e1 * t136 * t30 + 0.2e1 * t124 * t131 + 0.2e1 * t125 * t132 + t116 * t10 + 0.2e1 * t16 * t121 + 0.2e1 * t17 * t122 + 0.2e1 * t123 * t87 - t253 * t11 + (-0.2e1 * pkin(2) * t178 + (t279 + t354) * t240 + ((0.2e1 * Ifges(4,4) * t300 + Ifges(4,5) * t230 + t188 * t350) * t240 + (-0.2e1 * Ifges(4,4) * t301 + t189 * t350 + Ifges(5,5) * t184 + Ifges(6,5) * t134 - 0.2e1 * Ifges(4,6) * t230 + Ifges(5,6) * t183 + Ifges(6,6) * t254 + ((2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3) - Ifges(6,3)) * t300) * t235) * qJD(3)) * t228 - (t9 - t28) * t254 + t75 * t84 + 0.2e1 * t2 * t78 + 0.2e1 * t3 * t79 + 0.2e1 * t14 * t71 + 0.2e1 * t52 * t67 + t39 * t48 + t38 * t49 + 0.2e1 * t45 * t22 + 0.2e1 * t25 * t23 + 0.2e1 * t26 * t24 + (t14 * t45 + t2 * t26 + t25 * t3) * t352 + 0.2e1 * m(5) * (t124 * t81 + t125 * t80 + t174 * t181) + 0.2e1 * m(4) * (t180 * t189 - t181 * t188) + (-t83 + t47) * t76 + t230 * t216; -t109 * mrSges(4,2) + t21 * t138 + t5 * t142 + t6 * t143 + t43 * t85 + t44 * t86 - t255 * t77 + (t100 + t197) * t147 + m(7) * (t32 * t44 + t33 * t43 + t5 * t97 + t6 * t96 + t259) + m(6) * (t106 * t70 + t147 * t278 + t166 * t20 + t259) + (t151 * t255 - t152 * t70 - t193 * t20 + t194 * t21) * mrSges(6,3) + (m(5) * pkin(10) + mrSges(5,3)) * (-t59 * t234 + t60 * t239 + (-t113 * t239 - t114 * t234) * qJD(4)) + (m(6) * t223 + t155 + t358) * t108; t358 * t181 + ((Ifges(6,5) * t194 / 0.2e1 - Ifges(6,6) * t193 / 0.2e1 - Ifges(4,6) + Ifges(5,5) * t334 + Ifges(5,6) * t332) * t286 - (-Ifges(5,6) * t285 + t225 + t290) * t240 / 0.2e1) * t228 + (t47 / 0.2e1 - t83 / 0.2e1 - t324 * mrSges(6,3)) * t152 + m(6) * (t106 * t324 - t107 * t52 + t123 * t223 + t136 * t278 + t16 * t166 - t165 * t17) + t223 * t30 + t159 * t209 / 0.2e1 + t158 * t211 / 0.2e1 + t174 * t197 + t183 * t200 / 0.2e1 + t184 * t202 / 0.2e1 - t180 * mrSges(4,2) + t166 * t68 + t123 * t155 + t75 * t157 / 0.2e1 + t2 * t142 + t3 * t143 + t134 * t102 / 0.2e1 + t136 * t100 + t14 * t138 + t38 * t120 / 0.2e1 + t106 * t121 - t253 * t56 / 0.2e1 - (t54 / 0.2e1 - t101 / 0.2e1) * t254 - pkin(3) * t104 + t96 * t23 + t97 * t24 + t25 * t85 + t26 * t86 + t45 * t77 + t32 * t78 + t33 * t79 + t119 * t344 + t216 + t91 * t332 + t92 * t334 + t55 * t342 + (t9 / 0.2e1 - t28 / 0.2e1 - t16 * mrSges(6,3)) * t193 + m(7) * (t107 * t45 + t14 * t165 + t2 * t97 + t25 * t33 + t26 * t32 + t3 * t96) + (-t156 / 0.2e1 + t118 / 0.2e1) * t76 + (t129 * t332 + (-t128 / 0.2e1 + pkin(4) * t87) * t234) * qJD(4) - (t49 * t333 + t48 * t336 + t84 / 0.2e1 - t52 * mrSges(6,3)) * t151 + (t11 * t333 + t29 / 0.2e1 + t10 * t336 - t17 * mrSges(6,3) + (t49 * t336 - t237 * t48 / 0.2e1) * qJD(6)) * t194 + ((-t124 * t239 - t125 * t234) * qJD(4) + t257) * mrSges(5,3) + (m(5) * (-t124 * t284 - t125 * t285 + t257) + t239 * t132 - t234 * t131 - t163 * t285 - t164 * t284) * pkin(10) + t311 * t107 + t325 * t165; -0.2e1 * pkin(3) * t197 + 0.2e1 * t223 * t100 + t138 * t348 + 0.2e1 * t32 * t142 + 0.2e1 * t33 * t143 + t77 * t347 + t239 * t200 + t234 * t202 + 0.2e1 * t96 * t85 + 0.2e1 * t97 * t86 + (t239 * t211 + (t155 * t351 - t209) * t234) * qJD(4) + (t106 * t166 + t223 * t278 + t308) * t353 + (t32 * t97 + t33 * t96 + t308) * t352 + (t106 * t349 - t101 + t54) * t193 + (t166 * t349 + t118 - t156) * t152 - (mrSges(6,3) * t347 - t119 * t232 + t120 * t237 + t157) * t151 + (mrSges(6,3) * t348 - t232 * t55 + t237 * t56 + t102 + (-t119 * t237 - t120 * t232) * qJD(6)) * t194; m(7) * t21 * t222 - t60 * mrSges(5,2) + t59 * mrSges(5,1) + ((t20 * t233 - t21 * t238) * t346 + (m(7) * (t292 * t44 - t298 * t43 - t315) / 0.2e1 + (t238 * t70 - t315) * t346) * qJD(5)) * t351 + t246 + m(7) * (-t281 * t43 - t282 * t44 - t326 + t328) * t221; t243 - t80 * mrSges(5,2) + t81 * mrSges(5,1) + (m(6) * (t16 * t233 + t17 * t238) + t238 * t67 + t233 * t68 + (t311 * t233 + (-t232 * t79 + t237 * t78 + t121) * t238 + m(7) * (t233 * t45 - t25 * t298 + t26 * t292) + m(6) * (-t233 * t52 + t238 * t324)) * qJD(5)) * pkin(4) + t275 * t222 + t245 * t221 - t354; t242 + (m(6) * (t106 * t233 - t107 * t238) + (t151 * t238 - t152 * t233) * mrSges(6,3) + ((mrSges(6,3) * t194 + t138) * t233 + (-mrSges(6,3) * t193 + t142 * t237 - t143 * t232) * t238 + m(7) * (t292 * t97 - t298 * t96 + t306) + m(6) * (t166 * t238 + t306)) * qJD(5)) * pkin(4) + t244 * t221 + t268 * t222 + t225 + (-Ifges(5,6) * t234 + pkin(10) * t206) * qJD(4); 0.2e1 * t222 * t196 + (-0.2e1 * t313 - 0.2e1 * t316 + (t221 * t356 + t222 * t233) * t352 + 0.2e1 * t297 + 0.2e1 * t263) * t317 + t249; m(7) * (-pkin(5) * t21 + (t248 + t328) * pkin(12)) + t246; -pkin(5) * t275 + pkin(12) * t245 + t243; -pkin(5) * t268 + pkin(12) * t244 + t242; (t222 - pkin(5)) * t196 + (-t313 - t316 + m(7) * (-pkin(5) * t233 + pkin(12) * t356) + t297 + t263) * t317 + t249; -0.2e1 * pkin(5) * t196 + t249; mrSges(7,1) * t6 - mrSges(7,2) * t5; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t9; mrSges(7,1) * t33 - mrSges(7,2) * t32 + t54; t224 - t262 * pkin(4) * t283 + (t205 * t221 - t318) * qJD(6); t224 + (pkin(12) * t205 - t318) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
