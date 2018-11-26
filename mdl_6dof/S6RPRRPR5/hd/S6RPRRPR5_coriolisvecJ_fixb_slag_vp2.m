% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2018-11-23 16:18
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:18:13
% EndTime: 2018-11-23 16:18:22
% DurationCPUTime: 9.33s
% Computational Cost: add. (10006->506), mult. (26753->657), div. (0->0), fcn. (20705->8), ass. (0->245)
t355 = Ifges(5,4) + Ifges(6,6);
t180 = sin(pkin(10));
t181 = cos(pkin(10));
t291 = sin(qJ(3));
t293 = cos(qJ(3));
t214 = t180 * t291 - t181 * t293;
t292 = cos(qJ(4));
t205 = t292 * t214;
t151 = qJD(1) * t205;
t165 = t180 * t293 + t181 * t291;
t160 = t165 * qJD(1);
t183 = sin(qJ(4));
t132 = t160 * t183 + t151;
t179 = qJD(3) + qJD(4);
t182 = sin(qJ(6));
t184 = cos(qJ(6));
t115 = t132 * t182 + t179 * t184;
t304 = t115 / 0.2e1;
t360 = Ifges(7,5) * t304;
t295 = t179 / 0.2e1;
t211 = qJD(1) * t214;
t204 = t183 * t211;
t192 = t160 * t292 - t204;
t298 = t192 / 0.2e1;
t299 = -t192 / 0.2e1;
t336 = -Ifges(5,5) + Ifges(6,4);
t359 = -Ifges(5,1) * t298 + Ifges(6,2) * t299 + t336 * t295;
t114 = t132 * t184 - t179 * t182;
t358 = t114 / 0.2e1;
t128 = qJD(6) + t192;
t357 = t128 / 0.2e1;
t302 = -t132 / 0.2e1;
t354 = Ifges(6,3) + Ifges(5,2);
t343 = Ifges(7,6) * t358 + Ifges(7,3) * t357;
t241 = -t181 * pkin(2) - pkin(1);
t168 = qJD(1) * t241 + qJD(2);
t139 = pkin(3) * t211 + t168;
t307 = pkin(4) + pkin(9);
t337 = pkin(5) * t192;
t274 = pkin(7) + qJ(2);
t169 = t274 * t180;
t166 = qJD(1) * t169;
t170 = t274 * t181;
t167 = qJD(1) * t170;
t137 = -t293 * t166 - t167 * t291;
t112 = -t160 * pkin(8) + t137;
t110 = qJD(3) * pkin(3) + t112;
t138 = -t166 * t291 + t167 * t293;
t113 = -pkin(8) * t211 + t138;
t251 = t183 * t113;
t63 = -t292 * t110 + t251;
t218 = t63 + t337;
t345 = qJD(5) + t218;
t38 = -t179 * t307 + t345;
t65 = t132 * pkin(4) - qJ(5) * t192 + t139;
t46 = t132 * pkin(9) + t65;
t14 = -t182 * t46 + t184 * t38;
t15 = t182 * t38 + t184 * t46;
t351 = t302 * t355 + t360;
t330 = -qJD(5) - t63;
t61 = -pkin(4) * t179 - t330;
t344 = t61 * mrSges(6,1) + t14 * mrSges(7,1) + t139 * mrSges(5,2) - t15 * mrSges(7,2) + t63 * mrSges(5,3) - t65 * mrSges(6,3) + t343 + t351 - t359;
t352 = t343 + t344;
t301 = t132 / 0.2e1;
t296 = -t179 / 0.2e1;
t350 = mrSges(6,2) - mrSges(5,1);
t349 = mrSges(5,3) + mrSges(6,1);
t233 = qJD(4) * t292;
t67 = t112 * t292 - t251;
t348 = pkin(3) * t233 - t67;
t287 = t132 * pkin(5);
t335 = Ifges(6,5) - Ifges(5,6);
t333 = -qJD(5) - t348;
t257 = qJ(5) * t132;
t220 = t14 * t182 - t15 * t184;
t209 = qJD(3) * t214;
t202 = qJD(1) * t209;
t236 = qJD(2) * t291;
t227 = qJD(1) * t236;
t237 = qJD(2) * t293;
t228 = qJD(1) * t237;
t234 = qJD(3) * t291;
t235 = qJD(3) * t293;
t187 = pkin(8) * t202 + t166 * t234 - t167 * t235 - t180 * t228 - t181 * t227;
t109 = t292 * t113;
t64 = t183 * t110 + t109;
t107 = -t166 * t235 - t167 * t234 - t180 * t227 + t181 * t228;
t210 = qJD(3) * t165;
t203 = qJD(1) * t210;
t97 = -pkin(8) * t203 + t107;
t28 = qJD(4) * t64 + t183 * t97 - t292 * t187;
t196 = qJD(3) * t205;
t250 = qJD(4) * t183;
t87 = qJD(1) * t196 + qJD(4) * t151 + t160 * t250 + t183 * t203;
t11 = -t87 * pkin(5) + t28;
t198 = pkin(3) * t203;
t197 = t292 * t210;
t88 = qJD(1) * t197 - qJD(4) * t204 + t160 * t233 - t183 * t202;
t31 = t88 * pkin(4) + t87 * qJ(5) - qJD(5) * t192 + t198;
t13 = t88 * pkin(9) + t31;
t1 = qJD(6) * t14 + t11 * t182 + t13 * t184;
t256 = qJD(6) * t15;
t2 = t11 * t184 - t13 * t182 - t256;
t324 = -t1 * t182 - t184 * t2;
t193 = m(7) * (-qJD(6) * t220 - t324);
t54 = qJD(6) * t114 + t182 * t88;
t32 = -mrSges(7,1) * t87 - mrSges(7,3) * t54;
t55 = -qJD(6) * t115 + t184 * t88;
t33 = mrSges(7,2) * t87 + mrSges(7,3) * t55;
t347 = t182 * t33 + t184 * t32 + t193;
t226 = mrSges(7,1) * t184 - mrSges(7,2) * t182;
t62 = -t179 * qJ(5) - t64;
t39 = -t62 - t287;
t266 = t115 * Ifges(7,4);
t50 = t114 * Ifges(7,2) + t128 * Ifges(7,6) + t266;
t346 = t39 * t226 - t184 * t50 / 0.2e1;
t342 = (m(3) * qJ(2) + mrSges(3,3)) * (-t180 ^ 2 - t181 ^ 2);
t341 = t139 * mrSges(5,1) + t62 * mrSges(6,1) - t65 * mrSges(6,2) - t64 * mrSges(5,3) + Ifges(6,5) * t295 + Ifges(5,6) * t296 + t299 * t355 + t301 * t354;
t314 = t54 / 0.2e1;
t313 = t55 / 0.2e1;
t309 = -t87 / 0.2e1;
t340 = m(5) * t63;
t338 = pkin(4) * t192;
t332 = t337 - t333;
t331 = t192 * t307;
t116 = -mrSges(5,2) * t179 - mrSges(5,3) * t132;
t118 = mrSges(6,1) * t132 - mrSges(6,3) * t179;
t329 = t116 - t118;
t328 = -t179 * t350 - t192 * t349;
t248 = qJD(6) * t184;
t326 = -t192 * t184 - t248;
t323 = m(6) * t61 - t328;
t322 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t321 = -m(5) * t64 + m(6) * t62 - t329;
t221 = t14 * t184 + t15 * t182;
t320 = -m(7) * t221 - t323;
t223 = Ifges(7,5) * t182 + Ifges(7,6) * t184;
t269 = Ifges(7,4) * t182;
t224 = Ifges(7,2) * t184 + t269;
t268 = Ifges(7,4) * t184;
t225 = Ifges(7,1) * t182 + t268;
t294 = -t182 / 0.2e1;
t303 = -t128 / 0.2e1;
t305 = -t115 / 0.2e1;
t306 = -t114 / 0.2e1;
t111 = Ifges(7,4) * t114;
t51 = t115 * Ifges(7,1) + t128 * Ifges(7,5) + t111;
t316 = t223 * t303 + t224 * t306 + t225 * t305 + t51 * t294 - t341 + t346;
t315 = Ifges(7,1) * t314 + Ifges(7,4) * t313 + Ifges(7,5) * t309;
t290 = pkin(3) * t160;
t289 = pkin(3) * t183;
t140 = -t293 * t169 - t170 * t291;
t124 = -t165 * pkin(8) + t140;
t141 = -t291 * t169 + t293 * t170;
t125 = -pkin(8) * t214 + t141;
t74 = -t292 * t124 + t125 * t183;
t284 = t28 * t74;
t279 = t87 * mrSges(6,1);
t278 = t87 * mrSges(5,3);
t277 = t88 * mrSges(6,1);
t276 = t88 * mrSges(5,3);
t272 = mrSges(4,3) * t160;
t271 = Ifges(4,4) * t160;
t270 = Ifges(4,4) * t165;
t212 = t183 * t214;
t96 = -qJD(4) * t212 + t165 * t233 - t183 * t209 + t197;
t261 = t182 * t96;
t259 = t184 * t96;
t68 = -mrSges(7,1) * t114 + mrSges(7,2) * t115;
t258 = t118 - t68;
t255 = t192 * t182;
t135 = t165 * t183 + t205;
t253 = t135 * t182;
t252 = t135 * t184;
t249 = qJD(6) * t182;
t247 = Ifges(7,5) * t54 + Ifges(7,6) * t55 - Ifges(7,3) * t87;
t246 = t292 * pkin(3);
t239 = t88 * mrSges(5,1) - t87 * mrSges(5,2);
t238 = -t88 * mrSges(6,2) + t87 * mrSges(6,3);
t232 = -t249 / 0.2e1;
t66 = t112 * t183 + t109;
t176 = -t246 - pkin(4);
t222 = t290 + t257;
t136 = t165 * t292 - t212;
t147 = pkin(3) * t214 + t241;
t191 = -t136 * qJ(5) + t147;
t57 = t135 * t307 + t191;
t58 = pkin(5) * t136 + t74;
t30 = t182 * t58 + t184 * t57;
t29 = -t182 * t57 + t184 * t58;
t71 = -mrSges(7,2) * t128 + mrSges(7,3) * t114;
t72 = mrSges(7,1) * t128 - mrSges(7,3) * t115;
t219 = -t182 * t72 + t184 * t71;
t217 = t135 * t248 + t261;
t216 = t135 * t249 - t259;
t75 = t183 * t124 + t125 * t292;
t27 = t110 * t233 - t113 * t250 + t183 * t187 + t292 * t97;
t121 = -t169 * t235 - t170 * t234 - t180 * t236 + t181 * t237;
t103 = -pkin(8) * t210 + t121;
t188 = pkin(8) * t209 + t169 * t234 - t170 * t235 - t180 * t237 - t181 * t236;
t34 = -t292 * t103 - t124 * t233 + t125 * t250 - t183 * t188;
t213 = Ifges(4,4) * t214;
t208 = t165 * qJD(2);
t24 = -qJD(5) * t179 - t27;
t207 = pkin(3) * t210;
t206 = mrSges(4,3) * t211;
t200 = -t210 / 0.2e1;
t199 = -t209 / 0.2e1;
t35 = qJD(4) * t75 + t183 * t103 - t292 * t188;
t194 = mrSges(4,1) * t203 - mrSges(4,2) * t202;
t95 = qJD(4) * t205 + t165 * t250 + t183 * t210 + t196;
t36 = t96 * pkin(4) + t95 * qJ(5) - t136 * qJD(5) + t207;
t10 = -pkin(5) * t88 - t24;
t7 = t54 * Ifges(7,4) + t55 * Ifges(7,2) - t87 * Ifges(7,6);
t190 = -t27 * mrSges(5,2) - t24 * mrSges(6,3) + t184 * t315 + t7 * t294 + t14 * mrSges(7,3) * t249 + (Ifges(7,1) * t184 - t269) * t314 + (-Ifges(7,2) * t182 + t268) * t313 + t51 * t232 + (Ifges(7,5) * t184 - Ifges(7,6) * t182) * t309 + t10 * (mrSges(7,1) * t182 + mrSges(7,2) * t184) + t335 * t88 + t336 * t87 + t350 * t28 + t346 * qJD(6) - (t114 * t224 + t115 * t225 + t128 * t223) * qJD(6) / 0.2e1;
t189 = qJD(3) * (-Ifges(4,1) * t214 - t270) / 0.2e1;
t175 = qJ(5) + t289;
t156 = Ifges(4,4) * t211;
t146 = qJD(3) * mrSges(4,1) - t272;
t145 = -qJD(3) * mrSges(4,2) - t206;
t130 = Ifges(4,1) * t160 + Ifges(4,5) * qJD(3) - t156;
t129 = -Ifges(4,2) * t211 + Ifges(4,6) * qJD(3) + t271;
t122 = -qJD(3) * t141 - t208;
t108 = -qJD(1) * t208 - qJD(3) * t138;
t91 = -mrSges(6,2) * t132 - mrSges(6,3) * t192;
t90 = mrSges(5,1) * t132 + mrSges(5,2) * t192;
t89 = t257 + t338;
t73 = t135 * pkin(4) + t191;
t69 = t222 + t338;
t60 = t257 + t331;
t59 = -t135 * pkin(5) + t75;
t56 = t222 + t331;
t47 = t66 - t287;
t45 = t64 - t287;
t23 = t182 * t45 + t184 * t60;
t22 = -t182 * t60 + t184 * t45;
t21 = t96 * pkin(9) + t36;
t20 = -mrSges(7,1) * t55 + mrSges(7,2) * t54;
t19 = t182 * t47 + t184 * t56;
t18 = -t182 * t56 + t184 * t47;
t17 = -t95 * pkin(5) + t35;
t16 = -pkin(5) * t96 - t34;
t4 = -qJD(6) * t30 + t17 * t184 - t182 * t21;
t3 = qJD(6) * t29 + t17 * t182 + t184 * t21;
t5 = [t73 * t238 + t147 * t239 + (m(5) * t27 - m(6) * t24 - t276 - t277) * t75 + m(7) * (t1 * t30 + t10 * t59 + t14 * t4 + t15 * t3 + t16 * t39 + t2 * t29) + m(4) * (t107 * t141 + t108 * t140 + t121 * t138 + t122 * t137) + (qJD(6) * t51 + t7) * t252 / 0.2e1 + m(6) * (t31 * t73 + t36 * t65 + t284) + (-t107 * t214 - t108 * t165 + t137 * t209 - t138 * t210 + t140 * t202 - t141 * t203) * mrSges(4,3) + (-Ifges(5,4) * t302 + t301 * Ifges(6,6) - t352 + t359 - t360) * t95 + (-t278 - t279) * t74 + t51 * t261 / 0.2e1 + t50 * t259 / 0.2e1 + (t1 * t252 - t14 * t217 - t15 * t216 - t2 * t253) * mrSges(7,3) + (t87 * t135 + t299 * t96) * Ifges(6,6) + t130 * t199 + t129 * t200 - (-Ifges(4,2) * t165 - t213) * t202 + (Ifges(7,1) * t217 - Ifges(7,4) * t216) * t304 + (mrSges(5,1) * t198 + t24 * mrSges(6,1) - t31 * mrSges(6,2) - t27 * mrSges(5,3) - t10 * t226 + t223 * t309 + t224 * t313 + t225 * t314 + t50 * t232 + t354 * t88 + (t87 / 0.2e1 - t309) * Ifges(5,4)) * t135 + (t349 * t28 - t31 * mrSges(6,3) + mrSges(5,2) * t198 + Ifges(7,6) * t313 + Ifges(7,5) * t314 + (Ifges(7,3) + Ifges(5,1)) * t309 + t322 + t247 / 0.2e1 + (-Ifges(6,2) - Ifges(5,1) / 0.2e1) * t87 - t355 * t88) * t136 + t39 * (mrSges(7,1) * t216 + mrSges(7,2) * t217) + qJD(3) ^ 2 * (-Ifges(4,5) * t214 - Ifges(4,6) * t165) / 0.2e1 + (t323 + t340) * t35 + t121 * t145 + t122 * t146 + t90 * t207 + t321 * t34 + (-Ifges(5,4) * t298 - Ifges(5,2) * t302 + Ifges(6,3) * t301 + t335 * t295 + t341) * t96 + t36 * t91 + t3 * t71 + t4 * t72 + t16 * t68 + t59 * t20 + t29 * t32 + t30 * t33 + (-0.2e1 * t342 * qJD(2) + t165 * t189 + (-Ifges(4,2) * t214 + t270) * t200 + (Ifges(4,1) * t165 - t213) * t199) * qJD(1) + (Ifges(7,5) * t217 - Ifges(7,6) * t216) * t357 + (Ifges(7,4) * t217 - Ifges(7,2) * t216) * t358 + t241 * t194 + t160 * t189 + t253 * t315 + m(5) * (t139 * t207 + t147 * t198 + t284) + t168 * (mrSges(4,1) * t165 - mrSges(4,2) * t214) * qJD(3); t238 + t239 + m(6) * t31 + m(7) * (-qJD(6) * t221 + t1 * t184 - t182 * t2) + t184 * t33 - t182 * t32 + t145 * t211 - m(4) * (-t137 * t160 - t138 * t211) + t160 * t146 + m(5) * t198 + t194 + t326 * t72 + (-t255 - t249) * t71 + (t320 - t340) * t192 - (-m(7) * t39 + t321 - t68) * t132 + t342 * qJD(1) ^ 2; (t146 + t272) * t138 + (-t277 + t20) * t175 + (-Ifges(4,2) * t160 + t130 - t156) * t211 / 0.2e1 - Ifges(4,5) * t202 - Ifges(4,6) * t203 + (t182 * t71 + t184 * t72 - t320) * pkin(3) * t250 + (-t206 - t145) * t137 - t176 * t279 + t348 * t116 - t276 * t289 - t90 * t290 - t160 * (-Ifges(4,1) * t211 - t271) / 0.2e1 - t168 * (t160 * mrSges(4,1) - mrSges(4,2) * t211) - qJD(3) * (-Ifges(4,5) * t211 - Ifges(4,6) * t160) / 0.2e1 + (-Ifges(5,4) * t299 - Ifges(5,2) * t301 + Ifges(6,6) * t298 + Ifges(6,3) * t302 + t296 * t335 + t316) * t192 + t332 * t68 + (t10 * t175 - t14 * t18 - t15 * t19 + t332 * t39) * m(7) + t333 * t118 + (-t175 * t24 + t176 * t28 + t333 * t62 - t61 * t66 - t65 * t69) * m(6) + (t14 * t255 + t15 * t326 + t324) * mrSges(7,3) + t328 * t66 + t160 * t129 / 0.2e1 + (t248 * t71 - t249 * t72 + t347) * (-pkin(9) + t176) + (-Ifges(5,1) * t299 - Ifges(5,4) * t301 - Ifges(7,5) * t305 + Ifges(6,2) * t298 + Ifges(6,6) * t302 - Ifges(7,6) * t306 - Ifges(7,3) * t303 + t336 * t296 + t344) * t132 - t107 * mrSges(4,2) + t108 * mrSges(4,1) - t69 * t91 - t19 * t71 - t18 * t72 + (-t139 * t290 - t63 * t66 - t64 * t67 + (-t292 * t28 + t183 * t27 + (t183 * t63 + t292 * t64) * qJD(4)) * pkin(3)) * m(5) + t190 + t246 * t278; (-(qJD(6) * t71 + t32) * t307 + (-t2 - t256) * mrSges(7,3)) * t184 + t329 * t63 + ((-Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t179 + t351 + t352) * t132 + (pkin(4) * t87 - qJ(5) * t88) * mrSges(6,1) + (t316 + (-Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1) * t179 + (-Ifges(6,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t132 + t220 * mrSges(7,3) + (Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t192) * t192 + (-t1 * mrSges(7,3) - (-qJD(6) * t72 + t33) * t307) * t182 + t328 * t64 - t258 * qJD(5) - t89 * t91 - t23 * t71 - t22 * t72 + t218 * t68 + qJ(5) * t20 - t307 * t193 + t190 + (qJ(5) * t10 - t14 * t22 - t15 * t23 + t345 * t39) * m(7) + (-pkin(4) * t28 - qJ(5) * t24 + t330 * t62 - t61 * t64 - t65 * t89) * m(6); -t279 + t258 * t179 + t219 * qJD(6) + (t219 + t91) * t192 - m(7) * (t179 * t39 + t192 * t220) + (t179 * t62 + t192 * t65 + t28) * m(6) + t347; -t39 * (mrSges(7,1) * t115 + mrSges(7,2) * t114) + (Ifges(7,1) * t114 - t266) * t305 + t50 * t304 + (Ifges(7,5) * t114 - Ifges(7,6) * t115) * t303 - t14 * t71 + t15 * t72 + (t114 * t14 + t115 * t15) * mrSges(7,3) + t247 + (-Ifges(7,2) * t115 + t111 + t51) * t306 + t322;];
tauc  = t5(:);
