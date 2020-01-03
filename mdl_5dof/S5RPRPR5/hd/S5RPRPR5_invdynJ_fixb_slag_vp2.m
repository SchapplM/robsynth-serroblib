% Calculate vector of inverse dynamics joint torques for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:41:27
% EndTime: 2020-01-03 11:42:07
% DurationCPUTime: 11.53s
% Computational Cost: add. (4906->497), mult. (12066->690), div. (0->0), fcn. (8654->14), ass. (0->264)
t215 = cos(pkin(8));
t281 = qJD(1) * t215;
t188 = qJD(3) - t281;
t274 = qJDD(1) * qJ(2);
t276 = qJD(1) * qJD(2);
t183 = t274 + t276;
t357 = t183 + t276;
t356 = Ifges(4,3) + Ifges(5,3);
t212 = sin(pkin(9));
t214 = cos(pkin(9));
t218 = sin(qJ(3));
t221 = cos(qJ(3));
t166 = t212 * t221 + t214 * t218;
t153 = t166 * qJD(3);
t233 = qJD(1) * t166;
t339 = t215 * t233 - t153;
t237 = t212 * t218 - t214 * t221;
t338 = t188 * t237;
t213 = sin(pkin(8));
t209 = t213 ^ 2;
t355 = t357 * t209;
t322 = m(5) + m(6);
t354 = t322 + m(4) + m(3);
t326 = m(5) * pkin(3);
t272 = qJDD(1) * t215;
t186 = qJDD(3) - t272;
t316 = t186 / 0.2e1;
t275 = qJD(1) * qJD(3);
t147 = (qJDD(1) * t221 - t218 * t275) * t213;
t148 = (-qJDD(1) * t218 - t221 * t275) * t213;
t82 = t147 * t214 + t148 * t212;
t352 = Ifges(5,5) * t82;
t217 = sin(qJ(5));
t220 = cos(qJ(5));
t282 = qJD(1) * t213;
t259 = t221 * t282;
t260 = t218 * t282;
t123 = t212 * t260 - t214 * t259;
t124 = t213 * t233;
t250 = t123 * t217 - t124 * t220;
t81 = -t147 * t212 + t148 * t214;
t21 = qJD(5) * t250 + t217 * t81 + t220 * t82;
t351 = Ifges(6,5) * t21;
t350 = Ifges(5,6) * t81;
t64 = t123 * t220 + t124 * t217;
t22 = qJD(5) * t64 - t217 * t82 + t220 * t81;
t349 = Ifges(6,6) * t22;
t94 = -t166 * t217 - t220 * t237;
t348 = qJD(5) * t94 + t217 * t339 - t220 * t338;
t95 = t166 * t220 - t217 * t237;
t347 = -qJD(5) * t95 + t217 * t338 + t220 * t339;
t346 = Ifges(4,5) * t147;
t345 = Ifges(4,6) * t148;
t178 = qJDD(5) + t186;
t344 = Ifges(6,3) * t178;
t343 = -mrSges(4,1) - t326;
t312 = pkin(3) * t214;
t196 = pkin(4) + t312;
t313 = pkin(3) * t212;
t150 = t196 * t217 + t220 * t313;
t309 = pkin(7) * t124;
t245 = pkin(2) * t215 + pkin(6) * t213;
t174 = -pkin(1) - t245;
t151 = qJD(1) * t174 + qJD(2);
t263 = qJ(2) * t281;
t105 = t151 * t218 + t221 * t263;
t89 = -qJ(4) * t260 + t105;
t301 = t214 * t89;
t132 = t221 * t151;
t296 = t213 * t221;
t266 = qJ(4) * t296;
t295 = t215 * t218;
t267 = qJ(2) * t295;
t232 = -t266 - t267;
t88 = qJD(1) * t232 + t132;
t45 = -t212 * t88 - t301;
t29 = t45 + t309;
t310 = pkin(7) * t123;
t83 = t212 * t89;
t46 = t214 * t88 - t83;
t30 = t46 + t310;
t342 = -qJD(5) * t150 + t217 * t30 - t220 * t29;
t149 = t196 * t220 - t217 * t313;
t341 = qJD(5) * t149 - t217 * t29 - t220 * t30;
t219 = sin(qJ(1));
t222 = cos(qJ(1));
t340 = t213 * (-g(2) * t219 + g(3) * t222);
t293 = t215 * t221;
t187 = qJ(2) * t293;
t122 = t174 * t218 + t187;
t210 = t215 ^ 2;
t337 = t209 + t210;
t336 = t186 * t356 + t345 + t346 + t350 + t352;
t335 = 0.2e1 * t316;
t333 = t210 * t357;
t211 = qJ(3) + pkin(9);
t200 = sin(t211);
t311 = pkin(3) * t218;
t170 = pkin(4) * t200 + t311;
t332 = m(6) * t170 - mrSges(2,2) + mrSges(3,3);
t304 = Ifges(4,4) * t221;
t305 = Ifges(4,4) * t218;
t331 = (qJ(2) * (mrSges(4,1) * t221 - mrSges(4,2) * t218) - t218 * (-Ifges(4,2) * t221 - t305) / 0.2e1 + t221 * (-Ifges(4,1) * t218 - t304) / 0.2e1) * t209;
t146 = qJDD(1) * t174 + qJDD(2);
t131 = t221 * t146;
t262 = qJ(2) * qJD(3) * t215;
t277 = qJD(4) * t213;
t227 = qJD(1) * (-t262 - t277);
t36 = pkin(3) * t186 - qJ(4) * t147 + t131 + (-qJD(3) * t151 - t183 * t215) * t218 + t221 * t227;
t278 = qJD(3) * t221;
t265 = t146 * t218 + t151 * t278 + t183 * t293;
t41 = qJ(4) * t148 + t218 * t227 + t265;
t12 = -t212 * t41 + t214 * t36;
t6 = pkin(4) * t186 - pkin(7) * t82 + t12;
t13 = t212 * t36 + t214 * t41;
t7 = pkin(7) * t81 + t13;
t75 = pkin(3) * t188 + t88;
t42 = t214 * t75 - t83;
t25 = pkin(4) * t188 + t310 + t42;
t43 = t212 * t75 + t301;
t26 = t43 - t309;
t8 = -t217 * t26 + t220 * t25;
t2 = qJD(5) * t8 + t217 * t6 + t220 * t7;
t9 = t217 * t25 + t220 * t26;
t3 = -qJD(5) * t9 - t217 * t7 + t220 * t6;
t330 = -t3 * mrSges(6,1) + t2 * mrSges(6,2);
t201 = cos(t211);
t206 = t221 * pkin(3);
t171 = pkin(4) * t201 + t206;
t216 = -qJ(4) - pkin(6);
t242 = mrSges(3,1) * t215 - mrSges(3,2) * t213;
t329 = -m(4) * t245 - mrSges(2,1) - t242 + (-m(6) * (pkin(2) + t171) - m(5) * (t206 + pkin(2))) * t215 + (-m(6) * (pkin(7) - t216) - mrSges(6,3) - mrSges(4,3) + m(5) * t216 - mrSges(5,3)) * t213;
t234 = (-Ifges(4,2) * t218 + t304) * t213;
t235 = (Ifges(4,1) * t221 - t305) * t213;
t328 = t221 * (Ifges(4,6) * t188 + qJD(1) * t234) + t218 * (Ifges(4,5) * t188 + qJD(1) * t235);
t246 = t218 * t262;
t52 = -qJD(1) * t246 + t265;
t53 = -qJD(3) * t105 - t183 * t295 + t131;
t327 = -t53 * mrSges(4,1) - t12 * mrSges(5,1) + t52 * mrSges(4,2) + t13 * mrSges(5,2);
t223 = qJD(1) ^ 2;
t325 = -t250 / 0.2e1;
t324 = t64 / 0.2e1;
t323 = -t64 / 0.2e1;
t321 = t8 * mrSges(6,3);
t320 = t9 * mrSges(6,3);
t319 = -t123 / 0.2e1;
t180 = qJD(5) + t188;
t317 = -t180 / 0.2e1;
t314 = Ifges(6,4) * t64;
t308 = g(1) * t213;
t280 = qJD(2) * t215;
t285 = t174 * t278 + t221 * t280;
t72 = qJD(3) * t232 - t218 * t277 + t285;
t258 = t218 * t280;
t73 = -t258 - t221 * t277 + (-t187 + (qJ(4) * t213 - t174) * t218) * qJD(3);
t35 = t212 * t73 + t214 * t72;
t297 = t213 * t218;
t106 = -qJ(4) * t297 + t122;
t164 = t221 * t174;
t93 = -t266 + t164 + (-qJ(2) * t218 - pkin(3)) * t215;
t51 = t106 * t214 + t212 * t93;
t307 = mrSges(5,3) * t123;
t306 = mrSges(5,3) * t124;
t303 = Ifges(5,4) * t123;
t300 = qJ(2) * t223;
t202 = t213 * qJ(2);
t172 = t213 * t183;
t294 = t215 * t219;
t292 = t215 * t222;
t291 = t218 * t219;
t290 = t218 * t222;
t289 = t219 * t221;
t288 = t221 * t222;
t203 = qJ(5) + t211;
t194 = sin(t203);
t195 = cos(t203);
t117 = -t194 * t294 - t195 * t222;
t118 = -t194 * t222 + t195 * t294;
t287 = mrSges(6,1) * t117 - mrSges(6,2) * t118;
t119 = t194 * t292 - t195 * t219;
t120 = t194 * t219 + t195 * t292;
t286 = mrSges(6,1) * t119 + mrSges(6,2) * t120;
t284 = t355 * qJ(2);
t199 = t213 * qJD(2);
t257 = t213 * t278;
t161 = pkin(3) * t257 + t199;
t167 = pkin(3) * t297 + t202;
t279 = qJD(3) * t213;
t273 = qJDD(1) * t213;
t271 = t344 + t349 + t351;
t269 = mrSges(4,3) * t297;
t152 = pkin(3) * t260 + qJ(2) * t282 + qJD(4);
t256 = -mrSges(5,1) * t81 + mrSges(5,2) * t82;
t255 = -mrSges(6,1) * t22 + mrSges(6,2) * t21;
t34 = -t212 * t72 + t214 * t73;
t50 = -t106 * t212 + t214 * t93;
t249 = pkin(3) * t259;
t248 = mrSges(4,3) * t260;
t247 = mrSges(4,3) * t259;
t100 = -pkin(3) * t148 + qJDD(4) + t172;
t243 = -mrSges(3,1) * t272 + mrSges(3,2) * t273;
t241 = mrSges(4,1) * t218 + mrSges(4,2) * t221;
t240 = -mrSges(6,1) * t194 - mrSges(6,2) * t195;
t142 = t237 * t213;
t40 = -pkin(4) * t215 + pkin(7) * t142 + t50;
t141 = t166 * t213;
t44 = -pkin(7) * t141 + t51;
t14 = -t217 * t44 + t220 * t40;
t15 = t217 * t40 + t220 * t44;
t76 = -t141 * t220 + t142 * t217;
t77 = -t141 * t217 - t142 * t220;
t143 = -mrSges(4,2) * t188 - t248;
t144 = mrSges(4,1) * t188 - t247;
t238 = t143 * t221 - t144 * t218;
t236 = t271 - t330;
t157 = t215 * t290 - t289;
t155 = -t215 * t291 - t288;
t230 = t188 * t213 * (-Ifges(4,5) * t218 - Ifges(4,6) * t221);
t198 = -qJDD(1) * pkin(1) + qJDD(2);
t191 = t209 * t300;
t159 = t241 * t213;
t158 = t215 * t288 + t291;
t156 = t215 * t289 - t290;
t145 = t241 * t282;
t138 = t200 * t219 + t201 * t292;
t137 = t200 * t292 - t201 * t219;
t136 = -t200 * t222 + t201 * t294;
t135 = -t200 * t294 - t201 * t222;
t128 = t213 * t153;
t126 = t237 * t279;
t121 = t164 - t267;
t116 = Ifges(5,4) * t124;
t108 = -mrSges(4,2) * t186 + mrSges(4,3) * t148;
t107 = mrSges(4,1) * t186 - mrSges(4,3) * t147;
t104 = -t218 * t263 + t132;
t103 = -qJD(3) * t122 - t258;
t102 = -t246 + t285;
t101 = -pkin(4) * t123 + t249;
t98 = mrSges(5,1) * t188 + t307;
t97 = -mrSges(5,2) * t188 - t306;
t96 = pkin(4) * t141 + t167;
t90 = -pkin(4) * t126 + t161;
t87 = pkin(4) * t124 + t152;
t74 = mrSges(5,1) * t124 - mrSges(5,2) * t123;
t60 = Ifges(6,4) * t250;
t59 = mrSges(5,1) * t186 - mrSges(5,3) * t82;
t58 = -mrSges(5,2) * t186 + mrSges(5,3) * t81;
t57 = -t123 * Ifges(5,1) + t188 * Ifges(5,5) - t116;
t56 = -t124 * Ifges(5,2) + t188 * Ifges(5,6) - t303;
t55 = mrSges(6,1) * t180 + mrSges(6,3) * t64;
t54 = -mrSges(6,2) * t180 + mrSges(6,3) * t250;
t47 = -pkin(4) * t81 + t100;
t39 = -qJD(5) * t77 + t126 * t220 + t128 * t217;
t38 = qJD(5) * t76 + t126 * t217 - t128 * t220;
t33 = -mrSges(6,1) * t250 - mrSges(6,2) * t64;
t28 = -Ifges(6,1) * t64 + Ifges(6,5) * t180 + t60;
t27 = Ifges(6,2) * t250 + Ifges(6,6) * t180 - t314;
t24 = pkin(7) * t126 + t35;
t23 = pkin(7) * t128 + t34;
t17 = -mrSges(6,2) * t178 + mrSges(6,3) * t22;
t16 = mrSges(6,1) * t178 - mrSges(6,3) * t21;
t5 = -qJD(5) * t15 - t217 * t24 + t220 * t23;
t4 = qJD(5) * t14 + t217 * t23 + t220 * t24;
t1 = [-t198 * t242 - pkin(1) * t243 + (t274 * t337 + t333 + t355) * mrSges(3,3) - t328 * t279 / 0.2e1 + m(5) * (t100 * t167 + t12 * t50 + t13 * t51 + t152 * t161 + t34 * t42 + t35 * t43) + m(6) * (t14 * t3 + t15 * t2 + t4 * t9 + t47 * t96 + t5 * t8 + t87 * t90) + m(3) * (-pkin(1) * t198 + qJ(2) * t333 + t284) + t331 * t275 - (t100 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,1) * t82 + Ifges(5,4) * t81 + Ifges(5,5) * t335) * t142 - (-t100 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t82 + Ifges(5,2) * t81 + Ifges(5,6) * t335) * t141 - (t271 + t336) * t215 / 0.2e1 + (-t105 * t257 - t296 * t53) * mrSges(4,3) + (Ifges(3,4) * t272 + Ifges(3,1) * t273 + (Ifges(4,5) * t221 - Ifges(4,6) * t218) * t316) * t213 + t145 * t199 + (t126 * t43 + t128 * t42) * mrSges(5,3) + (-Ifges(5,1) * t128 + Ifges(5,4) * t126) * t319 - t124 * (-Ifges(5,4) * t128 + Ifges(5,2) * t126) / 0.2e1 + t152 * (-mrSges(5,1) * t126 - mrSges(5,2) * t128) + t188 * (-Ifges(5,5) * t128 + Ifges(5,6) * t126) / 0.2e1 + (-mrSges(6,1) * t47 + mrSges(6,3) * t2 + Ifges(6,4) * t21 + Ifges(6,2) * t22 + Ifges(6,6) * t178) * t76 + (-t291 * t326 - t158 * mrSges(4,1) - t138 * mrSges(5,1) - t120 * mrSges(6,1) + t157 * mrSges(4,2) + t137 * mrSges(5,2) + t119 * mrSges(6,2) - t354 * (pkin(1) * t222 + qJ(2) * t219) - t332 * t219 + t329 * t222) * g(2) + (t290 * t326 - t156 * mrSges(4,1) - t136 * mrSges(5,1) - t118 * mrSges(6,1) - t155 * mrSges(4,2) - t135 * mrSges(5,2) - t117 * mrSges(6,2) - t354 * (pkin(1) * t219 - qJ(2) * t222) + t332 * t222 + t329 * t219) * g(3) + (-mrSges(4,1) * t148 + mrSges(4,2) * t147) * t202 + t147 * t235 / 0.2e1 + t250 * (Ifges(6,4) * t38 + Ifges(6,2) * t39) / 0.2e1 - t52 * t269 + t39 * t320 + (Ifges(6,1) * t38 + Ifges(6,4) * t39) * t323 + (mrSges(6,2) * t47 - mrSges(6,3) * t3 + Ifges(6,1) * t21 + Ifges(6,4) * t22 + Ifges(6,5) * t178) * t77 + t14 * t16 + t15 * t17 + t38 * t28 / 0.2e1 + t39 * t27 / 0.2e1 + t96 * t255 + t167 * t256 + t4 * t54 + t5 * t55 + t51 * t58 + t50 * t59 + t87 * (-mrSges(6,1) * t39 + mrSges(6,2) * t38) + t90 * t33 + t35 * t97 + t34 * t98 + m(4) * (t102 * t105 + t103 * t104 + t121 * t53 + t122 * t52 + t284) + t121 * t107 + t122 * t108 + (t230 / 0.2e1 + t104 * t269) * qJD(3) + t126 * t56 / 0.2e1 + (Ifges(4,1) * t147 + Ifges(4,4) * t148 + Ifges(4,5) * t186) * t296 / 0.2e1 - t128 * t57 / 0.2e1 - (Ifges(4,4) * t147 + Ifges(4,2) * t148 + Ifges(4,6) * t186) * t297 / 0.2e1 + t159 * t172 + t102 * t143 + t103 * t144 + t161 * t74 + t180 * (Ifges(6,5) * t38 + Ifges(6,6) * t39) / 0.2e1 - t38 * t321 + t148 * t234 / 0.2e1 + Ifges(2,3) * qJDD(1) + (-t346 / 0.2e1 - t345 / 0.2e1 + Ifges(3,2) * t272 + Ifges(3,4) * t273 - t351 / 0.2e1 - t349 / 0.2e1 - t350 / 0.2e1 - t352 / 0.2e1 - t344 / 0.2e1 - t356 * t316 + t327 + t330) * t215; -t338 * t97 + t347 * t55 + t348 * t54 + t238 * qJD(3) + (-t238 * t215 + (-t145 - t33 - t74) * t213) * qJD(1) + t243 - t337 * t223 * mrSges(3,3) + t94 * t16 + t95 * t17 + t339 * t98 - t237 * t59 + t166 * t58 + t218 * t108 + t221 * t107 + (g(2) * t222 + g(3) * t219) * t354 + (t2 * t95 - t282 * t87 + t3 * t94 + t347 * t8 + t348 * t9) * m(6) + (-t12 * t237 + t13 * t166 - t152 * t282 - t338 * t43 + t339 * t42) * m(5) + (t52 * t218 + t53 * t221 - t191 + t188 * (-t104 * t218 + t105 * t221)) * m(4) + (-t210 * t300 - t191 + t198) * m(3); t328 * t282 / 0.2e1 + t341 * t54 - t331 * t223 - t327 + t336 - qJD(1) * t230 / 0.2e1 - t152 * (-mrSges(5,1) * t123 - mrSges(5,2) * t124) + t123 * (-Ifges(5,1) * t124 + t303) / 0.2e1 - t188 * (-Ifges(5,5) * t124 + Ifges(5,6) * t123) / 0.2e1 - t42 * t306 + t59 * t312 + t58 * t313 + t56 * t319 + t236 + (Ifges(5,2) * t123 - t116 + t57) * t124 / 0.2e1 + (Ifges(6,5) * t317 + t321 + Ifges(6,1) * t324 + Ifges(6,4) * t325 - t28 / 0.2e1 - t87 * mrSges(6,2)) * t250 + (t12 * t214 + t13 * t212) * t326 + (-t101 * t87 + t149 * t3 + t150 * t2 + t170 * t308 + t341 * t9 + t342 * t8) * m(6) + t342 * t55 + (-m(6) * (-t170 * t294 - t171 * t222) - t287 + mrSges(4,2) * t156 - t135 * mrSges(5,1) + t136 * mrSges(5,2) + t343 * t155) * g(2) + (-m(6) * (t170 * t292 - t171 * t219) - t286 - mrSges(4,2) * t158 - t137 * mrSges(5,1) - t138 * mrSges(5,2) + t343 * t157) * g(3) + (-t248 - t143) * t104 - t74 * t249 - m(5) * (t152 * t249 + t42 * t45 + t43 * t46) + (t247 + t144) * t105 - t46 * t97 - t45 * t98 - t101 * t33 + t149 * t16 + t150 * t17 + g(1) * t159 - t43 * t307 + (m(5) * t311 + mrSges(5,1) * t200 + mrSges(5,2) * t201 - t240) * t308 + (Ifges(6,6) * t317 + Ifges(6,4) * t324 + Ifges(6,2) * t325 - t27 / 0.2e1 + t87 * mrSges(6,1) - t320) * t64; t322 * t215 * g(1) - t123 * t98 + t124 * t97 - t250 * t54 - t64 * t55 + t255 + t256 + (-t250 * t9 - t64 * t8 + t340 + t47) * m(6) + (-t123 * t42 + t124 * t43 + t100 + t340) * m(5); -t87 * (-mrSges(6,1) * t64 + mrSges(6,2) * t250) + (Ifges(6,1) * t250 + t314) * t324 + t27 * t323 + (Ifges(6,5) * t250 + Ifges(6,6) * t64) * t317 - t8 * t54 + t9 * t55 - g(2) * t287 - g(3) * t286 - t240 * t308 + (t250 * t8 - t64 * t9) * mrSges(6,3) + t236 + (Ifges(6,2) * t64 + t28 + t60) * t325;];
tau = t1;
