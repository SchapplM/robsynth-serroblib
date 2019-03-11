% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:32
% EndTime: 2019-03-09 04:57:42
% DurationCPUTime: 5.14s
% Computational Cost: add. (9994->489), mult. (24480->678), div. (0->0), fcn. (17302->10), ass. (0->236)
t219 = cos(qJ(3));
t203 = sin(pkin(10)) * pkin(1) + pkin(7);
t195 = t203 * qJD(1);
t250 = pkin(8) * qJD(1) + t195;
t216 = sin(qJ(3));
t261 = t216 * qJD(2);
t157 = t219 * t250 + t261;
t215 = sin(qJ(4));
t151 = t215 * t157;
t209 = t219 * qJD(2);
t156 = -t250 * t216 + t209;
t218 = cos(qJ(4));
t111 = t218 * t156 - t151;
t210 = qJD(3) + qJD(4);
t214 = sin(qJ(6));
t217 = cos(qJ(6));
t192 = -t215 * t216 + t218 * t219;
t187 = t192 * qJD(1);
t193 = t215 * t219 + t218 * t216;
t188 = t193 * qJD(1);
t211 = sin(pkin(11));
t278 = cos(pkin(11));
t227 = t211 * t187 + t188 * t278;
t118 = t210 * t217 - t214 * t227;
t249 = t278 * t187 - t188 * t211;
t124 = qJD(6) - t249;
t244 = mrSges(7,1) * t214 + mrSges(7,2) * t217;
t153 = t218 * t157;
t299 = qJD(3) * pkin(3);
t154 = t156 + t299;
t107 = t154 * t215 + t153;
t275 = qJ(5) * t187;
t99 = t107 + t275;
t283 = t211 * t99;
t106 = t218 * t154 - t151;
t180 = t188 * qJ(5);
t98 = t106 - t180;
t90 = pkin(4) * t210 + t98;
t41 = t278 * t90 - t283;
t39 = -t210 * pkin(5) - t41;
t228 = t39 * t244;
t239 = Ifges(7,5) * t217 - Ifges(7,6) * t214;
t301 = Ifges(7,4) * t217;
t241 = -Ifges(7,2) * t214 + t301;
t302 = Ifges(7,4) * t214;
t243 = Ifges(7,1) * t217 - t302;
t316 = t217 / 0.2e1;
t317 = -t214 / 0.2e1;
t119 = t210 * t214 + t217 * t227;
t323 = t119 / 0.2e1;
t303 = Ifges(7,4) * t119;
t52 = Ifges(7,2) * t118 + Ifges(7,6) * t124 + t303;
t116 = Ifges(7,4) * t118;
t53 = t119 * Ifges(7,1) + t124 * Ifges(7,5) + t116;
t340 = t53 * t316 + t52 * t317 + t124 * t239 / 0.2e1 + t243 * t323 + t118 * t241 / 0.2e1 + t228;
t337 = mrSges(6,3) * t249;
t291 = t227 * Ifges(6,4);
t294 = t249 * Ifges(6,4);
t101 = -t180 + t111;
t110 = -t156 * t215 - t153;
t230 = t110 - t275;
t251 = t278 * t215;
t300 = pkin(3) * qJD(4);
t336 = t101 * t211 - t230 * t278 - (t211 * t218 + t251) * t300;
t268 = t211 * t215;
t174 = (t218 * t278 - t268) * t300;
t46 = t101 * t278 + t211 * t230;
t335 = -t46 + t174;
t279 = -mrSges(6,1) * t210 - mrSges(7,1) * t118 + mrSges(7,2) * t119 + mrSges(6,3) * t227;
t306 = pkin(8) + t203;
t190 = t306 * t216;
t191 = t306 * t219;
t136 = -t215 * t190 + t218 * t191;
t333 = -t195 * t216 + t209;
t149 = t210 * t192;
t137 = t149 * qJD(1);
t150 = t210 * t193;
t138 = t150 * qJD(1);
t94 = t137 * t278 - t211 * t138;
t56 = qJD(6) * t118 + t217 * t94;
t93 = t137 * t211 + t138 * t278;
t23 = mrSges(7,1) * t93 - mrSges(7,3) * t56;
t57 = -qJD(6) * t119 - t214 * t94;
t24 = -mrSges(7,2) * t93 + mrSges(7,3) * t57;
t332 = -t214 * t23 + t217 * t24;
t95 = t278 * t99;
t42 = t211 * t90 + t95;
t40 = pkin(9) * t210 + t42;
t258 = -cos(pkin(10)) * pkin(1) - pkin(2);
t194 = -pkin(3) * t219 + t258;
t189 = qJD(1) * t194;
t142 = -t187 * pkin(4) + qJD(5) + t189;
t67 = -pkin(5) * t249 - pkin(9) * t227 + t142;
t14 = -t214 * t40 + t217 * t67;
t15 = t214 * t67 + t217 * t40;
t331 = -t14 * t214 + t15 * t217;
t267 = qJD(1) * t216;
t207 = pkin(3) * t267;
t122 = pkin(4) * t138 + qJD(3) * t207;
t27 = pkin(5) * t93 - pkin(9) * t94 + t122;
t170 = t195 * t219 + t261;
t266 = qJD(1) * t219;
t221 = (-t215 * (-pkin(8) * t267 + t333) + t218 * (-pkin(8) * t266 - t170)) * qJD(3);
t233 = -t137 * qJ(5) - t188 * qJD(5);
t264 = qJD(4) * t218;
t265 = qJD(4) * t215;
t63 = qJD(3) * t111 + t154 * t264 - t157 * t265;
t36 = -qJ(5) * t138 + qJD(5) * t187 + t63;
t9 = t278 * t36 + (-t154 * t265 - t157 * t264 + t221 + t233) * t211;
t2 = qJD(6) * t14 + t214 * t27 + t217 * t9;
t274 = qJD(6) * t15;
t3 = -t214 * t9 + t217 * t27 - t274;
t330 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t56 + Ifges(7,6) * t57;
t329 = t56 / 0.2e1;
t328 = t57 / 0.2e1;
t327 = t93 / 0.2e1;
t117 = qJ(5) * t192 + t136;
t135 = -t218 * t190 - t191 * t215;
t229 = -qJ(5) * t193 + t135;
t68 = t117 * t211 - t229 * t278;
t64 = -qJD(4) * t107 + t221;
t8 = t211 * t36 - t278 * (t233 + t64);
t326 = t68 * t8;
t325 = -t118 / 0.2e1;
t324 = -t119 / 0.2e1;
t322 = -t124 / 0.2e1;
t320 = t187 / 0.2e1;
t319 = t188 / 0.2e1;
t315 = m(5) * t189;
t314 = pkin(4) * t188;
t313 = pkin(4) * t211;
t312 = t14 * mrSges(7,3);
t143 = -t192 * t278 + t193 * t211;
t311 = t143 * t8;
t310 = t2 * t217;
t309 = t214 * t3;
t308 = t93 * mrSges(6,3);
t307 = t94 * mrSges(6,3);
t305 = mrSges(5,3) * t187;
t304 = Ifges(4,4) * t216;
t298 = t118 * Ifges(7,6);
t297 = t119 * Ifges(7,5);
t296 = t124 * Ifges(7,3);
t295 = t227 * t42;
t293 = t249 * Ifges(6,2);
t292 = t227 * Ifges(6,1);
t288 = t188 * mrSges(5,3);
t287 = t188 * Ifges(5,4);
t197 = t258 * qJD(1);
t286 = t197 * mrSges(4,2);
t285 = t210 * Ifges(6,5);
t284 = t210 * Ifges(6,6);
t277 = Ifges(4,5) * qJD(3);
t276 = Ifges(4,6) * qJD(3);
t273 = t249 * t214;
t272 = t249 * t217;
t271 = t192 * t137;
t270 = t193 * t138;
t205 = pkin(3) * t218 + pkin(4);
t179 = pkin(3) * t251 + t211 * t205;
t263 = qJD(6) * t214;
t262 = qJD(6) * t217;
t206 = Ifges(4,4) * t266;
t257 = t278 * pkin(4);
t256 = t277 / 0.2e1;
t255 = -t276 / 0.2e1;
t254 = m(4) * t203 + mrSges(4,3);
t253 = t93 * mrSges(6,1) + t94 * mrSges(6,2);
t132 = pkin(4) * t150 + t216 * t299;
t252 = qJD(3) * t306;
t248 = t276 / 0.2e1 + (t219 * Ifges(4,2) + t304) * qJD(1) / 0.2e1 - t197 * mrSges(4,1);
t246 = -t2 * t214 - t217 * t3;
t245 = mrSges(7,1) * t217 - mrSges(7,2) * t214;
t242 = Ifges(7,1) * t214 + t301;
t240 = Ifges(7,2) * t217 + t302;
t238 = Ifges(7,5) * t214 + Ifges(7,6) * t217;
t237 = -t14 * t217 - t15 * t214;
t69 = t117 * t278 + t211 * t229;
t144 = t211 * t192 + t193 * t278;
t155 = -t192 * pkin(4) + t194;
t73 = t143 * pkin(5) - t144 * pkin(9) + t155;
t26 = t214 * t73 + t217 * t69;
t25 = -t214 * t69 + t217 * t73;
t74 = -mrSges(7,2) * t124 + mrSges(7,3) * t118;
t75 = mrSges(7,1) * t124 - mrSges(7,3) * t119;
t235 = -t214 * t75 + t217 * t74;
t120 = -mrSges(6,2) * t210 + t337;
t231 = -t120 - t235;
t181 = t216 * t252;
t182 = t219 * t252;
t91 = -t218 * t181 - t215 * t182 - t190 * t264 - t191 * t265;
t72 = pkin(5) * t227 - pkin(9) * t249 + t314;
t178 = -pkin(3) * t268 + t205 * t278;
t92 = -qJD(4) * t136 + t181 * t215 - t218 * t182;
t223 = qJD(6) * t237 - t309 + t310;
t222 = -qJ(5) * t149 - qJD(5) * t193 + t92;
t12 = t56 * Ifges(7,4) + t57 * Ifges(7,2) + t93 * Ifges(7,6);
t125 = t187 * Ifges(5,2) + t210 * Ifges(5,6) + t287;
t183 = Ifges(5,4) * t187;
t126 = t188 * Ifges(5,1) + t210 * Ifges(5,5) + t183;
t13 = t56 * Ifges(7,1) + t57 * Ifges(7,4) + t93 * Ifges(7,5);
t51 = t296 + t297 + t298;
t78 = t284 + t291 + t293;
t79 = t285 + t292 + t294;
t220 = -(Ifges(5,5) * t187 + Ifges(6,5) * t249 - Ifges(5,6) * t188 - Ifges(6,6) * t227) * t210 / 0.2e1 - (-Ifges(5,2) * t188 + t126 + t183) * t187 / 0.2e1 - (-Ifges(6,2) * t227 + t294 + t79) * t249 / 0.2e1 - (Ifges(6,1) * t249 - t291 + t51) * t227 / 0.2e1 + (-t245 - mrSges(6,1)) * t8 + t340 * qJD(6) - t53 * t272 / 0.2e1 + t52 * t273 / 0.2e1 + t238 * t327 + t240 * t328 + t242 * t329 + t12 * t316 + t125 * t319 + t106 * t305 + mrSges(7,3) * t310 + t107 * t288 + t41 * t337 - t188 * (Ifges(5,1) * t187 - t287) / 0.2e1 - t15 * (-mrSges(7,2) * t227 - mrSges(7,3) * t273) - t14 * (mrSges(7,1) * t227 - mrSges(7,3) * t272) + t227 * t78 / 0.2e1 + (Ifges(7,3) * t227 + t239 * t249) * t322 - t142 * (mrSges(6,1) * t227 + mrSges(6,2) * t249) + (Ifges(7,5) * t227 + t243 * t249) * t324 + (Ifges(7,6) * t227 + t241 * t249) * t325 - t189 * (mrSges(5,1) * t188 + mrSges(5,2) * t187) - t249 * t228 - t9 * mrSges(6,2) - t63 * mrSges(5,2) + t64 * mrSges(5,1) - Ifges(6,6) * t93 + Ifges(6,5) * t94 + Ifges(5,5) * t137 - Ifges(5,6) * t138 + t214 * t13 / 0.2e1;
t204 = -t257 - pkin(5);
t198 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t266;
t196 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t267;
t186 = Ifges(4,1) * t267 + t206 + t277;
t172 = pkin(9) + t179;
t171 = -pkin(5) - t178;
t162 = t170 * qJD(3);
t161 = t333 * qJD(3);
t160 = mrSges(5,1) * t210 - t288;
t159 = -mrSges(5,2) * t210 + t305;
t158 = t207 + t314;
t140 = -mrSges(5,1) * t187 + mrSges(5,2) * t188;
t103 = t149 * t278 - t211 * t150;
t102 = t149 * t211 + t150 * t278;
t87 = Ifges(7,3) * t93;
t85 = -mrSges(6,1) * t249 + mrSges(6,2) * t227;
t70 = t207 + t72;
t61 = -qJ(5) * t150 + qJD(5) * t192 + t91;
t44 = t278 * t98 - t283;
t43 = t211 * t98 + t95;
t37 = pkin(5) * t102 - pkin(9) * t103 + t132;
t22 = t211 * t222 + t278 * t61;
t21 = t211 * t61 - t222 * t278;
t20 = t214 * t70 + t217 * t46;
t19 = -t214 * t46 + t217 * t70;
t18 = t214 * t72 + t217 * t44;
t17 = -t214 * t44 + t217 * t72;
t16 = -mrSges(7,1) * t57 + mrSges(7,2) * t56;
t5 = -qJD(6) * t26 - t214 * t22 + t217 * t37;
t4 = qJD(6) * t25 + t214 * t37 + t217 * t22;
t1 = [(t254 * t161 + (-t203 * t196 + t186 / 0.2e1 + 0.3e1 / 0.2e1 * t206 + t256 - t254 * t333 + 0.2e1 * t286) * qJD(3)) * t219 + t189 * (mrSges(5,1) * t150 + mrSges(5,2) * t149) + t210 * (Ifges(5,5) * t149 - Ifges(5,6) * t150) / 0.2e1 + (t12 * t317 + t13 * t316 + Ifges(6,1) * t94 - Ifges(6,4) * t93 + t122 * mrSges(6,2) + t243 * t329 + t241 * t328 + t239 * t327 + (mrSges(6,3) + t244) * t8 + t246 * mrSges(7,3) + (t53 * t317 - t217 * t52 / 0.2e1 + t39 * t245 + t240 * t325 + t242 * t324 + t238 * t322 - t331 * mrSges(7,3)) * qJD(6)) * t144 + (-t106 * t149 - t107 * t150 - t135 * t137 - t136 * t138 + t192 * t63 - t193 * t64) * mrSges(5,3) + (-t192 * t138 - t150 * t320) * Ifges(5,2) + t194 * (mrSges(5,1) * t138 + mrSges(5,2) * t137) + (t149 * t320 - t150 * t319 - t270 + t271) * Ifges(5,4) + (t79 / 0.2e1 + t294 / 0.2e1 + t292 / 0.2e1 + t142 * mrSges(6,2) + t285 / 0.2e1 + t237 * mrSges(7,3) + t340) * t103 + m(5) * (t106 * t92 + t107 * t91 + t135 * t64 + t136 * t63) + (-t102 * t42 - t103 * t41 + t68 * t94 - t69 * t93) * mrSges(6,3) + (t193 * t137 + t149 * t319) * Ifges(5,1) + (t254 * t162 + (-t203 * t198 + t255 - t254 * t170 + (t258 * mrSges(4,1) - 0.3e1 / 0.2e1 * t304 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t219) * qJD(1) + (qJD(1) * (-mrSges(5,1) * t192 + mrSges(5,2) * t193) + t140 + 0.2e1 * t315) * pkin(3) - t248) * qJD(3)) * t216 + (t14 * mrSges(7,1) - t15 * mrSges(7,2) + t298 / 0.2e1 + t297 / 0.2e1 + t296 / 0.2e1 + t51 / 0.2e1 - t78 / 0.2e1 - t293 / 0.2e1 - t291 / 0.2e1 + t142 * mrSges(6,1) - t284 / 0.2e1) * t102 + t279 * t21 + t155 * t253 + (-t9 * mrSges(6,3) + t87 / 0.2e1 - Ifges(6,4) * t94 + t122 * mrSges(6,1) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t93 + t330) * t143 + m(7) * (t14 * t5 + t15 * t4 + t2 * t26 + t21 * t39 + t25 * t3 + t326) + m(6) * (t122 * t155 + t132 * t142 - t21 * t41 + t22 * t42 + t69 * t9 + t326) + t25 * t23 + t26 * t24 + t68 * t16 + t4 * t74 + t5 * t75 + t22 * t120 + t132 * t85 + t149 * t126 / 0.2e1 - t150 * t125 / 0.2e1 + t91 * t159 + t92 * t160; t149 * t159 - t150 * t160 + (t16 + t307) * t143 + t279 * t102 + (-t270 - t271) * mrSges(5,3) - t231 * t103 + (-t216 * t196 + t219 * t198 + (-t216 ^ 2 - t219 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3) + (-t308 + (-t214 * t74 - t217 * t75) * qJD(6) + t332) * t144 + m(4) * (t161 * t216 - t162 * t219 + (t170 * t219 - t216 * t333) * qJD(3)) + m(5) * (-t106 * t150 + t107 * t149 + t192 * t64 + t193 * t63) + m(6) * (-t102 * t41 + t103 * t42 + t144 * t9 + t311) + m(7) * (t102 * t39 + t103 * t331 + t144 * t223 + t311); (t172 * t24 + t174 * t74 + (-t172 * t75 - t312) * qJD(6)) * t217 + (-t174 * t75 + (-qJD(6) * t74 - t23) * t172 + (-t3 - t274) * mrSges(7,3)) * t214 + t220 + ((t256 - t286 + t333 * mrSges(4,3) - t186 / 0.2e1 - t206 / 0.2e1) * t219 + (t255 + t170 * mrSges(4,3) + (t304 / 0.2e1 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t219) * qJD(1) + (-t140 - t315) * pkin(3) + t248) * t216) * qJD(1) + (-t178 * t94 - t179 * t93 + t295) * mrSges(6,3) + ((t159 * t218 - t160 * t215) * qJD(4) + (-t137 * t218 - t138 * t215) * mrSges(5,3)) * pkin(3) + t335 * t120 - t20 * t74 - t19 * t75 - t158 * t85 - t111 * t159 - t110 * t160 - t161 * mrSges(4,2) - t162 * mrSges(4,1) + t171 * t16 + t170 * t196 - t333 * t198 - t279 * t336 + (-t142 * t158 - t178 * t8 + t179 * t9 + t335 * t42 + t336 * t41) * m(6) + ((t215 * t63 + t218 * t64 + (-t106 * t215 + t107 * t218) * qJD(4)) * pkin(3) - t106 * t110 - t107 * t111) * m(5) + (-t14 * t19 - t15 * t20 + t171 * t8 + t172 * t223 + t174 * t331 - t336 * t39) * m(7); mrSges(6,3) * t295 - t106 * t159 + t107 * t160 - t44 * t120 + t204 * t16 - t17 * t75 - t18 * t74 - t257 * t307 - t262 * t312 - t308 * t313 - t85 * t314 + t220 - t279 * t43 + (-t15 * t263 - t309) * mrSges(7,3) + (-t14 * t17 - t15 * t18 + t204 * t8 - t39 * t43) * m(7) + ((t211 * t9 - t278 * t8) * pkin(4) - t142 * t314 + t41 * t43 - t42 * t44) * m(6) + (m(7) * t223 - t262 * t75 - t263 * t74 + t332) * (pkin(9) + t313); t235 * qJD(6) - t279 * t227 + t231 * t249 + t214 * t24 + t217 * t23 + t253 + (t124 * t331 - t227 * t39 - t246) * m(7) + (t227 * t41 - t249 * t42 + t122) * m(6); t87 - t39 * (mrSges(7,1) * t119 + mrSges(7,2) * t118) + (Ifges(7,1) * t118 - t303) * t324 + t52 * t323 + (Ifges(7,5) * t118 - Ifges(7,6) * t119) * t322 - t14 * t74 + t15 * t75 + (t118 * t14 + t119 * t15) * mrSges(7,3) + (-Ifges(7,2) * t119 + t116 + t53) * t325 + t330;];
tauc  = t1(:);
