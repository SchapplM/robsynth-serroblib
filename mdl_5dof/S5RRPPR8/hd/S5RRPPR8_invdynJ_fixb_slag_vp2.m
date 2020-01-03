% Calculate vector of inverse dynamics joint torques for
% S5RRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:08
% EndTime: 2019-12-31 19:38:30
% DurationCPUTime: 14.04s
% Computational Cost: add. (3326->506), mult. (7270->663), div. (0->0), fcn. (4677->10), ass. (0->231)
t329 = Ifges(4,4) + Ifges(3,5);
t328 = Ifges(4,6) - Ifges(3,6);
t330 = qJD(2) - qJD(5);
t197 = sin(pkin(8));
t198 = cos(pkin(8));
t201 = sin(qJ(5));
t204 = cos(qJ(5));
t130 = -t197 * t201 + t198 * t204;
t327 = t330 * t130;
t133 = t197 * t204 + t198 * t201;
t326 = t330 * t133;
t205 = cos(qJ(2));
t261 = t205 * qJD(1);
t175 = Ifges(3,4) * t261;
t202 = sin(qJ(2));
t282 = Ifges(4,5) * t205;
t227 = t202 * Ifges(4,1) - t282;
t262 = t202 * qJD(1);
t325 = Ifges(3,1) * t262 + qJD(1) * t227 + qJD(2) * t329 + t175;
t252 = mrSges(4,2) * t262;
t324 = -mrSges(3,3) * t262 - t252 + (mrSges(3,1) + mrSges(4,1)) * qJD(2);
t251 = mrSges(4,2) * t261;
t154 = qJD(2) * mrSges(4,3) + t251;
t323 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t261 + t154;
t322 = t328 * t202 + t205 * t329;
t229 = t205 * mrSges(4,1) + t202 * mrSges(4,3);
t231 = t205 * mrSges(3,1) - t202 * mrSges(3,2);
t321 = t229 + t231;
t259 = qJD(1) * qJD(2);
t144 = qJDD(1) * t202 + t205 * t259;
t257 = t205 * qJDD(1);
t173 = pkin(6) * t257;
t246 = t202 * t259;
t127 = -pkin(6) * t246 + t173;
t128 = t144 * pkin(6);
t320 = t127 * t205 + t128 * t202;
t313 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t89 = t127 + t313;
t241 = qJDD(3) + t128;
t97 = -qJDD(2) * pkin(2) + t241;
t319 = t202 * t97 + t205 * t89;
t260 = m(4) + m(5) + m(6);
t318 = mrSges(2,1) + t321;
t317 = m(5) * qJ(4) - m(6) * (-pkin(7) - qJ(4)) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + mrSges(6,3);
t107 = -t197 * t262 - t198 * t261;
t294 = pkin(7) * t107;
t176 = pkin(6) * t262;
t139 = qJ(4) * t262 - t176;
t177 = pkin(6) * t261;
t142 = -qJ(4) * t261 + t177;
t67 = -t139 * t197 + t198 * t142;
t31 = t67 + t294;
t108 = t197 * t261 - t198 * t262;
t293 = pkin(7) * t108;
t68 = t198 * t139 + t197 * t142;
t32 = t68 - t293;
t207 = -pkin(2) - pkin(3);
t145 = -qJ(3) * t197 + t198 * t207;
t138 = -pkin(4) + t145;
t146 = t198 * qJ(3) + t197 * t207;
t69 = t138 * t204 - t146 * t201;
t316 = qJD(3) * t130 + qJD(5) * t69 - t201 * t31 - t204 * t32;
t70 = t138 * t201 + t146 * t204;
t315 = -qJD(3) * t133 - qJD(5) * t70 + t201 * t32 - t204 * t31;
t314 = t204 * t107 + t108 * t201;
t203 = sin(qJ(1));
t206 = cos(qJ(1));
t312 = g(1) * t206 + g(2) * t203;
t288 = mrSges(5,3) * t107;
t84 = qJD(2) * mrSges(5,2) + t288;
t287 = mrSges(5,3) * t108;
t85 = -qJD(2) * mrSges(5,1) + t287;
t311 = t197 * t85 - t198 * t84 - t154;
t310 = qJ(3) * t260;
t253 = t207 * t202;
t271 = t197 * t205;
t278 = t205 * mrSges(4,3);
t166 = pkin(4) * t198 + pkin(3);
t290 = -pkin(2) - t166;
t309 = -m(5) * t253 - m(6) * (pkin(4) * t271 + t202 * t290) - t278 - (-m(4) * pkin(2) - mrSges(4,1)) * t202;
t307 = -t314 / 0.2e1;
t221 = t107 * t201 - t108 * t204;
t306 = -t221 / 0.2e1;
t305 = t221 / 0.2e1;
t195 = qJD(2) * qJ(3);
t118 = t142 + t195;
t249 = t207 * qJD(2);
t90 = qJD(3) + t249 - t139;
t40 = -t118 * t197 + t198 * t90;
t26 = -qJD(2) * pkin(4) + t293 + t40;
t41 = t198 * t118 + t197 * t90;
t27 = t41 + t294;
t5 = -t201 * t27 + t204 * t26;
t304 = t5 * mrSges(6,3);
t6 = t201 * t26 + t204 * t27;
t303 = t6 * mrSges(6,3);
t302 = t107 / 0.2e1;
t299 = t330 / 0.2e1;
t298 = t202 / 0.2e1;
t297 = Ifges(6,4) * t221;
t296 = pkin(6) * t202;
t295 = pkin(6) * t205;
t187 = t205 * pkin(2);
t289 = pkin(6) - qJ(4);
t47 = -qJ(4) * t144 - qJD(4) * t262 + qJDD(2) * t207 + t241;
t143 = t246 - t257;
t265 = qJD(2) * t202;
t255 = pkin(6) * t265;
t263 = qJD(4) * t205;
t48 = qJ(4) * t143 + t173 + (-t255 - t263) * qJD(1) + t313;
t19 = t197 * t47 + t198 * t48;
t104 = -t265 * t289 - t263;
t157 = t289 * t205;
t106 = qJD(2) * t157 - qJD(4) * t202;
t43 = t198 * t104 + t197 * t106;
t286 = Ifges(3,4) * t202;
t285 = Ifges(3,4) * t205;
t284 = Ifges(5,4) * t107;
t283 = Ifges(4,5) * t202;
t196 = qJDD(1) * pkin(1);
t272 = t197 * t202;
t184 = t202 * qJ(3);
t270 = t202 * t206;
t269 = t205 * t206;
t155 = t289 * t202;
t75 = t197 * t155 + t198 * t157;
t183 = t202 * qJD(3);
t264 = qJD(2) * t205;
t268 = qJ(3) * t264 + t183;
t267 = t187 + t184;
t266 = t206 * pkin(1) + t203 * pkin(6);
t123 = -qJD(1) * pkin(1) - pkin(2) * t261 - qJ(3) * t262;
t250 = t205 * pkin(3) + t267;
t71 = t143 * t198 - t144 * t197;
t72 = t143 * t197 + t144 * t198;
t248 = -t71 * mrSges(5,1) + t72 * mrSges(5,2);
t12 = qJD(5) * t314 + t201 * t71 + t204 * t72;
t13 = -qJD(5) * t221 - t201 * t72 + t204 * t71;
t247 = -t13 * mrSges(6,1) + t12 * mrSges(6,2);
t240 = -pkin(1) - t184;
t18 = -t197 * t48 + t198 * t47;
t239 = -t259 / 0.2e1;
t216 = t198 * t205 + t272;
t217 = -t198 * t202 + t271;
t237 = mrSges(5,1) * t216 - mrSges(5,2) * t217;
t42 = -t104 * t197 + t198 * t106;
t74 = t198 * t155 - t157 * t197;
t125 = pkin(1) + t250;
t100 = -qJDD(2) * mrSges(4,1) + t144 * mrSges(4,2);
t61 = t143 * pkin(2) - t144 * qJ(3) - qJD(1) * t183 - t196;
t191 = pkin(8) + qJ(5);
t180 = sin(t191);
t181 = cos(t191);
t219 = t180 * t205 - t181 * t202;
t79 = t219 * t203;
t218 = t180 * t202 + t181 * t205;
t80 = t218 * t203;
t234 = -t79 * mrSges(6,1) - t80 * mrSges(6,2);
t81 = t180 * t269 - t181 * t270;
t82 = t218 * t206;
t233 = -t81 * mrSges(6,1) - t82 * mrSges(6,2);
t232 = g(1) * t203 - g(2) * t206;
t230 = mrSges(3,1) * t202 + mrSges(3,2) * t205;
t228 = -mrSges(6,1) * t218 + mrSges(6,2) * t219;
t226 = t205 * Ifges(3,2) + t286;
t223 = -t197 * t40 + t198 * t41;
t44 = pkin(7) * t217 + t74;
t45 = -pkin(7) * t216 + t75;
t16 = -t201 * t45 + t204 * t44;
t17 = t201 * t44 + t204 * t45;
t83 = pkin(3) * t261 + qJD(4) - t123;
t64 = t201 * t217 - t204 * t216;
t65 = -t201 * t216 - t204 * t217;
t148 = -qJD(2) * pkin(2) + qJD(3) + t176;
t152 = t177 + t195;
t220 = t148 * t205 - t152 * t202;
t170 = qJ(3) * t261;
t98 = qJD(1) * t253 + t170;
t215 = pkin(4) * t272 + t166 * t205;
t214 = pkin(1) * t230;
t213 = t123 * (t202 * mrSges(4,1) - t278);
t212 = t202 * (Ifges(3,1) * t205 - t286);
t211 = t205 * (Ifges(4,3) * t202 + t282);
t78 = t202 * t249 + t268;
t30 = -pkin(3) * t143 + qJDD(4) - t61;
t7 = -qJDD(2) * pkin(4) - pkin(7) * t72 + t18;
t8 = pkin(7) * t71 + t19;
t1 = qJD(5) * t5 + t201 * t7 + t204 * t8;
t190 = -qJDD(2) + qJDD(5);
t2 = -qJD(5) * t6 - t201 * t8 + t204 * t7;
t208 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t12 + Ifges(6,6) * t13 + Ifges(6,3) * t190;
t174 = Ifges(4,5) * t262;
t149 = -pkin(1) - t267;
t141 = pkin(2) * t262 - t170;
t140 = t229 * qJD(1);
t120 = Ifges(3,6) * qJD(2) + qJD(1) * t226;
t119 = Ifges(4,6) * qJD(2) - Ifges(4,3) * t261 + t174;
t114 = t216 * qJD(2);
t113 = t217 * qJD(2);
t105 = pkin(2) * t265 - t268;
t101 = -mrSges(4,2) * t143 + qJDD(2) * mrSges(4,3);
t99 = Ifges(5,4) * t108;
t96 = t216 * t206;
t95 = t217 * t206;
t94 = t216 * t203;
t93 = t217 * t203;
t73 = pkin(4) * t216 + t125;
t63 = pkin(4) * t108 + t98;
t62 = -mrSges(5,1) * t107 - mrSges(5,2) * t108;
t60 = -qJDD(2) * mrSges(5,1) - mrSges(5,3) * t72;
t59 = qJDD(2) * mrSges(5,2) + mrSges(5,3) * t71;
t52 = -t108 * Ifges(5,1) - Ifges(5,5) * qJD(2) + t284;
t51 = t107 * Ifges(5,2) - Ifges(5,6) * qJD(2) - t99;
t50 = -pkin(4) * t107 + t83;
t49 = pkin(4) * t113 + t78;
t46 = Ifges(6,4) * t314;
t36 = -mrSges(6,1) * t330 - mrSges(6,3) * t221;
t35 = mrSges(6,2) * t330 + mrSges(6,3) * t314;
t29 = -pkin(7) * t113 + t43;
t28 = -pkin(7) * t114 + t42;
t25 = -qJD(5) * t65 - t113 * t204 - t114 * t201;
t24 = qJD(5) * t64 - t113 * t201 + t114 * t204;
t23 = -mrSges(6,1) * t314 + mrSges(6,2) * t221;
t22 = -pkin(4) * t71 + t30;
t21 = Ifges(6,1) * t221 - Ifges(6,5) * t330 + t46;
t20 = Ifges(6,2) * t314 - Ifges(6,6) * t330 + t297;
t10 = -mrSges(6,2) * t190 + mrSges(6,3) * t13;
t9 = mrSges(6,1) * t190 - mrSges(6,3) * t12;
t4 = -qJD(5) * t17 - t201 * t29 + t204 * t28;
t3 = qJD(5) * t16 + t201 * t28 + t204 * t29;
t11 = [(mrSges(6,2) * t22 - mrSges(6,3) * t2 + Ifges(6,1) * t12 + Ifges(6,4) * t13 + Ifges(6,5) * t190) * t65 + t205 * (Ifges(3,4) * t144 + Ifges(3,6) * qJDD(2)) / 0.2e1 + (-mrSges(6,1) * t22 + mrSges(6,3) * t1 + Ifges(6,4) * t12 + Ifges(6,2) * t13 + Ifges(6,6) * t190) * t64 - t105 * t140 - t113 * t51 / 0.2e1 + t114 * t52 / 0.2e1 + t78 * t62 + t43 * t84 + t42 * t85 + t74 * t60 + t75 * t59 + (-t152 * t265 + t319) * mrSges(4,2) + (-t324 * pkin(6) + t148 * mrSges(4,2) + t325 / 0.2e1) * t264 + (t94 * mrSges(5,1) + t80 * mrSges(6,1) - t93 * mrSges(5,2) - t79 * mrSges(6,2) + (m(3) * pkin(1) - m(6) * (-pkin(1) + t290 * t205 + (-pkin(4) * t197 - qJ(3)) * t202) - m(4) * (t240 - t187) - m(5) * (t205 * t207 + t240) + t318) * t203 + ((-m(3) - t260) * pkin(6) + t317) * t206) * g(1) - t214 * t259 + t49 * t23 + t50 * (-mrSges(6,1) * t25 + mrSges(6,2) * t24) + t3 * t35 + t4 * t36 + t25 * t20 / 0.2e1 + (-qJDD(2) * mrSges(3,1) + t100) * t296 + t24 * t21 / 0.2e1 + t16 * t9 + t17 * t10 + t211 * t239 - t149 * mrSges(4,3) * t144 + (-m(5) * pkin(3) * t269 - m(3) * t266 - t96 * mrSges(5,1) - t82 * mrSges(6,1) + t95 * mrSges(5,2) + t81 * mrSges(6,2) - t260 * (pkin(2) * t269 + qJ(3) * t270 + t266) + (-m(6) * t215 - t318) * t206 + t317 * t203) * g(2) + (t212 + t202 * (Ifges(4,1) * t205 + t283) + t205 * (-Ifges(3,2) * t202 + t285)) * t259 / 0.2e1 + (t202 * Ifges(3,1) + t227 + t285) * t144 / 0.2e1 + (-t113 * t41 - t114 * t40) * mrSges(5,3) + m(5) * (t125 * t30 + t18 * t74 + t19 * t75 + t40 * t42 + t41 * t43 + t78 * t83) + m(6) * (t1 * t17 + t16 * t2 + t22 * t73 + t3 * t6 + t4 * t5 + t49 * t50) - t61 * t229 - t108 * (Ifges(5,1) * t114 - Ifges(5,4) * t113) / 0.2e1 + t83 * (mrSges(5,1) * t113 + mrSges(5,2) * t114) - qJD(2) * (Ifges(5,5) * t114 - Ifges(5,6) * t113) / 0.2e1 + (Ifges(5,4) * t114 - Ifges(5,2) * t113) * t302 - t205 * (Ifges(4,5) * t144 + Ifges(4,6) * qJDD(2)) / 0.2e1 + (-t120 / 0.2e1 + t119 / 0.2e1) * t265 + (-mrSges(5,3) * t19 - Ifges(5,4) * t72 - Ifges(5,2) * t71 + Ifges(5,6) * qJDD(2)) * t216 + (mrSges(5,3) * t18 - Ifges(5,1) * t72 - Ifges(5,4) * t71 + Ifges(5,5) * qJDD(2)) * t217 + t231 * t196 + t314 * (Ifges(6,4) * t24 + Ifges(6,2) * t25) / 0.2e1 + (-pkin(1) * mrSges(3,1) + t149 * mrSges(4,1) - t226 / 0.2e1 + t283 / 0.2e1 - mrSges(3,3) * t295 + (Ifges(4,5) - Ifges(3,4)) * t298 + (-Ifges(4,3) - Ifges(3,2) / 0.2e1) * t205) * t143 + t73 * t247 + t30 * t237 - t330 * (Ifges(6,5) * t24 + Ifges(6,6) * t25) / 0.2e1 + ((Ifges(3,1) + Ifges(4,1)) * t144 + t329 * qJDD(2)) * t298 + (t202 * t329 - t328 * t205) * qJDD(2) / 0.2e1 + t125 * t248 + (-pkin(1) * t144 - qJDD(2) * t295) * mrSges(3,2) + m(4) * (t105 * t123 + t149 * t61 + (qJD(2) * t220 + t319) * pkin(6)) + (t144 * t296 + t320) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t320) + t322 * qJD(2) ^ 2 / 0.2e1 - t323 * t255 - t24 * t304 + t101 * t295 + qJD(2) * t213 + t25 * t303 + (Ifges(6,1) * t24 + Ifges(6,4) * t25) * t305 + Ifges(2,3) * qJDD(1); (-pkin(2) * t97 + qJ(3) * t89 - t123 * t141) * m(4) + (-m(4) * t220 * pkin(6) - t213 + (-t212 / 0.2e1 + t211 / 0.2e1 + t214) * qJD(1)) * qJD(1) + t145 * t60 + t146 * t59 + t141 * t140 - t127 * mrSges(3,2) - t128 * mrSges(3,1) - t83 * (mrSges(5,1) * t108 - mrSges(5,2) * t107) + qJD(2) * (-Ifges(5,5) * t107 - Ifges(5,6) * t108) / 0.2e1 - t97 * mrSges(4,1) - t98 * t62 - pkin(2) * t100 + qJ(3) * t101 + t89 * mrSges(4,3) - t68 * t84 - t67 * t85 + t69 * t9 + t70 * t10 - Ifges(5,6) * t71 - Ifges(5,5) * t72 + (t145 * t18 + t146 * t19 - t40 * t67 - t41 * t68 - t83 * t98) * m(5) - t63 * t23 - t18 * mrSges(5,1) + t19 * mrSges(5,2) + (-t20 / 0.2e1 + t50 * mrSges(6,1) - t303 + Ifges(6,6) * t299 + Ifges(6,4) * t306 + Ifges(6,2) * t307) * t221 - t148 * t251 - (Ifges(4,1) * t261 + t119 + t174) * t262 / 0.2e1 + (-Ifges(5,1) * t107 + t51 - t99) * t108 / 0.2e1 + (Ifges(5,3) + Ifges(4,2) + Ifges(3,3)) * qJDD(2) - (-t50 * mrSges(6,2) - t21 / 0.2e1 + Ifges(6,5) * t299 + t304 + Ifges(6,1) * t306 + Ifges(6,4) * t307) * t314 + (-t95 * mrSges(5,1) - t96 * mrSges(5,2) + t206 * t309 - t269 * t310 + t233) * g(1) + (-t93 * mrSges(5,1) - t94 * mrSges(5,2) + t234 + (-t205 * t310 + t309) * t203) * g(2) + t328 * t143 + t329 * t144 + t152 * t252 - t208 + (-m(5) * t250 - t237 - m(6) * (t215 + t267) + t228 - m(4) * t267 - t321) * g(3) + t322 * t239 + t323 * t176 + t324 * t177 - (-Ifges(3,2) * t262 + t175 + t325) * t261 / 0.2e1 + t120 * t262 / 0.2e1 + t315 * t36 + t316 * t35 + (t1 * t70 + t2 * t69 + t315 * t5 + t316 * t6 - t50 * t63) * m(6) - t107 * (-Ifges(5,2) * t108 - t284) / 0.2e1 - t40 * t288 + (m(4) * t152 + m(5) * t223 - t311) * qJD(3) + t41 * t287 + t52 * t302 + t312 * t230; t133 * t10 + t130 * t9 + t197 * t59 + t198 * t60 + t326 * t36 - t327 * t35 + t311 * qJD(2) + t260 * t205 * g(3) + ((-t140 - t23 - t62) * qJD(1) - t312 * t260) * t202 + t100 + (t1 * t133 + t130 * t2 - t262 * t50 + t326 * t5 - t327 * t6) * m(6) + (-qJD(2) * t223 + t18 * t198 + t19 * t197 - t262 * t83) * m(5) + (-qJD(2) * t152 + t123 * t262 + t97) * m(4); -t107 * t84 - t108 * t85 - t314 * t35 + t221 * t36 + t247 + t248 + (t221 * t5 - t314 * t6 + t22 + t232) * m(6) + (-t107 * t41 - t108 * t40 + t232 + t30) * m(5); -t50 * (mrSges(6,1) * t221 + mrSges(6,2) * t314) + (Ifges(6,1) * t314 - t297) * t306 + t20 * t305 + (Ifges(6,5) * t314 - Ifges(6,6) * t221) * t299 - t5 * t35 + t6 * t36 - g(1) * t233 - g(2) * t234 - g(3) * t228 + (t221 * t6 + t314 * t5) * mrSges(6,3) + t208 + (-Ifges(6,2) * t221 + t21 + t46) * t307;];
tau = t11;
