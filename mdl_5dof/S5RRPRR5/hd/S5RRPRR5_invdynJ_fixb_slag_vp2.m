% Calculate vector of inverse dynamics joint torques for
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
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
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:12
% EndTime: 2022-01-20 11:02:25
% DurationCPUTime: 4.33s
% Computational Cost: add. (5951->373), mult. (9273->501), div. (0->0), fcn. (6441->16), ass. (0->180)
t222 = cos(qJ(2));
t270 = pkin(1) * qJD(1);
t250 = t222 * t270;
t235 = qJD(3) - t250;
t213 = sin(pkin(9));
t217 = sin(qJ(4));
t214 = cos(pkin(9));
t221 = cos(qJ(4));
t255 = t214 * t221;
t156 = -t213 * t217 + t255;
t215 = -pkin(7) - qJ(3);
t167 = t215 * t213;
t201 = t214 * pkin(7);
t168 = qJ(3) * t214 + t201;
t252 = qJD(4) * t221;
t289 = t167 * t252 + qJD(3) * t255 + (-qJD(3) * t213 - qJD(4) * t168) * t217 - t156 * t250;
t113 = t217 * t167 + t221 * t168;
t157 = t213 * t221 + t214 * t217;
t288 = -t113 * qJD(4) - t235 * t157;
t209 = pkin(9) + qJ(4);
t198 = qJ(5) + t209;
t187 = sin(t198);
t188 = cos(t198);
t305 = t188 * mrSges(6,1) - t187 * mrSges(6,2);
t234 = -t214 * mrSges(4,1) + t213 * mrSges(4,2);
t196 = sin(t209);
t197 = cos(t209);
t300 = -t197 * mrSges(5,1) + t196 * mrSges(5,2) - t305;
t304 = -mrSges(3,1) + t234 + t300;
t143 = t156 * qJD(4);
t277 = pkin(8) * t143;
t303 = -t277 + t288;
t144 = t157 * qJD(4);
t141 = t144 * pkin(8);
t302 = t141 - t289;
t204 = qJDD(4) + qJDD(5);
t210 = qJD(4) + qJD(5);
t220 = cos(qJ(5));
t216 = sin(qJ(5));
t211 = qJD(1) + qJD(2);
t131 = t156 * t211;
t218 = sin(qJ(2));
t248 = t218 * t270;
t162 = qJ(3) * t211 + t248;
t245 = pkin(7) * t211 + t162;
t119 = t245 * t213;
t120 = t245 * t214;
t71 = -t119 * t217 + t120 * t221;
t53 = pkin(8) * t131 + t71;
t263 = t216 * t53;
t132 = t157 * t211;
t259 = t120 * t217;
t70 = -t221 * t119 - t259;
t52 = -pkin(8) * t132 + t70;
t50 = qJD(4) * pkin(4) + t52;
t22 = t220 * t50 - t263;
t262 = t220 * t53;
t23 = t216 * t50 + t262;
t239 = t220 * t131 - t132 * t216;
t86 = t131 * t216 + t132 * t220;
t279 = Ifges(6,4) * t86;
t269 = pkin(1) * qJD(2);
t247 = qJD(1) * t269;
t260 = pkin(1) * qJDD(1);
t152 = t218 * t260 + t222 * t247;
t205 = qJDD(1) + qJDD(2);
t118 = qJ(3) * t205 + qJD(3) * t211 + t152;
t246 = pkin(7) * t205 + t118;
t100 = t246 * t213;
t101 = t246 * t214;
t37 = -t71 * qJD(4) - t221 * t100 - t101 * t217;
t91 = t211 * t143 + t157 * t205;
t18 = qJDD(4) * pkin(4) - pkin(8) * t91 + t37;
t36 = -qJD(4) * t259 - t217 * t100 + t221 * t101 - t119 * t252;
t92 = -t211 * t144 + t156 * t205;
t19 = pkin(8) * t92 + t36;
t3 = t22 * qJD(5) + t18 * t216 + t19 * t220;
t34 = t239 * qJD(5) + t216 * t92 + t220 * t91;
t35 = -t86 * qJD(5) - t216 * t91 + t220 * t92;
t4 = -t23 * qJD(5) + t18 * t220 - t19 * t216;
t78 = Ifges(6,4) * t239;
t43 = Ifges(6,1) * t86 + Ifges(6,5) * t210 + t78;
t190 = pkin(3) * t214 + pkin(2);
t129 = -t190 * t211 + t235;
t93 = -pkin(4) * t131 + t129;
t301 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t34 + Ifges(6,6) * t35 + Ifges(6,3) * t204 - (Ifges(6,5) * t239 - Ifges(6,6) * t86) * t210 / 0.2e1 + (t22 * t239 + t23 * t86) * mrSges(6,3) - (-Ifges(6,2) * t86 + t43 + t78) * t239 / 0.2e1 - t93 * (mrSges(6,1) * t86 + mrSges(6,2) * t239) - (Ifges(6,1) * t239 - t279) * t86 / 0.2e1;
t253 = t213 ^ 2 + t214 ^ 2;
t243 = t253 * t118;
t212 = qJ(1) + qJ(2);
t199 = sin(t212);
t200 = cos(t212);
t286 = g(1) * t200 + g(2) * t199;
t298 = mrSges(3,2) - mrSges(6,3) - mrSges(5,3) - mrSges(4,3);
t42 = Ifges(6,2) * t239 + Ifges(6,6) * t210 + t279;
t296 = t42 / 0.2e1;
t112 = t221 * t167 - t168 * t217;
t276 = pkin(8) * t157;
t94 = t112 - t276;
t149 = t156 * pkin(8);
t95 = t149 + t113;
t45 = -t216 * t95 + t220 * t94;
t295 = t45 * qJD(5) + t303 * t216 - t302 * t220;
t46 = t216 * t94 + t220 * t95;
t294 = -t46 * qJD(5) + t302 * t216 + t303 * t220;
t189 = pkin(1) * t218 + qJ(3);
t146 = (-pkin(7) - t189) * t213;
t147 = t189 * t214 + t201;
t103 = t217 * t146 + t221 * t147;
t283 = t86 / 0.2e1;
t281 = t132 / 0.2e1;
t278 = pkin(1) * t222;
t273 = t144 * pkin(4);
t271 = Ifges(5,4) * t132;
t261 = mrSges(5,1) * t131 - mrSges(5,2) * t132 - t234 * t211;
t257 = t205 * t213;
t256 = t205 * t214;
t254 = t200 * pkin(2) + t199 * qJ(3);
t251 = m(4) + m(5) + m(6);
t249 = t218 * t269;
t47 = -t92 * mrSges(5,1) + t91 * mrSges(5,2);
t9 = -t35 * mrSges(6,1) + t34 * mrSges(6,2);
t244 = t253 * mrSges(4,3);
t242 = t253 * t162;
t174 = t222 * t269 + qJD(3);
t241 = t253 * t174;
t240 = t253 * t205;
t102 = t221 * t146 - t147 * t217;
t159 = pkin(4) * t197 + t190;
t206 = pkin(8) - t215;
t238 = t200 * t159 + t199 * t206;
t237 = t200 * t190 - t199 * t215;
t236 = qJD(3) * t253;
t139 = -mrSges(4,1) * t256 + mrSges(4,2) * t257;
t151 = -t218 * t247 + t222 * t260;
t233 = mrSges(6,1) * t187 + mrSges(6,2) * t188;
t76 = t102 - t276;
t77 = t149 + t103;
t38 = -t216 * t77 + t220 * t76;
t39 = t216 * t76 + t220 * t77;
t105 = t156 * t220 - t157 * t216;
t106 = t156 * t216 + t157 * t220;
t128 = -pkin(4) * t156 - t190;
t230 = qJDD(3) - t151;
t60 = t146 * t252 + t174 * t255 + (-qJD(4) * t147 - t174 * t213) * t217;
t110 = -t190 * t205 + t230;
t61 = -t103 * qJD(4) - t157 * t174;
t226 = t298 * t199 + t304 * t200;
t225 = (m(4) * pkin(2) + m(5) * t190 + m(6) * t159 - t304) * t199 + (-m(4) * qJ(3) + m(5) * t215 - m(6) * t206 + t298) * t200;
t130 = -pkin(2) * t205 + t230;
t56 = t105 * qJD(5) + t143 * t220 - t144 * t216;
t57 = -t106 * qJD(5) - t143 * t216 - t144 * t220;
t58 = -pkin(4) * t92 + t110;
t81 = t131 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t271;
t127 = Ifges(5,4) * t131;
t82 = t132 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t127;
t224 = (t58 * mrSges(6,2) - t4 * mrSges(6,3) + Ifges(6,1) * t34 + Ifges(6,4) * t35 + Ifges(6,5) * t204) * t106 + (-t58 * mrSges(6,1) + t3 * mrSges(6,3) + Ifges(6,4) * t34 + Ifges(6,2) * t35 + Ifges(6,6) * t204) * t105 + (-t70 * t143 - t71 * t144) * mrSges(5,3) + (Ifges(4,4) * t213 + Ifges(4,2) * t214) * t256 + (Ifges(4,1) * t213 + Ifges(4,4) * t214) * t257 + t130 * t234 + (-t22 * t56 + t23 * t57) * mrSges(6,3) + mrSges(4,3) * t243 + t57 * t296 + (t110 * mrSges(5,2) - t37 * mrSges(5,3) + Ifges(5,1) * t91 + Ifges(5,4) * t92 + Ifges(5,5) * qJDD(4)) * t157 + (-t110 * mrSges(5,1) + t36 * mrSges(5,3) + Ifges(5,4) * t91 + Ifges(5,2) * t92 + Ifges(5,6) * qJDD(4)) * t156 + (Ifges(5,1) * t143 - Ifges(5,4) * t144) * t281 + qJD(4) * (Ifges(5,5) * t143 - Ifges(5,6) * t144) / 0.2e1 + t129 * (mrSges(5,1) * t144 + mrSges(5,2) * t143) + t131 * (Ifges(5,4) * t143 - Ifges(5,2) * t144) / 0.2e1 + t239 * (Ifges(6,4) * t56 + Ifges(6,2) * t57) / 0.2e1 + t56 * t43 / 0.2e1 + (Ifges(6,1) * t56 + Ifges(6,4) * t57) * t283 + t93 * (-mrSges(6,1) * t57 + mrSges(6,2) * t56) + t143 * t82 / 0.2e1 - t144 * t81 / 0.2e1 + t151 * mrSges(3,1) - t152 * mrSges(3,2) + Ifges(3,3) * t205 + t210 * (Ifges(6,5) * t56 + Ifges(6,6) * t57) / 0.2e1;
t223 = cos(qJ(1));
t219 = sin(qJ(1));
t202 = t223 * pkin(1);
t192 = -pkin(2) - t278;
t166 = -t190 - t278;
t158 = -pkin(2) * t211 + t235;
t121 = t249 + t273;
t117 = t128 - t278;
t116 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t132;
t115 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t131;
t80 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t92;
t79 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t91;
t73 = mrSges(6,1) * t210 - mrSges(6,3) * t86;
t72 = -mrSges(6,2) * t210 + mrSges(6,3) * t239;
t49 = t61 - t277;
t48 = -t141 + t60;
t44 = -mrSges(6,1) * t239 + mrSges(6,2) * t86;
t29 = -mrSges(6,2) * t204 + mrSges(6,3) * t35;
t28 = mrSges(6,1) * t204 - mrSges(6,3) * t34;
t25 = t220 * t52 - t263;
t24 = -t216 * t52 - t262;
t8 = -t39 * qJD(5) - t216 * t48 + t220 * t49;
t7 = t38 * qJD(5) + t216 * t49 + t220 * t48;
t1 = [t224 + m(6) * (t117 * t58 + t121 * t93 + t22 * t8 + t23 * t7 + t3 * t39 + t38 * t4) + m(5) * (t102 * t37 + t103 * t36 + t110 * t166 + t129 * t249 + t60 * t71 + t61 * t70) + m(4) * (t130 * t192 + t158 * t249 + t162 * t241 + t189 * t243) + (-m(6) * (t202 + t238) - m(4) * (t202 + t254) - m(5) * (t202 + t237) - mrSges(2,1) * t223 + t219 * mrSges(2,2) + t226) * g(2) + (t189 * t240 + t211 * t241) * mrSges(4,3) + m(3) * (t151 * t222 + t152 * t218) * pkin(1) + (-g(2) * m(3) * t223 + (mrSges(3,1) * t222 - mrSges(3,2) * t218) * t205 + (-mrSges(3,2) * t211 * t222 + (-mrSges(3,1) * t211 - t261) * t218) * qJD(2)) * pkin(1) + t38 * t28 + t39 * t29 + t7 * t72 + t8 * t73 + (mrSges(2,2) * t223 + (mrSges(2,1) + (m(3) + t251) * pkin(1)) * t219 + t225) * g(1) + t102 * t79 + t103 * t80 + t60 * t115 + t61 * t116 + t117 * t9 + t121 * t44 + Ifges(2,3) * qJDD(1) + t166 * t47 + t192 * t139; (qJ(3) * t240 + t211 * t236) * mrSges(4,3) + t225 * g(1) + t226 * g(2) + t224 + t288 * t116 + t289 * t115 + t294 * t73 + t295 * t72 + t44 * t273 + ((-t44 + t261) * t218 + (mrSges(3,1) * t218 + (mrSges(3,2) - t244) * t222) * t211) * t270 + t45 * t28 + t46 * t29 + t112 * t79 + t113 * t80 + t128 * t9 - pkin(2) * t139 - t190 * t47 + (-t238 * g(2) + t128 * t58 + t3 * t46 + t4 * t45 + (-t248 + t273) * t93 + t295 * t23 + t294 * t22) * m(6) + (-t237 * g(2) - t110 * t190 + t112 * t37 + t113 * t36 - t129 * t248 + t288 * t70 + t289 * t71) * m(5) + (-t254 * g(2) - pkin(2) * t130 + qJ(3) * t243 + t162 * t236 - (t158 * t218 + t222 * t242) * t270) * m(4); -t211 ^ 2 * t244 - t131 * t115 + t132 * t116 - t239 * t72 + t86 * t73 + t139 + t47 + t9 + (-g(1) * t199 + g(2) * t200) * t251 + (t22 * t86 - t23 * t239 + t58) * m(6) + (-t131 * t71 + t132 * t70 + t110) * m(5) + (-t211 * t242 + t130) * m(4); -(-Ifges(5,2) * t132 + t127 + t82) * t131 / 0.2e1 + t86 * t296 + t301 + (-t132 * t44 + t216 * t29 + t220 * t28 + (-g(3) * t197 - t132 * t93 + t286 * t196 + t216 * t3 + t220 * t4) * m(6) + (-t216 * t73 + t220 * t72 + (-t216 * t22 + t220 * t23) * m(6)) * qJD(5)) * pkin(4) + t300 * g(3) + t286 * (mrSges(5,1) * t196 + mrSges(5,2) * t197 + t233) - m(6) * (t22 * t24 + t23 * t25) - t132 * (Ifges(5,1) * t131 - t271) / 0.2e1 - t36 * mrSges(5,2) + t37 * mrSges(5,1) - t25 * t72 - t24 * t73 + t81 * t281 + Ifges(5,5) * t91 + Ifges(5,6) * t92 - t70 * t115 + t71 * t116 - qJD(4) * (Ifges(5,5) * t131 - Ifges(5,6) * t132) / 0.2e1 - t129 * (mrSges(5,1) * t132 + mrSges(5,2) * t131) + (t131 * t70 + t132 * t71) * mrSges(5,3) + Ifges(5,3) * qJDD(4); -g(3) * t305 - t22 * t72 + t23 * t73 + t286 * t233 + t42 * t283 + t301;];
tau = t1;
