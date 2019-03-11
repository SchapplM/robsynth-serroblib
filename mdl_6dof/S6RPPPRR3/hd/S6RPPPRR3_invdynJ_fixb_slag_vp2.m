% Calculate vector of inverse dynamics joint torques for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPPRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:28
% EndTime: 2019-03-09 01:33:41
% DurationCPUTime: 7.48s
% Computational Cost: add. (5099->468), mult. (9708->628), div. (0->0), fcn. (6297->12), ass. (0->203)
t160 = sin(pkin(9));
t162 = cos(pkin(9));
t159 = sin(pkin(10));
t161 = cos(pkin(10));
t166 = sin(qJ(5));
t168 = cos(qJ(5));
t113 = t159 * t168 + t161 * t166;
t274 = t113 * qJD(5);
t111 = t159 * t166 - t168 * t161;
t275 = t111 * qJD(1);
t284 = -t160 * t274 + t162 * t275;
t197 = mrSges(5,1) * t161 - mrSges(5,2) * t159;
t297 = -m(5) * pkin(3) - mrSges(4,1) - t197;
t295 = t159 ^ 2 + t161 ^ 2;
t200 = t295 * mrSges(5,3);
t296 = mrSges(4,2) - t200;
t102 = qJD(6) - t275;
t68 = qJD(1) * t274 + qJDD(1) * t111;
t294 = -m(5) - m(4);
t293 = -m(6) - m(7);
t165 = sin(qJ(6));
t167 = cos(qJ(6));
t145 = t161 * qJD(3);
t169 = -pkin(1) - pkin(2);
t126 = qJD(1) * t169 + qJD(2);
t149 = t162 * qJ(2);
t100 = qJD(1) * t149 + t160 * t126;
t86 = -qJD(1) * qJ(4) + t100;
t59 = t145 + (pkin(7) * qJD(1) - t86) * t159;
t223 = qJD(1) * t161;
t72 = t159 * qJD(3) + t161 * t86;
t60 = -pkin(7) * t223 + t72;
t24 = t166 * t59 + t168 * t60;
t22 = qJD(5) * pkin(8) + t24;
t107 = t113 * qJD(1);
t224 = qJD(1) * t160;
t99 = -qJ(2) * t224 + t126 * t162;
t85 = qJD(1) * pkin(3) + qJD(4) - t99;
t76 = pkin(4) * t223 + t85;
t35 = -pkin(5) * t275 + pkin(8) * t107 + t76;
t11 = -t165 * t22 + t167 * t35;
t292 = t11 * mrSges(7,1);
t12 = t165 * t35 + t167 * t22;
t291 = t12 * mrSges(7,2);
t215 = m(5) - t293;
t290 = m(4) + t215;
t289 = mrSges(3,1) + mrSges(2,1);
t288 = -mrSges(3,3) + mrSges(2,2);
t67 = qJD(5) * t275 - qJDD(1) * t113;
t81 = qJD(5) * t167 + t107 * t165;
t33 = qJD(6) * t81 + qJDD(5) * t165 + t167 * t67;
t82 = qJD(5) * t165 - t107 * t167;
t34 = -qJD(6) * t82 + qJDD(5) * t167 - t165 * t67;
t13 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t247 = -qJDD(5) * mrSges(6,1) + mrSges(6,3) * t67 + t13;
t98 = t111 * t160;
t74 = -t162 * t167 + t165 * t98;
t287 = qJD(6) * t74 - t165 * t224 + t167 * t284;
t183 = t162 * t165 + t167 * t98;
t286 = qJD(6) * t183 - t165 * t284 - t167 * t224;
t240 = t107 * mrSges(6,3);
t246 = -qJD(5) * mrSges(6,1) - mrSges(7,1) * t81 + mrSges(7,2) * t82 - t240;
t220 = qJD(5) * t168;
t221 = qJD(5) * t166;
t282 = -t159 * t221 + t161 * t220;
t285 = t162 * t107 - t282 * t160;
t211 = qJDD(1) * t161;
t212 = qJDD(1) * t159;
t110 = mrSges(5,1) * t211 - mrSges(5,2) * t212;
t27 = -t68 * mrSges(6,1) + t67 * mrSges(6,2);
t283 = -t110 - t27;
t214 = qJD(1) * qJD(2);
t127 = qJDD(1) * qJ(2) + t214;
t143 = t161 * qJDD(3);
t125 = qJDD(1) * t169 + qJDD(2);
t84 = t160 * t125 + t162 * t127;
t73 = -qJDD(1) * qJ(4) - qJD(1) * qJD(4) + t84;
t50 = -t159 * t73 + t143;
t51 = t159 * qJDD(3) + t161 * t73;
t280 = t50 * t159 - t51 * t161;
t63 = qJDD(6) - t68;
t16 = mrSges(7,1) * t63 - mrSges(7,3) * t33;
t17 = -mrSges(7,2) * t63 + mrSges(7,3) * t34;
t279 = -t165 * t16 + t167 * t17;
t278 = t161 * t72 - t159 * (-t159 * t86 + t145);
t83 = t125 * t162 - t160 * t127;
t77 = qJDD(1) * pkin(3) + qJDD(4) - t83;
t65 = pkin(4) * t211 + t77;
t18 = -pkin(5) * t68 - pkin(8) * t67 + t65;
t43 = t143 + (pkin(7) * qJDD(1) - t73) * t159;
t44 = -pkin(7) * t211 + t51;
t7 = t166 * t43 + t168 * t44 + t59 * t220 - t221 * t60;
t5 = qJDD(5) * pkin(8) + t7;
t1 = qJD(6) * t11 + t165 * t18 + t167 * t5;
t2 = -qJD(6) * t12 - t165 * t5 + t167 * t18;
t277 = t1 * t167 - t165 * t2;
t8 = -qJD(5) * t24 - t166 * t44 + t168 * t43;
t23 = -t166 * t60 + t168 * t59;
t21 = -qJD(5) * pkin(5) - t23;
t271 = -m(7) * t21 - t246;
t270 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t269 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t157 = pkin(10) + qJ(5);
t146 = sin(t157);
t147 = cos(t157);
t194 = -mrSges(7,1) * t167 + mrSges(7,2) * t165;
t175 = m(7) * pkin(5) - t194;
t196 = mrSges(6,1) * t147 - mrSges(6,2) * t146;
t209 = m(7) * pkin(8) + mrSges(7,3);
t268 = t146 * t209 + t147 * t175 + t196;
t267 = m(5) * t85 + m(6) * t76 - mrSges(6,1) * t275 - mrSges(6,2) * t107 + t197 * qJD(1);
t266 = t33 / 0.2e1;
t265 = t34 / 0.2e1;
t264 = t63 / 0.2e1;
t263 = -t81 / 0.2e1;
t262 = -t82 / 0.2e1;
t261 = t82 / 0.2e1;
t260 = -t102 / 0.2e1;
t258 = -t107 / 0.2e1;
t256 = cos(qJ(1));
t255 = sin(qJ(1));
t250 = t82 * Ifges(7,4);
t120 = t160 * t169 + t149;
t115 = -qJ(4) + t120;
t248 = pkin(7) - t115;
t244 = mrSges(6,3) * t275;
t243 = Ifges(6,4) * t107;
t242 = Ifges(7,4) * t165;
t241 = Ifges(7,4) * t167;
t233 = t275 * t165;
t232 = t275 * t167;
t231 = t113 * t165;
t230 = t113 * t167;
t229 = t147 * t165;
t228 = t147 * t167;
t227 = t167 * t282;
t226 = t256 * pkin(1) + t255 * qJ(2);
t222 = qJD(2) * t160;
t219 = qJD(6) * t165;
t218 = qJD(6) * t167;
t217 = qJDD(1) * mrSges(3,1);
t216 = qJDD(1) * mrSges(4,1);
t210 = Ifges(7,5) * t33 + Ifges(7,6) * t34 + Ifges(7,3) * t63;
t29 = t81 * Ifges(7,2) + t102 * Ifges(7,6) + t250;
t207 = -t165 * t29 / 0.2e1;
t206 = t256 * pkin(2) + t226;
t201 = t218 / 0.2e1;
t119 = -t160 * qJ(2) + t162 * t169;
t117 = pkin(3) - t119;
t199 = -pkin(1) * t255 + t256 * qJ(2);
t193 = mrSges(7,1) * t165 + mrSges(7,2) * t167;
t192 = Ifges(7,1) * t167 - t242;
t191 = -Ifges(7,2) * t165 + t241;
t190 = Ifges(7,5) * t167 - Ifges(7,6) * t165;
t189 = t11 * t165 - t12 * t167;
t91 = t248 * t159;
t92 = t248 * t161;
t40 = t166 * t91 - t168 * t92;
t151 = t161 * pkin(4);
t103 = t117 + t151;
t45 = -pkin(5) * t111 + pkin(8) * t113 + t103;
t20 = t165 * t45 + t167 * t40;
t19 = -t165 * t40 + t167 * t45;
t46 = -mrSges(7,2) * t102 + mrSges(7,3) * t81;
t47 = mrSges(7,1) * t102 - mrSges(7,3) * t82;
t186 = -t165 * t47 + t167 * t46;
t39 = -t166 * t92 - t168 * t91;
t184 = t100 * t162 - t160 * t99;
t93 = -qJD(5) * mrSges(6,2) + t244;
t182 = -t186 - t93;
t181 = t21 * t193;
t180 = t113 * t218 + t165 * t282;
t179 = t113 * t219 - t227;
t174 = -pkin(2) * t255 + t199;
t172 = (-t11 * t167 - t12 * t165) * qJD(6) + t277;
t170 = qJD(1) ^ 2;
t164 = -pkin(7) - qJ(4);
t141 = -qJDD(1) * pkin(1) + qJDD(2);
t137 = t151 + pkin(3);
t129 = qJD(2) * t162 - qJD(4);
t114 = t160 * t256 - t162 * t255;
t112 = -t160 * t255 - t162 * t256;
t101 = Ifges(6,4) * t275;
t97 = t113 * t160;
t78 = Ifges(7,4) * t81;
t64 = -pkin(5) * t107 - pkin(8) * t275;
t57 = -t112 * t228 + t114 * t165;
t56 = t112 * t229 + t114 * t167;
t55 = -t107 * Ifges(6,1) + Ifges(6,5) * qJD(5) + t101;
t54 = Ifges(6,2) * t275 + Ifges(6,6) * qJD(5) - t243;
t53 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t68;
t49 = -pkin(5) * t274 + pkin(8) * t282 + t222;
t30 = Ifges(7,1) * t82 + Ifges(7,5) * t102 + t78;
t28 = t82 * Ifges(7,5) + t81 * Ifges(7,6) + t102 * Ifges(7,3);
t25 = -qJD(5) * t39 - t111 * t129;
t15 = t165 * t64 + t167 * t23;
t14 = -t165 * t23 + t167 * t64;
t10 = t33 * Ifges(7,1) + t34 * Ifges(7,4) + t63 * Ifges(7,5);
t9 = t33 * Ifges(7,4) + t34 * Ifges(7,2) + t63 * Ifges(7,6);
t6 = -qJDD(5) * pkin(5) - t8;
t4 = -qJD(6) * t20 - t165 * t25 + t167 * t49;
t3 = qJD(6) * t19 + t165 * t49 + t167 * t25;
t26 = [(t160 * t214 - t83) * mrSges(4,1) - (-Ifges(5,4) * t159 - Ifges(5,2) * t161) * t211 - (-Ifges(5,1) * t159 - Ifges(5,4) * t161) * t212 + m(3) * (-pkin(1) * t141 + (t127 + t214) * qJ(2)) - t119 * t216 + t77 * t197 + (-m(6) * t8 + m(7) * t6 + t247) * t39 + m(4) * (qJD(2) * t184 + t119 * t83 + t120 * t84) + t21 * (-mrSges(7,1) * t180 + mrSges(7,2) * t179) + pkin(1) * t217 + (-t210 / 0.2e1 - Ifges(7,3) * t264 - Ifges(7,6) * t265 - Ifges(7,5) * t266 + t7 * mrSges(6,3) + Ifges(6,4) * t67 + Ifges(6,2) * t68 + Ifges(6,6) * qJDD(5) - t65 * mrSges(6,1) + t270) * t111 + (-m(6) * t23 - t271) * (qJD(5) * t40 + t113 * t129) + t267 * t222 + (t162 * t214 + t84) * mrSges(4,2) + (qJD(6) * t30 + t9) * t231 / 0.2e1 + m(6) * (t103 * t65 + t24 * t25 + t40 * t7) + m(7) * (t1 * t20 + t11 * t4 + t12 * t3 + t19 * t2) - t274 * t292 + t81 * (Ifges(7,4) * t179 + Ifges(7,2) * t180 - Ifges(7,6) * t274) / 0.2e1 + t102 * (Ifges(7,5) * t179 + Ifges(7,6) * t180 - Ifges(7,3) * t274) / 0.2e1 + (Ifges(7,1) * t179 + Ifges(7,4) * t180 - Ifges(7,5) * t274) * t261 + t274 * t291 + t274 * t54 / 0.2e1 - t282 * t207 + (-Ifges(6,1) * t282 + Ifges(6,4) * t274) * t258 + t275 * (-Ifges(6,4) * t282 + Ifges(6,2) * t274) / 0.2e1 + (t23 * t282 + t24 * t274) * mrSges(6,3) + t76 * (-mrSges(6,1) * t274 - mrSges(6,2) * t282) + qJD(5) * (-Ifges(6,5) * t282 + Ifges(6,6) * t274) / 0.2e1 - t282 * t55 / 0.2e1 + m(5) * (-t115 * t280 + t117 * t77 + t129 * t278) + t3 * t46 + t4 * t47 + t19 * t16 + t20 * t17 + (t120 * mrSges(4,2) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1) - t10 * t230 / 0.2e1 - t30 * t227 / 0.2e1 + (t1 * t231 - t11 * t179 + t12 * t180 + t2 * t230) * mrSges(7,3) + (t280 + t295 * (-qJD(1) * t129 - qJDD(1) * t115)) * mrSges(5,3) - t274 * t28 / 0.2e1 + (-t65 * mrSges(6,2) + t8 * mrSges(6,3) - Ifges(6,1) * t67 - Ifges(6,4) * t68 - Ifges(6,5) * qJDD(5) - t190 * t264 - t191 * t265 - t192 * t266 - t193 * t6 + t201 * t29) * t113 + (-m(3) * t226 - t57 * mrSges(7,1) - t56 * mrSges(7,2) - t289 * t256 + t288 * t255 + t294 * t206 + t293 * (-t112 * t137 - t114 * t164 + t206) + t269 * t114 + (t196 - m(7) * (-pkin(5) * t147 - pkin(8) * t146) + t146 * mrSges(7,3) - t297) * t112) * g(2) + (-m(3) * t199 + t288 * t256 + t289 * t255 + t294 * t174 + t293 * (-t112 * t164 + t114 * t137 + t174) + (-t268 + t297) * t114 + (-t193 + t269) * t112) * g(1) + t40 * t53 + t25 * t93 + t103 * t27 + t117 * t110 + 0.2e1 * t127 * mrSges(3,3) - t141 * mrSges(3,1); -t217 + t74 * t16 - t183 * t17 - t98 * t53 + t247 * t97 + t284 * t93 + t286 * t47 + t287 * t46 + (-m(3) * qJ(2) - mrSges(3,3)) * t170 + (-t170 * t296 - t216 + t283) * t162 + (-t170 * mrSges(4,1) + qJDD(1) * t296) * t160 + m(3) * t141 + m(4) * (t160 * t84 + t162 * t83) + m(5) * (-t160 * t280 - t162 * t77) - t246 * t285 + (-t1 * t183 + t11 * t286 + t12 * t287 + t2 * t74 - t21 * t285 + t6 * t97) * m(7) + (-t162 * t65 + t23 * t285 + t24 * t284 - t7 * t98 - t8 * t97) * m(6) + (-m(5) * t162 * t278 - m(4) * t184 - t160 * t267) * qJD(1) + (-t255 * g(1) + t256 * g(2)) * (m(3) + t290); m(4) * qJDD(3) + t247 * t111 + t246 * t274 - t182 * t282 + (t53 + (-t165 * t46 - t167 * t47) * qJD(6) + t279) * t113 + m(5) * (t159 * t51 + t161 * t50) + m(6) * (-t111 * t8 + t113 * t7 - t23 * t274 + t24 * t282) + m(7) * (t111 * t6 + t113 * t172 - t189 * t282 + t21 * t274) + t290 * g(3); t186 * qJD(6) + t182 * t275 + t246 * t107 + t167 * t16 + t165 * t17 - t170 * t200 + (-g(1) * t114 + g(2) * t112) * t215 + (t1 * t165 - t102 * t189 + t107 * t21 + t2 * t167) * m(7) + (-t107 * t23 - t24 * t275 + t65) * m(6) + (qJD(1) * t278 + t77) * m(5) - t283; (t244 - t93) * t23 - t107 * t291 + t6 * t194 + (-t232 / 0.2e1 + t201) * t30 + t54 * t258 + (Ifges(7,5) * t165 + Ifges(7,6) * t167) * t264 + (Ifges(7,2) * t167 + t242) * t265 + (Ifges(7,1) * t165 + t241) * t266 + ((-t219 + t233) * t12 + (-t218 + t232) * t11 + t277) * mrSges(7,3) + (m(7) * t172 - t218 * t47 - t219 * t46 + t279) * pkin(8) + (-t240 + t271) * t24 + t268 * g(3) + (t181 + t207) * qJD(6) + (t102 * t190 + t191 * t81 + t192 * t82) * qJD(6) / 0.2e1 + t107 * t292 + ((mrSges(6,2) - t209) * t147 + (mrSges(6,1) + t175) * t146) * (-g(1) * t112 - g(2) * t114) - t275 * t181 + (-Ifges(7,3) * t107 + t190 * t275) * t260 + (-Ifges(7,5) * t107 + t192 * t275) * t262 + (-Ifges(7,6) * t107 + t191 * t275) * t263 + (Ifges(6,1) * t275 + t243 + t28) * t107 / 0.2e1 - (Ifges(6,2) * t107 + t101 + t55) * t275 / 0.2e1 - qJD(5) * (Ifges(6,5) * t275 + Ifges(6,6) * t107) / 0.2e1 - t76 * (-mrSges(6,1) * t107 + mrSges(6,2) * t275) - t15 * t46 - t14 * t47 - pkin(5) * t13 - t7 * mrSges(6,2) + t8 * mrSges(6,1) + t29 * t233 / 0.2e1 + Ifges(6,3) * qJDD(5) + Ifges(6,5) * t67 + Ifges(6,6) * t68 + (-pkin(5) * t6 - t11 * t14 - t12 * t15) * m(7) + t165 * t10 / 0.2e1 + t167 * t9 / 0.2e1; -t21 * (mrSges(7,1) * t82 + mrSges(7,2) * t81) + (Ifges(7,1) * t81 - t250) * t262 + t29 * t261 + (Ifges(7,5) * t81 - Ifges(7,6) * t82) * t260 - t11 * t46 + t12 * t47 - g(1) * (mrSges(7,1) * t56 - mrSges(7,2) * t57) - g(2) * ((-t112 * t167 + t114 * t229) * mrSges(7,1) + (t112 * t165 + t114 * t228) * mrSges(7,2)) - g(3) * t193 * t146 + (t11 * t81 + t12 * t82) * mrSges(7,3) + t210 + (-Ifges(7,2) * t82 + t30 + t78) * t263 - t270;];
tau  = t26;
