% Calculate vector of inverse dynamics joint torques for
% S5RRPPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:57
% EndTime: 2019-12-31 19:35:24
% DurationCPUTime: 11.95s
% Computational Cost: add. (3664->523), mult. (8252->659), div. (0->0), fcn. (5494->10), ass. (0->249)
t266 = m(5) + m(6);
t303 = m(4) + t266;
t302 = mrSges(4,2) - mrSges(5,3);
t301 = -mrSges(5,2) + mrSges(4,1);
t152 = sin(qJ(1));
t251 = g(2) * t152;
t295 = Ifges(5,5) - Ifges(4,6);
t151 = sin(qJ(2));
t154 = cos(qJ(2));
t124 = -mrSges(3,1) * t154 + mrSges(3,2) * t151;
t146 = qJ(2) + pkin(8);
t144 = sin(t146);
t145 = cos(t146);
t300 = t302 * t144 - t301 * t145 + t124;
t155 = cos(qJ(1));
t284 = g(1) * t155 + t251;
t150 = sin(qJ(5));
t153 = cos(qJ(5));
t206 = qJD(1) * qJD(2);
t197 = t151 * t206;
t205 = qJDD(1) * t154;
t118 = -t197 + t205;
t119 = qJDD(1) * t151 + t154 * t206;
t148 = sin(pkin(8));
t228 = cos(pkin(8));
t77 = -t118 * t228 + t119 * t148;
t192 = t228 * t154;
t214 = qJD(1) * t151;
t100 = -qJD(1) * t192 + t148 * t214;
t82 = -qJD(2) * t150 + t100 * t153;
t30 = qJD(5) * t82 + qJDD(2) * t153 + t150 * t77;
t273 = t30 / 0.2e1;
t83 = qJD(2) * t153 + t100 * t150;
t31 = -qJD(5) * t83 - qJDD(2) * t150 + t153 * t77;
t272 = t31 / 0.2e1;
t78 = t148 * t118 + t119 * t228;
t74 = qJDD(5) + t78;
t271 = t74 / 0.2e1;
t299 = t118 / 0.2e1;
t298 = Ifges(6,6) * t82;
t114 = t148 * t154 + t151 * t228;
t102 = t114 * qJD(1);
t93 = qJD(5) + t102;
t297 = Ifges(6,3) * t93;
t246 = qJD(2) / 0.2e1;
t296 = Ifges(4,5) - Ifges(5,4);
t92 = Ifges(4,4) * t100;
t294 = Ifges(4,1) * t102 + Ifges(4,5) * qJD(2) + Ifges(6,5) * t83 + t297 + t298 - t92;
t240 = mrSges(4,3) * t100;
t84 = -qJD(2) * mrSges(4,2) - t240;
t242 = mrSges(5,1) * t100;
t86 = -qJD(2) * mrSges(5,3) + t242;
t293 = -t84 + t86;
t239 = mrSges(4,3) * t102;
t241 = mrSges(5,1) * t102;
t292 = qJD(2) * t301 - t239 - t241;
t39 = -mrSges(6,1) * t82 + mrSges(6,2) * t83;
t291 = -t86 + t39;
t229 = qJDD(2) / 0.2e1;
t290 = t144 * t284;
t134 = t144 * qJ(4);
t219 = t145 * t155;
t289 = pkin(3) * t219 + t155 * t134;
t140 = pkin(6) * t205;
t111 = -pkin(6) * t197 + t140;
t112 = t119 * pkin(6);
t286 = t111 * t154 + t112 * t151;
t15 = mrSges(6,1) * t74 - mrSges(6,3) * t30;
t16 = -mrSges(6,2) * t74 + mrSges(6,3) * t31;
t285 = t153 * t15 + t150 * t16;
t283 = 0.2e1 * t229;
t180 = mrSges(6,1) * t153 - mrSges(6,2) * t150;
t258 = Ifges(6,4) * t83;
t26 = Ifges(6,2) * t82 + Ifges(6,6) * t93 + t258;
t254 = pkin(4) * t100;
t149 = -qJ(3) - pkin(6);
t123 = t149 * t151;
t116 = qJD(1) * t123;
t110 = qJD(2) * pkin(2) + t116;
t125 = t149 * t154;
t117 = qJD(1) * t125;
t193 = t228 * t117;
t69 = t148 * t110 - t193;
t57 = -qJD(2) * qJ(4) - t69;
t38 = -t57 - t254;
t282 = -t153 * t26 / 0.2e1 + t180 * t38;
t281 = -m(3) * pkin(1) - mrSges(2,1) + t300;
t211 = qJD(3) * t151;
t64 = qJDD(2) * pkin(2) - qJ(3) * t119 - qJD(1) * t211 - t112;
t212 = qJD(2) * t151;
t202 = pkin(6) * t212;
t210 = qJD(3) * t154;
t71 = qJ(3) * t118 + t140 + (-t202 + t210) * qJD(1);
t22 = -t148 * t71 + t228 * t64;
t169 = qJDD(4) - t22;
t265 = pkin(3) + pkin(7);
t11 = t78 * pkin(4) - qJDD(2) * t265 + t169;
t227 = qJDD(1) * pkin(1);
t90 = -pkin(2) * t118 + qJDD(3) - t227;
t158 = -qJ(4) * t78 - qJD(4) * t102 + t90;
t8 = t265 * t77 + t158;
t255 = pkin(2) * t154;
t139 = pkin(1) + t255;
t120 = -qJD(1) * t139 + qJD(3);
t160 = -qJ(4) * t102 + t120;
t32 = t100 * t265 + t160;
t104 = t148 * t117;
t68 = t110 * t228 + t104;
t167 = qJD(4) - t68;
t248 = t102 * pkin(4);
t33 = -qJD(2) * t265 + t167 + t248;
t9 = -t150 * t32 + t153 * t33;
t1 = qJD(5) * t9 + t11 * t150 + t153 * t8;
t10 = t150 * t33 + t153 * t32;
t2 = -qJD(5) * t10 + t11 * t153 - t150 * t8;
t280 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t279 = (-g(1) * t219 - t145 * t251) * qJ(4);
t45 = pkin(3) * t100 + t160;
t278 = -t120 * mrSges(4,1) + t45 * mrSges(5,2);
t277 = -m(3) * pkin(6) - m(6) * (pkin(4) - t149) - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t276 = t9 * mrSges(6,1) + t120 * mrSges(4,2) - t10 * mrSges(6,2) - t45 * mrSges(5,3);
t274 = Ifges(6,1) * t273 + Ifges(6,4) * t272 + Ifges(6,5) * t271;
t270 = -t82 / 0.2e1;
t269 = -t83 / 0.2e1;
t268 = t83 / 0.2e1;
t267 = -t93 / 0.2e1;
t264 = -t100 / 0.2e1;
t263 = t100 / 0.2e1;
t262 = -t102 / 0.2e1;
t261 = t102 / 0.2e1;
t257 = pkin(2) * t148;
t253 = pkin(6) * t154;
t250 = g(3) * t145;
t249 = t1 * t150;
t136 = t145 * pkin(3);
t247 = -qJD(2) / 0.2e1;
t23 = t148 * t64 + t228 * t71;
t238 = mrSges(6,3) * t153;
t237 = Ifges(3,4) * t151;
t236 = Ifges(3,4) * t154;
t235 = Ifges(6,4) * t150;
t234 = Ifges(6,4) * t153;
t233 = t102 * Ifges(4,4);
t232 = t102 * Ifges(5,6);
t101 = t114 * qJD(2);
t226 = t101 * t150;
t225 = t101 * t153;
t224 = t102 * t150;
t113 = t148 * t151 - t192;
t221 = t113 * t150;
t220 = t113 * t153;
t218 = t150 * t155;
t217 = t152 * t150;
t216 = t152 * t153;
t215 = t153 * t155;
t213 = qJD(1) * t154;
t209 = qJD(5) * t150;
t208 = qJD(5) * t153;
t204 = Ifges(6,5) * t30 + Ifges(6,6) * t31 + Ifges(6,3) * t74;
t203 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t214) * t253;
t143 = pkin(2) * t212;
t142 = pkin(2) * t214;
t199 = t228 * pkin(2);
t195 = -t209 / 0.2e1;
t61 = t78 * mrSges(5,1) + qJDD(2) * mrSges(5,2);
t194 = qJD(2) * t149;
t98 = t151 * t194 + t210;
t99 = t154 * t194 - t211;
t47 = t148 * t98 - t228 * t99;
t190 = qJ(4) * t100 + t142;
t189 = -t139 - t134;
t75 = t116 * t148 - t193;
t79 = -t228 * t123 - t125 * t148;
t130 = t155 * t139;
t188 = -t152 * t149 + t130;
t187 = -m(6) * t265 - mrSges(6,3);
t186 = t134 + t136 + t255;
t138 = -t199 - pkin(3);
t184 = t10 * t153 - t9 * t150;
t182 = mrSges(3,1) * t151 + mrSges(3,2) * t154;
t179 = t150 * mrSges(6,1) + t153 * mrSges(6,2);
t177 = Ifges(6,1) * t150 + t234;
t176 = t154 * Ifges(3,2) + t237;
t175 = Ifges(6,2) * t153 + t235;
t174 = Ifges(3,5) * t154 - Ifges(3,6) * t151;
t173 = Ifges(6,5) * t150 + Ifges(6,6) * t153;
t170 = -qJ(4) * t114 - t139;
t40 = t113 * t265 + t170;
t49 = pkin(4) * t114 + t79;
t19 = t150 * t49 + t153 * t40;
t18 = -t150 * t40 + t153 * t49;
t43 = -mrSges(6,2) * t93 + mrSges(6,3) * t82;
t44 = mrSges(6,1) * t93 - mrSges(6,3) * t83;
t172 = -t150 * t44 + t153 * t43;
t171 = -t150 * t43 - t153 * t44;
t168 = pkin(1) * t182;
t48 = t148 * t99 + t228 * t98;
t166 = t113 * t208 + t226;
t165 = t113 * t209 - t225;
t103 = qJD(2) * t192 - t148 * t212;
t164 = -qJ(4) * t103 - qJD(4) * t114 + t143;
t163 = t151 * (Ifges(3,1) * t154 - t237);
t76 = t116 * t228 + t104;
t80 = t148 * t123 - t125 * t228;
t20 = -qJDD(2) * qJ(4) - qJD(2) * qJD(4) - t23;
t157 = qJD(5) * t184 + t153 * t2 + t249;
t141 = Ifges(3,4) * t213;
t135 = qJ(4) + t257;
t122 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t213;
t108 = Ifges(3,1) * t214 + Ifges(3,5) * qJD(2) + t141;
t107 = Ifges(3,6) * qJD(2) + qJD(1) * t176;
t97 = -t144 * t217 + t215;
t96 = t144 * t216 + t218;
t95 = t144 * t218 + t216;
t94 = t144 * t215 - t217;
t91 = Ifges(5,6) * t100;
t81 = Ifges(6,4) * t82;
t73 = t78 * mrSges(5,3);
t72 = t78 * mrSges(4,2);
t67 = pkin(3) * t113 + t170;
t66 = -mrSges(5,2) * t100 - mrSges(5,3) * t102;
t65 = mrSges(4,1) * t100 + mrSges(4,2) * t102;
t60 = mrSges(5,1) * t77 - qJDD(2) * mrSges(5,3);
t59 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t78;
t58 = -qJDD(2) * mrSges(4,2) - mrSges(4,3) * t77;
t55 = -Ifges(4,2) * t100 + Ifges(4,6) * qJD(2) + t233;
t54 = Ifges(5,4) * qJD(2) - Ifges(5,2) * t102 + t91;
t53 = Ifges(5,5) * qJD(2) + Ifges(5,3) * t100 - t232;
t52 = -qJD(2) * pkin(3) + t167;
t50 = -t113 * pkin(4) + t80;
t46 = pkin(3) * t102 + t190;
t42 = t76 - t248;
t41 = t75 - t254;
t37 = pkin(3) * t101 + t164;
t36 = -t101 * pkin(4) + t48;
t35 = pkin(4) * t103 + t47;
t34 = t102 * t265 + t190;
t27 = t83 * Ifges(6,1) + t93 * Ifges(6,5) + t81;
t24 = t101 * t265 + t164;
t21 = -qJDD(2) * pkin(3) + t169;
t17 = pkin(3) * t77 + t158;
t14 = t150 * t41 + t153 * t34;
t13 = -t150 * t34 + t153 * t41;
t12 = -pkin(4) * t77 - t20;
t7 = -mrSges(6,1) * t31 + mrSges(6,2) * t30;
t5 = t30 * Ifges(6,4) + t31 * Ifges(6,2) + t74 * Ifges(6,6);
t4 = -qJD(5) * t19 - t150 * t24 + t153 * t35;
t3 = qJD(5) * t18 + t150 * t35 + t153 * t24;
t6 = [t78 * Ifges(4,1) * t114 + Ifges(3,6) * t154 * t229 + t176 * t299 + t119 * t236 / 0.2e1 + t93 * (Ifges(6,5) * t166 - Ifges(6,6) * t165) / 0.2e1 + t37 * t66 + t50 * t7 + t36 * t39 + t3 * t43 + t4 * t44 + t18 * t15 + t19 * t16 + (t21 * mrSges(5,1) + t90 * mrSges(4,2) - t22 * mrSges(4,3) - t17 * mrSges(5,3) + Ifges(4,5) * t229 + Ifges(6,5) * t273 + Ifges(6,6) * t272 + Ifges(6,3) * t271 + Ifges(5,2) * t78 + (-Ifges(4,4) / 0.2e1 - Ifges(5,6)) * t77 - t283 * Ifges(5,4) + t280) * t114 + (-Ifges(4,4) * t77 + Ifges(4,5) * qJDD(2) + t204) * t114 / 0.2e1 + t38 * (mrSges(6,1) * t165 + mrSges(6,2) * t166) + (-m(4) * t68 + m(5) * t52 - t292) * t47 + (m(4) * t69 - m(5) * t57 - t293) * t48 + (-t97 * mrSges(6,1) + t96 * mrSges(6,2) + ((m(4) + m(5)) * t149 + t277) * t155 + (-m(5) * (t189 - t136) - m(6) * t189 - t145 * t187 + m(4) * t139 - t281) * t152) * g(1) + m(4) * (t120 * t143 - t139 * t90) + m(5) * (t17 * t67 + t37 * t45) + t65 * t143 + t67 * (-t77 * mrSges(5,2) - t73) + (-m(5) * (t188 + t289) - m(4) * t188 - m(6) * (pkin(7) * t219 + t130 + t289) - t95 * mrSges(6,1) - t94 * mrSges(6,2) - mrSges(6,3) * t219 + t281 * t155 + t277 * t152) * g(2) - t107 * t212 / 0.2e1 + (t118 * t253 + t286) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t286) + t26 * t225 / 0.2e1 + t27 * t226 / 0.2e1 - t124 * t227 + t82 * (Ifges(6,4) * t166 - Ifges(6,2) * t165) / 0.2e1 - t122 * t202 + (Ifges(6,1) * t166 - Ifges(6,4) * t165) * t268 - t168 * t206 + (m(4) * t23 - m(5) * t20 + t58 - t60) * t80 + (t1 * t220 - t10 * t165 - t166 * t9 - t2 * t221) * mrSges(6,3) + (t154 * (-Ifges(3,2) * t151 + t236) + t163) * t206 / 0.2e1 + (qJD(5) * t27 + t5) * t220 / 0.2e1 + (-m(4) * t22 + m(5) * t21 - t59 + t61) * t79 + (t294 / 0.2e1 - t68 * mrSges(4,3) + t52 * mrSges(5,1) - t54 / 0.2e1 + t298 / 0.2e1 + t297 / 0.2e1 + Ifges(4,1) * t261 - Ifges(5,2) * t262 - Ifges(5,6) * t263 + Ifges(4,4) * t264 + Ifges(6,5) * t268 + t296 * t246 + t276) * t103 + Ifges(2,3) * qJDD(1) + t154 * t108 * t246 + (t57 * mrSges(5,1) - t69 * mrSges(4,3) - t55 / 0.2e1 + t53 / 0.2e1 - Ifges(4,4) * t261 + Ifges(5,6) * t262 + Ifges(5,3) * t263 - Ifges(4,2) * t264 + t295 * t246 - t278) * t101 + (t174 * t246 - t203) * qJD(2) + (-pkin(6) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t119) + Ifges(3,1) * t119 + Ifges(3,4) * t299 + t283 * Ifges(3,5)) * t151 + (t90 * mrSges(4,1) + t20 * mrSges(5,1) - t17 * mrSges(5,2) - t23 * mrSges(4,3) - t12 * t180 + t173 * t271 + t175 * t272 + t177 * t273 + t26 * t195 + (-Ifges(4,4) - Ifges(5,6)) * t78 + (Ifges(4,2) + Ifges(5,3)) * t77 + t295 * t283) * t113 + m(6) * (t1 * t19 + t10 * t3 + t12 * t50 + t18 * t2 + t36 * t38 + t4 * t9) + t221 * t274 - qJDD(2) * mrSges(3,2) * t253 - pkin(1) * (-mrSges(3,1) * t118 + mrSges(3,2) * t119) - t139 * (t77 * mrSges(4,1) + t72) + t154 * (Ifges(3,4) * t119 + Ifges(3,2) * t118 + Ifges(3,6) * qJDD(2)) / 0.2e1; Ifges(3,6) * t118 - t111 * mrSges(3,2) - t112 * mrSges(3,1) + (-t10 * t208 - t249 + (t209 + t224) * t9) * mrSges(6,3) - t46 * t66 - t42 * t39 - t14 * t43 - t13 * t44 - t23 * mrSges(4,2) - t20 * mrSges(5,3) + t21 * mrSges(5,2) + t22 * mrSges(4,1) + (t7 - t60) * t135 + t291 * qJD(4) + t292 * t75 + t293 * t76 + (-Ifges(4,2) * t102 + t294 - t92) * t263 - t174 * t206 / 0.2e1 + (t203 + (-t163 / 0.2e1 + t168) * qJD(1)) * qJD(1) + (m(6) * t157 + t208 * t43 - t209 * t44 + t285) * (-pkin(7) + t138) + (t182 + t303 * pkin(2) * t151 + (-t179 + t302) * t145 + (m(5) * pkin(3) - t187 + t301) * t144) * t284 + (-m(4) * t255 - m(6) * (pkin(7) * t145 + t186) - t145 * mrSges(6,3) - t179 * t144 - m(5) * t186 + t300) * g(3) - t65 * t142 - t68 * t240 - t57 * t241 + (t107 / 0.2e1 + pkin(6) * t122) * t214 + ((t148 * t23 + t22 * t228) * pkin(2) - t120 * t142 + t68 * t75 - t69 * t76) * m(4) + (-t233 + t53) * t262 + (t232 + t55) * t261 - t2 * t238 + t12 * t179 - (-Ifges(3,2) * t214 + t108 + t141) * t213 / 0.2e1 - (t173 * t93 + t175 * t82 + t177 * t83) * qJD(5) / 0.2e1 + (-t224 / 0.2e1 + t195) * t27 + t295 * t77 + (Ifges(5,3) * t264 - t10 * t238 + t173 * t267 + t175 * t270 + t177 * t269 + t247 * t295 + t278 + t282) * t102 + t296 * t78 + (-Ifges(4,1) * t262 - Ifges(6,5) * t269 + Ifges(5,2) * t261 - Ifges(6,6) * t270 - Ifges(6,3) * t267 - t247 * t296 + t276) * t100 + t59 * t199 + t69 * t239 + t52 * t242 + t58 * t257 + (Ifges(6,5) * t153 - Ifges(6,6) * t150) * t271 + (-Ifges(6,2) * t150 + t234) * t272 + (Ifges(6,1) * t153 - t235) * t273 + t153 * t274 + (Ifges(4,3) + Ifges(5,1) + Ifges(3,3)) * qJDD(2) + (t91 + t54) * t264 + Ifges(3,5) * t119 + t138 * t61 - t150 * t5 / 0.2e1 + (-t135 * t20 + t138 * t21 - t45 * t46 - t52 * t75 + (t76 - qJD(4)) * t57 + t279) * m(5) + (-t10 * t14 + t12 * t135 - t13 * t9 + (qJD(4) - t42) * t38 + t279) * m(6) + t282 * qJD(5); -t150 * t15 + t153 * t16 + t72 - t73 + t301 * t77 + t171 * qJD(5) + (t84 + t291) * t100 + (t171 + t292) * t102 + (-g(1) * t152 + g(2) * t155) * t303 + (t1 * t153 + t100 * t38 - t2 * t150 - t93 * (t10 * t150 + t9 * t153)) * m(6) + (-t100 * t57 - t102 * t52 + t17) * m(5) + (t100 * t69 + t102 * t68 + t90) * m(4); t172 * qJD(5) - t291 * qJD(2) + t266 * t250 + (t172 + t66) * t102 + t61 + (-qJD(2) * t38 + t102 * t184 + t157 - t290) * m(6) + (qJD(2) * t57 + t102 * t45 + t21 - t290) * m(5) + t285; -t38 * (mrSges(6,1) * t83 + mrSges(6,2) * t82) + (Ifges(6,1) * t82 - t258) * t269 + t26 * t268 + (Ifges(6,5) * t82 - Ifges(6,6) * t83) * t267 - t9 * t43 + t10 * t44 - g(1) * (mrSges(6,1) * t94 - mrSges(6,2) * t95) - g(2) * (mrSges(6,1) * t96 + mrSges(6,2) * t97) + t180 * t250 + (t10 * t83 + t82 * t9) * mrSges(6,3) + t204 + (-Ifges(6,2) * t83 + t27 + t81) * t270 + t280;];
tau = t6;
