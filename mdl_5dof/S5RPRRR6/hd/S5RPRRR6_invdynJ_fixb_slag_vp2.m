% Calculate vector of inverse dynamics joint torques for
% S5RPRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:00:57
% EndTime: 2019-12-31 19:01:13
% DurationCPUTime: 7.87s
% Computational Cost: add. (4928->476), mult. (10327->661), div. (0->0), fcn. (6728->14), ass. (0->219)
t174 = qJ(3) + qJ(4);
t167 = sin(t174);
t168 = cos(t174);
t177 = sin(qJ(5));
t279 = mrSges(6,2) * t177;
t330 = -t167 * t279 + t168 * (-m(6) * pkin(8) - mrSges(6,3));
t316 = t168 * pkin(4) + t167 * pkin(8);
t329 = m(6) * t316;
t178 = sin(qJ(4));
t179 = sin(qJ(3));
t182 = cos(qJ(4));
t183 = cos(qJ(3));
t132 = t178 * t183 + t179 * t182;
t124 = t132 * qJD(1);
t172 = qJD(3) + qJD(4);
t181 = cos(qJ(5));
t101 = -t124 * t177 + t172 * t181;
t102 = t124 * t181 + t172 * t177;
t277 = mrSges(5,3) * t124;
t328 = mrSges(5,1) * t172 + mrSges(6,1) * t101 - mrSges(6,2) * t102 - t277;
t327 = -t168 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t167;
t240 = qJD(1) * qJD(3);
t133 = qJDD(1) * t183 - t179 * t240;
t134 = qJDD(1) * t179 + t183 * t240;
t195 = t132 * qJD(4);
t74 = -qJD(1) * t195 + t133 * t182 - t134 * t178;
t171 = qJDD(3) + qJDD(4);
t131 = t178 * t179 - t182 * t183;
t194 = t131 * qJD(4);
t73 = -qJD(1) * t194 + t133 * t178 + t134 * t182;
t33 = qJD(5) * t101 + t171 * t177 + t181 * t73;
t71 = qJDD(5) - t74;
t17 = mrSges(6,1) * t71 - mrSges(6,3) * t33;
t34 = -qJD(5) * t102 + t171 * t181 - t177 * t73;
t18 = -mrSges(6,2) * t71 + mrSges(6,3) * t34;
t205 = -t177 * t17 + t181 * t18;
t241 = qJD(5) * t181;
t242 = qJD(5) * t177;
t123 = t131 * qJD(1);
t120 = qJD(5) + t123;
t75 = -mrSges(6,2) * t120 + mrSges(6,3) * t101;
t76 = mrSges(6,1) * t120 - mrSges(6,3) * t102;
t326 = -t76 * t241 - t75 * t242 + t205;
t303 = t33 / 0.2e1;
t302 = t34 / 0.2e1;
t272 = t102 * Ifges(6,4);
t47 = t101 * Ifges(6,2) + t120 * Ifges(6,6) + t272;
t323 = -t47 / 0.2e1;
t300 = t71 / 0.2e1;
t322 = -m(4) - m(3);
t321 = -m(5) - m(6);
t320 = t133 / 0.2e1;
t13 = -mrSges(6,1) * t34 + mrSges(6,2) * t33;
t59 = mrSges(5,1) * t171 - mrSges(5,3) * t73;
t319 = t13 - t59;
t212 = mrSges(6,1) * t177 + mrSges(6,2) * t181;
t166 = t183 * qJD(2);
t175 = sin(pkin(9));
t151 = pkin(1) * t175 + pkin(6);
t143 = t151 * qJD(1);
t220 = pkin(7) * qJD(1) + t143;
t103 = -t179 * t220 + t166;
t100 = qJD(3) * pkin(3) + t103;
t247 = qJD(2) * t179;
t104 = t183 * t220 + t247;
t260 = t104 * t178;
t54 = t100 * t182 - t260;
t50 = -pkin(4) * t172 - t54;
t318 = t212 * t50;
t176 = cos(pkin(9));
t152 = -pkin(1) * t176 - pkin(2);
t169 = t183 * pkin(3);
t140 = t152 - t169;
t252 = t168 * t181;
t253 = t168 * t177;
t315 = -mrSges(6,1) * t252 + mrSges(6,2) * t253 + t327;
t96 = -qJD(3) * t131 - t194;
t199 = t132 * t241 + t177 * t96;
t141 = t151 * qJDD(1);
t246 = qJD(3) * t179;
t84 = qJD(3) * t166 + t179 * qJDD(2) + t183 * t141 - t143 * t246;
t114 = t143 * t183 + t247;
t85 = -t114 * qJD(3) + t183 * qJDD(2) - t141 * t179;
t313 = -t179 * t85 + t183 * t84;
t278 = mrSges(5,3) * t123;
t106 = -mrSges(5,2) * t172 - t278;
t310 = -t177 * t76 + t181 * t75 + t106;
t251 = t182 * t104;
t55 = t100 * t178 + t251;
t51 = pkin(8) * t172 + t55;
t125 = t140 * qJD(1);
t72 = pkin(4) * t123 - pkin(8) * t124 + t125;
t20 = -t177 * t51 + t181 * t72;
t21 = t177 * t72 + t181 * t51;
t243 = qJD(4) * t182;
t244 = qJD(4) * t178;
t62 = qJDD(3) * pkin(3) - pkin(7) * t134 + t85;
t65 = pkin(7) * t133 + t84;
t15 = t100 * t243 - t104 * t244 + t178 * t62 + t182 * t65;
t11 = pkin(8) * t171 + t15;
t142 = t152 * qJDD(1);
t105 = -pkin(3) * t133 + t142;
t19 = -pkin(4) * t74 - pkin(8) * t73 + t105;
t3 = -qJD(5) * t21 - t11 * t177 + t181 * t19;
t285 = t177 * t3;
t309 = -t20 * t241 - t21 * t242 - t285;
t16 = -qJD(4) * t55 - t178 * t65 + t182 * t62;
t147 = -mrSges(4,1) * t183 + mrSges(4,2) * t179;
t308 = m(4) * pkin(2) + mrSges(3,1) - t147 - t327;
t2 = qJD(5) * t20 + t11 * t181 + t177 * t19;
t307 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t306 = -m(4) * pkin(6) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t305 = Ifges(6,1) * t303 + Ifges(6,4) * t302 + Ifges(6,5) * t300;
t304 = m(6) * pkin(4);
t299 = -t101 / 0.2e1;
t298 = -t102 / 0.2e1;
t297 = t102 / 0.2e1;
t296 = -t120 / 0.2e1;
t293 = t124 / 0.2e1;
t292 = t181 / 0.2e1;
t180 = sin(qJ(1));
t291 = pkin(1) * t180;
t290 = pkin(3) * t178;
t289 = pkin(3) * t179;
t288 = pkin(3) * t182;
t284 = t181 * t2;
t184 = cos(qJ(1));
t170 = t184 * pkin(1);
t282 = pkin(7) + t151;
t281 = mrSges(6,1) * t181;
t280 = mrSges(5,2) * t168;
t276 = Ifges(4,4) * t179;
t275 = Ifges(4,4) * t183;
t274 = Ifges(6,4) * t177;
t273 = Ifges(6,4) * t181;
t271 = t124 * Ifges(5,4);
t259 = t123 * t177;
t258 = t123 * t181;
t257 = t132 * t177;
t256 = t132 * t181;
t249 = qJD(1) * t179;
t248 = qJD(1) * t183;
t245 = qJD(3) * t183;
t239 = Ifges(6,5) * t33 + Ifges(6,6) * t34 + Ifges(6,3) * t71;
t237 = pkin(3) * t249;
t236 = pkin(3) * t246;
t233 = mrSges(4,3) * t249;
t232 = mrSges(4,3) * t248;
t98 = Ifges(6,4) * t101;
t48 = t102 * Ifges(6,1) + t120 * Ifges(6,5) + t98;
t227 = t48 * t292;
t224 = -t242 / 0.2e1;
t223 = qJD(3) * t282;
t173 = qJ(1) + pkin(9);
t162 = sin(t173);
t219 = t330 * t162;
t163 = cos(t173);
t218 = t330 * t163;
t217 = t183 * t223;
t93 = pkin(4) * t124 + pkin(8) * t123;
t214 = mrSges(4,1) * t179 + mrSges(4,2) * t183;
t211 = Ifges(6,1) * t181 - t274;
t210 = t183 * Ifges(4,2) + t276;
t209 = -Ifges(6,2) * t177 + t273;
t208 = Ifges(4,5) * t183 - Ifges(4,6) * t179;
t207 = Ifges(6,5) * t181 - Ifges(6,6) * t177;
t204 = -t20 * t177 + t21 * t181;
t86 = pkin(4) * t131 - pkin(8) * t132 + t140;
t126 = t282 * t179;
t127 = t282 * t183;
t90 = -t126 * t178 + t127 * t182;
t37 = -t177 * t90 + t181 * t86;
t38 = t177 * t86 + t181 * t90;
t202 = -t182 * t126 - t127 * t178;
t198 = t132 * t242 - t181 * t96;
t197 = t152 * qJD(1) * t214;
t196 = t179 * (Ifges(4,1) * t183 - t276);
t191 = -t285 + (-t21 * t177 - t20 * t181) * qJD(5);
t190 = m(6) * (-pkin(4) * t167 - t289) - t167 * t281;
t189 = t280 + (mrSges(5,1) + t281 + t304) * t167;
t188 = t191 + t284;
t119 = Ifges(5,4) * t123;
t12 = -pkin(4) * t171 - t16;
t46 = t102 * Ifges(6,5) + t101 * Ifges(6,6) + t120 * Ifges(6,3);
t8 = t33 * Ifges(6,4) + t34 * Ifges(6,2) + t71 * Ifges(6,6);
t81 = -t123 * Ifges(5,2) + t172 * Ifges(5,6) + t271;
t82 = t124 * Ifges(5,1) + t172 * Ifges(5,5) - t119;
t187 = Ifges(5,5) * t73 + Ifges(5,6) * t74 - t15 * mrSges(5,2) + t16 * mrSges(5,1) + (t101 * t209 + t102 * t211 + t120 * t207) * qJD(5) / 0.2e1 - (-Ifges(5,1) * t123 - t271 + t46) * t124 / 0.2e1 + (-Ifges(5,2) * t124 - t119 + t82) * t123 / 0.2e1 - t125 * (mrSges(5,1) * t124 - mrSges(5,2) * t123) - t172 * (-Ifges(5,5) * t123 - Ifges(5,6) * t124) / 0.2e1 + (Ifges(6,3) * t124 - t123 * t207) * t296 + (Ifges(6,5) * t124 - t123 * t211) * t298 + (Ifges(6,6) * t124 - t123 * t209) * t299 - t20 * (mrSges(6,1) * t124 + mrSges(6,3) * t258) - t21 * (-mrSges(6,2) * t124 + mrSges(6,3) * t259) + (t227 + t318) * qJD(5) + (Ifges(6,5) * t177 + Ifges(6,6) * t181) * t300 + (Ifges(6,2) * t181 + t274) * t302 + (Ifges(6,1) * t177 + t273) * t303 + t177 * t305 + t47 * t224 + Ifges(5,3) * t171 + t123 * t318 + t259 * t323 - t54 * t278 + mrSges(6,3) * t284 + t8 * t292 + t81 * t293 + t48 * t258 / 0.2e1 + t12 * (t279 - t281);
t185 = -pkin(7) - pkin(6);
t160 = Ifges(4,4) * t248;
t159 = t169 + pkin(2);
t158 = -pkin(4) - t288;
t146 = -qJD(3) * mrSges(4,2) + t232;
t144 = qJD(3) * mrSges(4,1) - t233;
t122 = Ifges(4,1) * t249 + Ifges(4,5) * qJD(3) + t160;
t121 = Ifges(4,6) * qJD(3) + qJD(1) * t210;
t118 = t179 * t223;
t117 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t134;
t116 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t133;
t113 = -t143 * t179 + t166;
t111 = t162 * t177 + t163 * t252;
t110 = t162 * t181 - t163 * t253;
t109 = -t162 * t252 + t163 * t177;
t108 = t162 * t253 + t163 * t181;
t97 = qJD(3) * t132 + t195;
t92 = mrSges(5,1) * t123 + mrSges(5,2) * t124;
t80 = t93 + t237;
t64 = t103 * t182 - t260;
t63 = t103 * t178 + t251;
t60 = -mrSges(5,2) * t171 + mrSges(5,3) * t74;
t49 = pkin(4) * t97 - pkin(8) * t96 + t236;
t41 = qJD(4) * t202 - t182 * t118 - t178 * t217;
t25 = t177 * t93 + t181 * t54;
t24 = -t177 * t54 + t181 * t93;
t23 = t177 * t80 + t181 * t64;
t22 = -t177 * t64 + t181 * t80;
t5 = -qJD(5) * t38 - t177 * t41 + t181 * t49;
t4 = qJD(5) * t37 + t177 * t49 + t181 * t41;
t1 = [t41 * t106 - t97 * t81 / 0.2e1 + t96 * t82 / 0.2e1 + t97 * t46 / 0.2e1 + t90 * t60 + t4 * t75 + t5 * t76 + t37 * t17 + t38 * t18 + (m(4) * t152 + t147) * t142 + m(6) * (t2 * t38 + t20 * t5 + t21 * t4 + t3 * t37) + m(5) * (t105 * t140 + t125 * t236 + t15 * t90 + t41 * t55) + t183 * (Ifges(4,4) * t134 + Ifges(4,2) * t133) / 0.2e1 + t20 * mrSges(6,1) * t97 - t21 * mrSges(6,2) * t97 + (-m(5) * t54 + m(6) * t50 - t328) * (qJD(4) * t90 - t118 * t178 + t182 * t217) + (t198 * t20 - t199 * t21 - t2 * t257 - t256 * t3) * mrSges(6,3) + (t105 * mrSges(5,2) - t16 * mrSges(5,3) + Ifges(5,1) * t73 + Ifges(5,4) * t74 + Ifges(5,5) * t171 + t12 * t212 + t207 * t300 + t209 * t302 + t211 * t303 + t224 * t48) * t132 + (mrSges(2,1) * t180 - t109 * mrSges(6,1) + mrSges(2,2) * t184 - t108 * mrSges(6,2) - t322 * t291 + t321 * (-t163 * t185 - t291) + t306 * t163 + (m(5) * t159 - m(6) * (-t159 - t316) + t308) * t162) * g(1) + t92 * t236 + (Ifges(2,3) + Ifges(3,3) + (0.2e1 * mrSges(3,1) * t176 - 0.2e1 * mrSges(3,2) * t175 + m(3) * (t175 ^ 2 + t176 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t50 * (mrSges(6,1) * t199 - mrSges(6,2) * t198) + t120 * (-Ifges(6,5) * t198 - Ifges(6,6) * t199 + Ifges(6,3) * t97) / 0.2e1 + t101 * (-Ifges(6,4) * t198 - Ifges(6,2) * t199 + Ifges(6,6) * t97) / 0.2e1 + (t183 * (-Ifges(4,2) * t179 + t275) + t196) * t240 / 0.2e1 + (-t15 * mrSges(5,3) + Ifges(6,3) * t300 + Ifges(6,6) * t302 + Ifges(6,5) * t303 + t105 * mrSges(5,1) - Ifges(5,4) * t73 - Ifges(5,2) * t74 - Ifges(5,6) * t171 + t239 / 0.2e1 + t307) * t131 + t256 * t305 + (-t54 * t96 - t55 * t97) * mrSges(5,3) + t96 * t227 + (-t144 * t245 - t146 * t246 - t179 * t117 + t183 * t116 + m(4) * ((-t113 * t183 - t114 * t179) * qJD(3) + t313)) * t151 + (-t113 * t245 - t114 * t246 + t313) * mrSges(4,3) + (Ifges(4,1) * t134 + Ifges(4,4) * t320) * t179 - (-m(5) * t16 + m(6) * t12 + t319) * t202 + t134 * t275 / 0.2e1 + (t197 + t208 * qJD(3) / 0.2e1) * qJD(3) - t123 * (Ifges(5,4) * t96 - Ifges(5,2) * t97) / 0.2e1 + t125 * (mrSges(5,1) * t97 + mrSges(5,2) * t96) + t140 * (-mrSges(5,1) * t74 + mrSges(5,2) * t73) + (-mrSges(2,1) * t184 - t111 * mrSges(6,1) + mrSges(2,2) * t180 - t110 * mrSges(6,2) + t321 * (t163 * t159 - t162 * t185 + t170) + t322 * t170 + t306 * t162 + (-t308 - t329) * t163) * g(2) + t152 * (-mrSges(4,1) * t133 + mrSges(4,2) * t134) + t172 * (Ifges(5,5) * t96 - Ifges(5,6) * t97) / 0.2e1 + qJDD(3) * (Ifges(4,5) * t179 + Ifges(4,6) * t183) + t210 * t320 + t199 * t323 + (Ifges(5,1) * t96 - Ifges(5,4) * t97) * t293 + (-Ifges(6,1) * t198 - Ifges(6,4) * t199 + Ifges(6,5) * t97) * t297 + t122 * t245 / 0.2e1 - t121 * t246 / 0.2e1 - t8 * t257 / 0.2e1; m(3) * qJDD(2) + t179 * t116 + t183 * t117 - t328 * t97 + t319 * t131 + (-t144 * t179 + t146 * t183) * qJD(3) + t310 * t96 + (t60 + (-t177 * t75 - t181 * t76) * qJD(5) + t205) * t132 + m(6) * (t12 * t131 + t132 * t188 + t204 * t96 + t50 * t97) + m(5) * (-t131 * t16 + t132 * t15 - t54 * t97 + t55 * t96) + m(4) * (t179 * t84 + t183 * t85 + (-t113 * t179 + t114 * t183) * qJD(3)) + (t321 + t322) * g(3); -t64 * t106 - t84 * mrSges(4,2) + t85 * mrSges(4,1) - t23 * t75 - t22 * t76 + (t233 + t144) * t114 + ((t15 * t178 + t16 * t182 + (-t178 * t54 + t182 * t55) * qJD(4)) * pkin(3) - t125 * t237 + t54 * t63 - t55 * t64) * m(5) + t187 + (m(6) * t188 + t326) * (pkin(8) + t290) + (m(5) * t289 + mrSges(5,1) * t167 + t214 + t280) * (g(1) * t163 + g(2) * t162) + (t147 - m(6) * (t169 + t316) - m(5) * t169 + t315) * g(3) - (-Ifges(4,2) * t249 + t122 + t160) * t248 / 0.2e1 + (t232 - t146) * t113 + t309 * mrSges(6,3) + t328 * (-pkin(3) * t244 + t63) + (-t20 * t22 - t21 * t23 - t50 * t63 + t12 * t158 + (t178 * t50 + t182 * t204) * qJD(4) * pkin(3)) * m(6) + t310 * pkin(3) * t243 + (-t197 - t196 * qJD(1) / 0.2e1) * qJD(1) + Ifges(4,6) * t133 + Ifges(4,5) * t134 - g(1) * (t163 * t190 - t218) - g(2) * (t162 * t190 - t219) + t158 * t13 + Ifges(4,3) * qJDD(3) + t55 * t277 + t59 * t288 + t60 * t290 - t92 * t237 - t208 * t240 / 0.2e1 + t121 * t249 / 0.2e1; t191 * mrSges(6,3) + (t163 * t189 + t218) * g(1) - t54 * t106 - t25 * t75 - t24 * t76 - pkin(4) * t13 + (t162 * t189 + t219) * g(2) - m(6) * (t20 * t24 + t21 * t25 + t50 * t55) + t187 + (m(6) * (t284 + t309) + t326) * pkin(8) + (t328 + t277) * t55 + (t315 - t329) * g(3) - t12 * t304; -t50 * (mrSges(6,1) * t102 + mrSges(6,2) * t101) + (Ifges(6,1) * t101 - t272) * t298 + t47 * t297 + (Ifges(6,5) * t101 - Ifges(6,6) * t102) * t296 - t20 * t75 + t21 * t76 - g(1) * (mrSges(6,1) * t110 - mrSges(6,2) * t111) - g(2) * (-mrSges(6,1) * t108 + mrSges(6,2) * t109) + g(3) * t212 * t167 + (t101 * t20 + t102 * t21) * mrSges(6,3) + t239 + (-Ifges(6,2) * t102 + t48 + t98) * t299 + t307;];
tau = t1;
