% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:41:16
% EndTime: 2019-03-09 02:41:37
% DurationCPUTime: 13.91s
% Computational Cost: add. (5199->587), mult. (10885->731), div. (0->0), fcn. (7143->14), ass. (0->273)
t294 = m(6) + m(7);
t234 = m(5) + t294;
t331 = -mrSges(5,2) + mrSges(6,3);
t330 = -mrSges(6,2) + mrSges(5,1);
t166 = qJ(1) + pkin(9);
t157 = sin(t166);
t280 = g(2) * t157;
t324 = Ifges(6,5) - Ifges(5,6);
t169 = sin(pkin(9));
t148 = pkin(1) * t169 + pkin(7);
t134 = t148 * qJDD(1);
t329 = qJD(2) * qJD(3) + t134;
t173 = sin(qJ(3));
t176 = cos(qJ(3));
t140 = -mrSges(4,1) * t176 + mrSges(4,2) * t173;
t165 = qJ(3) + pkin(10);
t156 = sin(t165);
t158 = cos(t165);
t328 = -t156 * t331 - t330 * t158 + t140;
t159 = cos(t166);
t310 = g(1) * t159 + t280;
t172 = sin(qJ(6));
t175 = cos(qJ(6));
t233 = qJD(1) * qJD(3);
t126 = qJDD(1) * t176 - t173 * t233;
t127 = qJDD(1) * t173 + t176 * t233;
t168 = sin(pkin(10));
t255 = cos(pkin(10));
t81 = -t126 * t255 + t127 * t168;
t219 = t255 * t176;
t241 = qJD(1) * t173;
t114 = -qJD(1) * t219 + t168 * t241;
t91 = -qJD(3) * t172 + t114 * t175;
t37 = qJD(6) * t91 + qJDD(3) * t175 + t172 * t81;
t300 = t37 / 0.2e1;
t92 = qJD(3) * t175 + t114 * t172;
t38 = -qJD(6) * t92 - qJDD(3) * t172 + t175 * t81;
t299 = t38 / 0.2e1;
t82 = t168 * t126 + t127 * t255;
t80 = qJDD(6) + t82;
t298 = t80 / 0.2e1;
t327 = m(3) + m(4);
t326 = t91 * Ifges(7,6);
t274 = qJD(3) / 0.2e1;
t325 = Ifges(5,5) - Ifges(6,4);
t109 = Ifges(5,4) * t114;
t124 = t168 * t176 + t173 * t255;
t116 = t124 * qJD(1);
t110 = qJD(6) + t116;
t321 = t110 * Ifges(7,3);
t323 = t116 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t92 * Ifges(7,5) - t109 + t321 + t326;
t71 = -qJDD(3) * mrSges(5,2) - mrSges(5,3) * t81;
t73 = mrSges(6,1) * t81 - qJDD(3) * mrSges(6,3);
t322 = t71 - t73;
t257 = qJDD(3) / 0.2e1;
t59 = -mrSges(7,2) * t110 + mrSges(7,3) * t91;
t60 = mrSges(7,1) * t110 - mrSges(7,3) * t92;
t197 = -t172 * t60 + t175 * t59;
t22 = mrSges(7,1) * t80 - mrSges(7,3) * t37;
t23 = -mrSges(7,2) * t80 + mrSges(7,3) * t38;
t312 = t172 * t23 + t175 * t22;
t74 = t82 * mrSges(6,1) + qJDD(3) * mrSges(6,2);
t320 = qJD(6) * t197 + t312 + t74;
t272 = mrSges(6,1) * t114;
t103 = -qJD(3) * mrSges(6,3) + t272;
t53 = -mrSges(7,1) * t91 + mrSges(7,2) * t92;
t319 = t53 - t103;
t318 = t310 * t156;
t270 = mrSges(5,3) * t114;
t101 = -qJD(3) * mrSges(5,2) - t270;
t317 = -t101 + t103;
t269 = mrSges(5,3) * t116;
t271 = mrSges(6,1) * t116;
t316 = -qJD(3) * t330 + t269 + t271;
t144 = t156 * qJ(5);
t246 = t158 * t159;
t315 = pkin(4) * t246 + t159 * t144;
t170 = cos(pkin(9));
t150 = -pkin(1) * t170 - pkin(2);
t163 = t176 * pkin(3);
t133 = t150 - t163;
t136 = t148 * qJD(1);
t238 = qJD(3) * t173;
t67 = t173 * qJDD(2) - t136 * t238 + t176 * t329;
t239 = qJD(2) * t173;
t106 = t136 * t176 + t239;
t161 = t176 * qJDD(2);
t68 = -t106 * qJD(3) - t134 * t173 + t161;
t311 = -t173 * t68 + t176 * t67;
t309 = 0.2e1 * t257;
t206 = mrSges(7,1) * t175 - mrSges(7,2) * t172;
t282 = pkin(5) * t114;
t214 = qJ(4) * qJD(1) + t136;
t94 = t176 * t214 + t239;
t222 = t255 * t94;
t162 = t176 * qJD(2);
t93 = -t173 * t214 + t162;
t87 = qJD(3) * pkin(3) + t93;
t47 = t168 * t87 + t222;
t41 = -qJD(3) * qJ(5) - t47;
t25 = -t41 - t282;
t276 = t92 * Ifges(7,4);
t30 = t91 * Ifges(7,2) + t110 * Ifges(7,6) + t276;
t308 = t25 * t206 - t175 * t30 / 0.2e1;
t307 = -m(4) * pkin(2) - mrSges(3,1) + t328;
t135 = t150 * qJDD(1);
t90 = -pkin(3) * t126 + qJDD(4) + t135;
t180 = -qJ(5) * t82 - qJD(5) * t116 + t90;
t293 = pkin(4) + pkin(8);
t11 = t293 * t81 + t180;
t232 = qJD(1) * qJD(4);
t237 = qJD(3) * t176;
t44 = -t136 * t237 + qJDD(3) * pkin(3) - qJ(4) * t127 + t161 + (-t232 - t329) * t173;
t48 = qJ(4) * t126 + t176 * t232 + t67;
t17 = -t168 * t48 + t255 * t44;
t194 = qJDD(5) - t17;
t7 = t82 * pkin(5) - qJDD(3) * t293 + t194;
t84 = t168 * t94;
t46 = t255 * t87 - t84;
t195 = qJD(5) - t46;
t277 = t116 * pkin(5);
t24 = -qJD(3) * t293 + t195 + t277;
t113 = qJD(1) * t133 + qJD(4);
t184 = -qJ(5) * t116 + t113;
t34 = t114 * t293 + t184;
t9 = -t172 * t34 + t175 * t24;
t1 = qJD(6) * t9 + t11 * t175 + t172 * t7;
t10 = t172 * t24 + t175 * t34;
t2 = -qJD(6) * t10 - t11 * t172 + t175 * t7;
t306 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t305 = (-g(1) * t246 - t158 * t280) * qJ(5);
t56 = pkin(4) * t114 + t184;
t304 = -t113 * mrSges(5,1) + t56 * mrSges(6,2);
t303 = -m(4) * pkin(7) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t302 = t9 * mrSges(7,1) + t113 * mrSges(5,2) - t10 * mrSges(7,2) - t56 * mrSges(6,3);
t301 = Ifges(7,1) * t300 + Ifges(7,4) * t299 + Ifges(7,5) * t298;
t297 = -t91 / 0.2e1;
t296 = -t92 / 0.2e1;
t295 = t92 / 0.2e1;
t292 = -t110 / 0.2e1;
t291 = -t114 / 0.2e1;
t290 = t114 / 0.2e1;
t289 = -t116 / 0.2e1;
t288 = t116 / 0.2e1;
t174 = sin(qJ(1));
t285 = pkin(1) * t174;
t284 = pkin(3) * t168;
t279 = g(3) * t158;
t278 = t1 * t172;
t146 = t158 * pkin(4);
t177 = cos(qJ(1));
t164 = t177 * pkin(1);
t275 = -qJD(3) / 0.2e1;
t18 = t168 * t44 + t255 * t48;
t268 = mrSges(7,3) * t175;
t267 = Ifges(4,4) * t173;
t266 = Ifges(4,4) * t176;
t265 = Ifges(7,4) * t172;
t264 = Ifges(7,4) * t175;
t263 = t116 * Ifges(5,4);
t262 = t116 * Ifges(6,6);
t115 = t124 * qJD(3);
t253 = t115 * t172;
t252 = t115 * t175;
t251 = t116 * t172;
t123 = t168 * t173 - t219;
t250 = t123 * t172;
t249 = t123 * t175;
t248 = t157 * t172;
t247 = t157 * t175;
t245 = t159 * t172;
t244 = t159 * t175;
t243 = qJ(4) + t148;
t152 = t163 + pkin(2);
t242 = t159 * t152 + t164;
t240 = qJD(1) * t176;
t236 = qJD(6) * t172;
t235 = qJD(6) * t175;
t230 = Ifges(7,5) * t37 + Ifges(7,6) * t38 + Ifges(7,3) * t80;
t155 = pkin(3) * t238;
t154 = pkin(3) * t241;
t229 = t101 + t319;
t228 = mrSges(4,3) * t241;
t227 = mrSges(4,3) * t240;
t225 = t146 + t144 + t163;
t224 = t255 * pkin(3);
t220 = -t236 / 0.2e1;
t50 = t168 * t93 + t222;
t215 = qJD(3) * t243;
t95 = qJD(4) * t176 - t173 * t215;
t96 = -qJD(4) * t173 - t176 * t215;
t54 = t168 * t95 - t255 * t96;
t217 = qJ(5) * t114 + t154;
t216 = -t152 - t144;
t121 = t243 * t173;
t122 = t243 * t176;
t69 = t255 * t121 + t122 * t168;
t213 = -m(7) * t293 - mrSges(7,3);
t149 = -t224 - pkin(4);
t171 = -qJ(4) - pkin(7);
t211 = -t157 * t171 + t242;
t210 = t10 * t175 - t9 * t172;
t209 = t10 * t172 + t9 * t175;
t208 = mrSges(4,1) * t173 + mrSges(4,2) * t176;
t205 = mrSges(7,1) * t172 + mrSges(7,2) * t175;
t203 = Ifges(7,1) * t172 + t264;
t202 = t176 * Ifges(4,2) + t267;
t201 = Ifges(7,2) * t175 + t265;
t200 = Ifges(4,5) * t176 - Ifges(4,6) * t173;
t199 = Ifges(7,5) * t172 + Ifges(7,6) * t175;
t188 = -qJ(5) * t124 + t133;
t52 = t123 * t293 + t188;
t57 = pkin(5) * t124 + t69;
t20 = t172 * t57 + t175 * t52;
t19 = -t172 * t52 + t175 * t57;
t196 = -t172 * t59 - t175 * t60;
t51 = t255 * t93 - t84;
t55 = t168 * t96 + t255 * t95;
t193 = t123 * t235 + t253;
t192 = t123 * t236 - t252;
t117 = qJD(3) * t219 - t168 * t238;
t191 = -qJ(5) * t117 - qJD(5) * t124 + t155;
t190 = t150 * qJD(1) * t208;
t189 = t173 * (Ifges(4,1) * t176 - t267);
t70 = -t168 * t121 + t122 * t255;
t186 = -t196 + t316;
t15 = -qJDD(3) * qJ(5) - qJD(3) * qJD(5) - t18;
t182 = qJD(6) * t210 + t175 * t2 + t278;
t153 = Ifges(4,4) * t240;
t145 = qJ(5) + t284;
t139 = -qJD(3) * mrSges(4,2) + t227;
t137 = qJD(3) * mrSges(4,1) - t228;
t119 = Ifges(4,1) * t241 + Ifges(4,5) * qJD(3) + t153;
t118 = Ifges(4,6) * qJD(3) + qJD(1) * t202;
t112 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t127;
t111 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t126;
t108 = Ifges(6,6) * t114;
t105 = -t136 * t173 + t162;
t100 = -t156 * t248 + t244;
t99 = t156 * t247 + t245;
t98 = t156 * t245 + t247;
t97 = t156 * t244 - t248;
t88 = Ifges(7,4) * t91;
t79 = t82 * mrSges(6,3);
t78 = t82 * mrSges(5,2);
t76 = -mrSges(6,2) * t114 - mrSges(6,3) * t116;
t75 = mrSges(5,1) * t114 + mrSges(5,2) * t116;
t72 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t82;
t65 = -t114 * Ifges(5,2) + Ifges(5,6) * qJD(3) + t263;
t64 = Ifges(6,4) * qJD(3) - t116 * Ifges(6,2) + t108;
t63 = Ifges(6,5) * qJD(3) + t114 * Ifges(6,3) - t262;
t62 = pkin(4) * t123 + t188;
t61 = pkin(4) * t116 + t217;
t58 = -t123 * pkin(5) + t70;
t49 = pkin(4) * t115 + t191;
t45 = t116 * t293 + t217;
t40 = -qJD(3) * pkin(4) + t195;
t33 = -t115 * pkin(5) + t55;
t32 = pkin(5) * t117 + t54;
t31 = Ifges(7,1) * t92 + Ifges(7,5) * t110 + t88;
t28 = t51 - t277;
t27 = t50 - t282;
t26 = t115 * t293 + t191;
t21 = pkin(4) * t81 + t180;
t16 = -qJDD(3) * pkin(4) + t194;
t14 = -mrSges(7,1) * t38 + mrSges(7,2) * t37;
t13 = t172 * t27 + t175 * t45;
t12 = -t172 * t45 + t175 * t27;
t8 = -pkin(5) * t81 - t15;
t5 = t37 * Ifges(7,4) + t38 * Ifges(7,2) + t80 * Ifges(7,6);
t4 = -qJD(6) * t20 - t172 * t26 + t175 * t32;
t3 = qJD(6) * t19 + t172 * t32 + t175 * t26;
t6 = [t250 * t301 + (m(4) * ((-t105 * t176 - t106 * t173) * qJD(3) + t311) - t137 * t237 - t139 * t238 - t173 * t112 + t176 * t111) * t148 + (-t105 * t237 - t106 * t238 + t311) * mrSges(4,3) + (-m(5) * t46 + m(6) * t40 + t316) * t54 + (m(5) * t47 - m(6) * t41 - t317) * t55 + t62 * (-t81 * mrSges(6,2) - t79) + (Ifges(7,1) * t193 - Ifges(7,4) * t192) * t295 + t91 * (Ifges(7,4) * t193 - Ifges(7,2) * t192) / 0.2e1 + m(7) * (t1 * t20 + t10 * t3 + t19 * t2 + t25 * t33 + t4 * t9 + t58 * t8) + (-m(5) * t17 + m(6) * t16 - t72 + t74) * t69 + (m(5) * t18 - m(6) * t15 + t322) * t70 + (Ifges(6,6) * t289 + Ifges(6,3) * t290 - Ifges(5,2) * t291 - Ifges(5,4) * t288 + t41 * mrSges(6,1) - t47 * mrSges(5,3) + t63 / 0.2e1 - t65 / 0.2e1 + t324 * t274 - t304) * t115 + (t200 * t274 + t190) * qJD(3) + (mrSges(2,1) * t174 - t100 * mrSges(7,1) + mrSges(2,2) * t177 + t99 * mrSges(7,2) + t327 * t285 - t234 * (-t159 * t171 - t285) + (-m(7) * pkin(5) + t303) * t159 + (-m(6) * (t216 - t146) - m(7) * t216 - t158 * t213 + m(5) * t152 - t307) * t157) * g(1) + (-Ifges(5,4) * t81 + Ifges(5,5) * qJDD(3) + t230) * t124 / 0.2e1 + (t16 * mrSges(6,1) + t90 * mrSges(5,2) - t17 * mrSges(5,3) - t21 * mrSges(6,3) + Ifges(5,5) * t257 + Ifges(7,5) * t300 + Ifges(7,6) * t299 + Ifges(7,3) * t298 + Ifges(6,2) * t82 + (-Ifges(5,4) / 0.2e1 - Ifges(6,6)) * t81 - t309 * Ifges(6,4) + t306) * t124 + (Ifges(4,5) * t173 + Ifges(4,6) * t176) * t257 + (m(3) * (t169 ^ 2 + t170 ^ 2) * pkin(1) ^ 2 + Ifges(3,3) + Ifges(2,3) + 0.2e1 * (mrSges(3,1) * t170 - mrSges(3,2) * t169) * pkin(1)) * qJDD(1) + (t189 + t176 * (-Ifges(4,2) * t173 + t266)) * t233 / 0.2e1 + (qJD(6) * t31 + t5) * t249 / 0.2e1 + t75 * t155 + (m(4) * t150 + t140) * t135 + t127 * t266 / 0.2e1 + t173 * (Ifges(4,4) * t126 + Ifges(4,5) * qJDD(3)) / 0.2e1 + m(5) * (t113 * t155 + t133 * t90) + m(6) * (t21 * t62 + t49 * t56) + (-m(7) * (pkin(8) * t246 + t242 + t315) - t98 * mrSges(7,1) - t97 * mrSges(7,2) - mrSges(7,3) * t246 - mrSges(2,1) * t177 + mrSges(2,2) * t174 - m(6) * (t211 + t315) - m(5) * t211 - t327 * t164 + t307 * t159 + (-m(7) * (pkin(5) - t171) + t303) * t157) * g(2) + (t90 * mrSges(5,1) + t15 * mrSges(6,1) - t21 * mrSges(6,2) - t18 * mrSges(5,3) + t199 * t298 + t201 * t299 + t203 * t300 - t8 * t206 + t30 * t220 + (-Ifges(6,6) - Ifges(5,4)) * t82 + (Ifges(6,3) + Ifges(5,2)) * t81 + t324 * t309) * t123 + t110 * (Ifges(7,5) * t193 - Ifges(7,6) * t192) / 0.2e1 + t82 * Ifges(5,1) * t124 + t30 * t252 / 0.2e1 + t31 * t253 / 0.2e1 + t119 * t237 / 0.2e1 - t118 * t238 / 0.2e1 + (t1 * t249 - t10 * t192 - t193 * t9 - t2 * t250) * mrSges(7,3) + t19 * t22 + t20 * t23 + t33 * t53 + t58 * t14 + t3 * t59 + t4 * t60 + t49 * t76 + t127 * t173 * Ifges(4,1) + t133 * (t81 * mrSges(5,1) + t78) + t150 * (-mrSges(4,1) * t126 + mrSges(4,2) * t127) + t25 * (mrSges(7,1) * t192 + mrSges(7,2) * t193) + t176 * (Ifges(4,4) * t127 + Ifges(4,2) * t126 + Ifges(4,6) * qJDD(3)) / 0.2e1 + (t323 / 0.2e1 - Ifges(6,2) * t289 - Ifges(6,6) * t290 + Ifges(5,4) * t291 + Ifges(7,5) * t295 + Ifges(5,1) * t288 - t46 * mrSges(5,3) + t40 * mrSges(6,1) - t64 / 0.2e1 + t321 / 0.2e1 + t326 / 0.2e1 + t325 * t274 + t302) * t117 + t126 * t202 / 0.2e1; m(3) * qJDD(2) + t173 * t111 + t176 * t112 + (-t137 * t173 + t139 * t176) * qJD(3) + (t14 + t322) * t124 + t229 * t117 + t186 * t115 + (-t234 - t327) * g(3) + (-t72 + t320) * t123 + m(5) * (-t115 * t46 + t117 * t47 - t123 * t17 + t124 * t18) + m(6) * (t115 * t40 - t117 * t41 + t123 * t16 - t124 * t15) + m(7) * (t115 * t209 + t117 * t25 + t123 * t182 + t124 * t8) + m(4) * (t173 * t67 + t176 * t68 + (-t105 * t173 + t106 * t176) * qJD(3)); (Ifges(7,5) * t175 - Ifges(7,6) * t172) * t298 + (-Ifges(7,2) * t172 + t264) * t299 + (Ifges(7,1) * t175 - t265) * t300 + t175 * t301 + (t64 + t108) * t291 + (m(7) * t182 + t235 * t59 - t236 * t60 + t312) * (-pkin(8) + t149) - t316 * t50 + t317 * t51 + (t208 + t234 * pkin(3) * t173 + (-t205 - t331) * t158 + (m(6) * pkin(4) - t213 + t330) * t156) * t310 + t308 * qJD(6) + (-t190 - t189 * qJD(1) / 0.2e1) * qJD(1) + (-Ifges(5,2) * t116 - t109 + t323) * t290 + t324 * t81 + (Ifges(6,3) * t291 - t10 * t268 + t199 * t292 + t201 * t297 + t203 * t296 + t275 * t324 + t304 + t308) * t116 + t325 * t82 + (-Ifges(5,1) * t289 - Ifges(7,5) * t296 + Ifges(6,2) * t288 - Ifges(7,6) * t297 - Ifges(7,3) * t292 - t275 * t325 + t302) * t114 + (t14 - t73) * t145 + t47 * t269 + t40 * t272 + (-t113 * t154 + t46 * t50 - t47 * t51 + (t168 * t18 + t17 * t255) * pkin(3)) * m(5) + (-t145 * t15 + t149 * t16 - t40 * t50 - t56 * t61 + (-qJD(5) + t51) * t41 + t305) * m(6) + (-t10 * t13 - t12 * t9 + t145 * t8 + (qJD(5) - t28) * t25 + t305) * m(7) - (-Ifges(4,2) * t241 + t119 + t153) * t240 / 0.2e1 - (t110 * t199 + t201 * t91 + t203 * t92) * qJD(6) / 0.2e1 + t71 * t284 + (t227 - t139) * t105 + t72 * t224 + t319 * qJD(5) + (t228 + t137) * t106 + (-t10 * t235 - t278 + (t236 + t251) * t9) * mrSges(7,3) + (t63 - t263) * t289 + (Ifges(5,3) + Ifges(4,3) + Ifges(6,1)) * qJDD(3) + (t262 + t65) * t288 + (t220 - t251 / 0.2e1) * t31 - t41 * t271 - t46 * t270 - t2 * t268 + t118 * t241 / 0.2e1 - t200 * t233 / 0.2e1 - t75 * t154 - t15 * mrSges(6,3) + t16 * mrSges(6,2) + t17 * mrSges(5,1) - t18 * mrSges(5,2) + (-m(5) * t163 - m(7) * (pkin(8) * t158 + t225) - t158 * mrSges(7,3) - t205 * t156 - m(6) * t225 + t328) * g(3) - t28 * t53 - t13 * t59 - t12 * t60 - t67 * mrSges(4,2) + t68 * mrSges(4,1) - t61 * t76 + Ifges(4,6) * t126 + Ifges(4,5) * t127 + t149 * t74 - t172 * t5 / 0.2e1 + t8 * t205; -t172 * t22 + t175 * t23 + t78 - t79 + t330 * t81 + t196 * qJD(6) + t229 * t114 - t186 * t116 + (-g(1) * t157 + g(2) * t159) * t234 + (t1 * t175 - t110 * t209 + t114 * t25 - t2 * t172) * m(7) + (-t114 * t41 - t116 * t40 + t21) * m(6) + (t114 * t47 + t116 * t46 + t90) * m(5); -t319 * qJD(3) + t294 * t279 + (t197 + t76) * t116 + (-qJD(3) * t25 + t116 * t210 + t182 - t318) * m(7) + (qJD(3) * t41 + t116 * t56 + t16 - t318) * m(6) + t320; -t25 * (mrSges(7,1) * t92 + mrSges(7,2) * t91) + (Ifges(7,1) * t91 - t276) * t296 + t30 * t295 + (Ifges(7,5) * t91 - Ifges(7,6) * t92) * t292 - t9 * t59 + t10 * t60 - g(1) * (mrSges(7,1) * t97 - mrSges(7,2) * t98) - g(2) * (mrSges(7,1) * t99 + mrSges(7,2) * t100) + t206 * t279 + (t10 * t92 + t9 * t91) * mrSges(7,3) + t230 + (-Ifges(7,2) * t92 + t31 + t88) * t297 + t306;];
tau  = t6;
