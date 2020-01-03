% Calculate vector of inverse dynamics joint torques for
% S5RPRRR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR12_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR12_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:29
% EndTime: 2019-12-31 19:12:43
% DurationCPUTime: 6.93s
% Computational Cost: add. (4821->465), mult. (9248->644), div. (0->0), fcn. (5756->10), ass. (0->225)
t187 = cos(qJ(3));
t348 = t187 / 0.2e1;
t181 = sin(qJ(5));
t185 = cos(qJ(5));
t178 = qJD(3) + qJD(4);
t190 = -pkin(1) - pkin(6);
t147 = qJD(1) * t190 + qJD(2);
t133 = t187 * t147;
t244 = t187 * qJD(1);
t107 = -pkin(7) * t244 + t133;
t105 = qJD(3) * pkin(3) + t107;
t182 = sin(qJ(4));
t183 = sin(qJ(3));
t106 = (-pkin(7) * qJD(1) + t147) * t183;
t186 = cos(qJ(4));
t259 = t186 * t106;
t69 = t105 * t182 + t259;
t60 = pkin(8) * t178 + t69;
t125 = t182 * t187 + t183 * t186;
t119 = t125 * qJD(1);
t253 = qJD(1) * t183;
t120 = -t182 * t253 + t186 * t244;
t138 = pkin(3) * t253 + qJD(1) * qJ(2);
t71 = pkin(4) * t119 - pkin(8) * t120 + t138;
t20 = -t181 * t60 + t185 * t71;
t347 = t20 * mrSges(6,1);
t21 = t181 * t71 + t185 * t60;
t346 = t21 * mrSges(6,2);
t284 = t120 * mrSges(5,3);
t96 = -t120 * t181 + t178 * t185;
t97 = t120 * t185 + t178 * t181;
t337 = mrSges(5,1) * t178 + mrSges(6,1) * t96 - mrSges(6,2) * t97 - t284;
t180 = qJ(3) + qJ(4);
t171 = sin(t180);
t172 = cos(t180);
t345 = t171 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t172;
t290 = mrSges(6,2) * t181;
t292 = mrSges(6,1) * t185;
t344 = t290 - t292;
t248 = qJD(4) * t186;
t249 = qJD(4) * t182;
t242 = qJD(1) * qJD(3);
t129 = qJDD(1) * t187 - t183 * t242;
t146 = qJDD(1) * t190 + qJDD(2);
t252 = qJD(3) * t183;
t98 = t187 * t146 - t147 * t252;
t75 = qJDD(3) * pkin(3) - pkin(7) * t129 + t98;
t130 = -qJDD(1) * t183 - t187 * t242;
t251 = qJD(3) * t187;
t99 = t183 * t146 + t147 * t251;
t83 = pkin(7) * t130 + t99;
t17 = t105 * t248 - t106 * t249 + t182 * t75 + t186 * t83;
t177 = qJDD(3) + qJDD(4);
t14 = pkin(8) * t177 + t17;
t243 = qJD(1) * qJD(2);
t148 = qJDD(1) * qJ(2) + t243;
t103 = -pkin(3) * t130 + t148;
t64 = -qJD(4) * t119 + t129 * t186 + t130 * t182;
t258 = t186 * t187;
t124 = t182 * t183 - t258;
t65 = qJD(1) * qJD(4) * t124 - t129 * t182 + t130 * t186;
t19 = -pkin(4) * t65 - pkin(8) * t64 + t103;
t2 = qJD(5) * t20 + t14 * t185 + t181 * t19;
t3 = -qJD(5) * t21 - t14 * t181 + t185 * t19;
t343 = -t181 * t3 + t185 * t2;
t115 = qJD(5) + t119;
t296 = t97 * Ifges(6,4);
t41 = t96 * Ifges(6,2) + t115 * Ifges(6,6) + t296;
t342 = -t41 / 0.2e1;
t95 = Ifges(6,4) * t96;
t42 = t97 * Ifges(6,1) + t115 * Ifges(6,5) + t95;
t341 = t42 / 0.2e1;
t311 = -m(3) - m(4);
t340 = -m(6) - m(5);
t225 = -pkin(4) * t171 + t172 * pkin(8);
t339 = m(6) * t225;
t217 = mrSges(6,1) * t181 + mrSges(6,2) * t185;
t272 = t106 * t182;
t68 = t105 * t186 - t272;
t59 = -pkin(4) * t178 - t68;
t338 = t217 * t59;
t188 = cos(qJ(1));
t237 = t172 * t290;
t291 = mrSges(5,2) * t171;
t335 = (-t237 - t291) * t188;
t184 = sin(qJ(1));
t233 = t172 * t292;
t266 = t172 * t184;
t267 = t171 * t184;
t334 = -mrSges(5,1) * t266 - mrSges(6,3) * t267 - (t233 - t237) * t184;
t333 = -t171 * t344 + t345;
t206 = t183 * t99 + t187 * t98;
t31 = qJD(5) * t96 + t177 * t181 + t185 * t64;
t63 = qJDD(5) - t65;
t12 = mrSges(6,1) * t63 - mrSges(6,3) * t31;
t32 = -qJD(5) * t97 + t177 * t185 - t181 * t64;
t13 = -mrSges(6,2) * t63 + mrSges(6,3) * t32;
t332 = -t181 * t12 + t185 * t13;
t331 = mrSges(5,1) * t172 + t171 * mrSges(6,3) + t233;
t330 = -g(1) * t184 + g(2) * t188;
t220 = mrSges(4,1) * t187 - mrSges(4,2) * t183;
t287 = Ifges(4,4) * t187;
t329 = qJ(2) * t220 + (-Ifges(4,1) * t183 - t287) * t348;
t289 = mrSges(5,3) * t119;
t101 = -mrSges(5,2) * t178 - t289;
t66 = -mrSges(6,2) * t115 + mrSges(6,3) * t96;
t67 = mrSges(6,1) * t115 - mrSges(6,3) * t97;
t328 = -t181 * t67 + t185 * t66 + t101;
t10 = -mrSges(6,1) * t32 + mrSges(6,2) * t31;
t18 = -qJD(4) * t69 - t182 * t83 + t186 * t75;
t15 = -pkin(4) * t177 - t18;
t327 = -m(6) * t15 - t10;
t326 = mrSges(3,2) - mrSges(2,1) - mrSges(4,3) - mrSges(5,3);
t136 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t253;
t262 = t183 * (qJD(3) * mrSges(4,1) - mrSges(4,3) * t244);
t325 = (t187 * t136 - t262) * qJD(3);
t91 = -t182 * t251 - t183 * t248 - t186 * t252 - t187 * t249;
t92 = t178 * t258 - t182 * t252 - t183 * t249;
t324 = t124 * t18 - t125 * t17 - t68 * t91 - t69 * t92;
t323 = -m(6) * t59 + t337;
t322 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t321 = m(4) * t206 + t187 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t129) + t183 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t130);
t219 = mrSges(4,1) * t183 + mrSges(4,2) * t187;
t320 = t339 + mrSges(2,2) - t219 - mrSges(3,3) - t345;
t210 = t21 * t181 + t20 * t185;
t193 = -qJD(5) * t210 + t343;
t246 = qJD(5) * t185;
t247 = qJD(5) * t181;
t319 = m(6) * t193 - t67 * t246 - t66 * t247 + t332;
t318 = qJD(1) ^ 2;
t317 = t31 / 0.2e1;
t316 = t32 / 0.2e1;
t315 = t63 / 0.2e1;
t314 = -t96 / 0.2e1;
t313 = -t97 / 0.2e1;
t312 = t97 / 0.2e1;
t310 = -t115 / 0.2e1;
t307 = t120 / 0.2e1;
t306 = pkin(3) * t182;
t305 = pkin(3) * t186;
t304 = pkin(3) * t187;
t175 = t183 * pkin(3);
t294 = pkin(7) - t190;
t288 = Ifges(4,4) * t183;
t286 = Ifges(6,4) * t181;
t285 = Ifges(6,4) * t185;
t283 = t120 * Ifges(5,4);
t282 = t124 * t15;
t274 = t185 * t91;
t271 = t119 * t181;
t270 = t119 * t185;
t269 = t124 * t181;
t268 = t124 * t185;
t265 = t181 * t184;
t264 = t181 * t188;
t261 = t184 * t185;
t260 = t185 * t188;
t160 = qJ(2) + t175;
t255 = pkin(4) * t266 + pkin(8) * t267;
t254 = t188 * pkin(1) + t184 * qJ(2);
t245 = qJDD(1) * mrSges(3,2);
t149 = pkin(3) * t251 + qJD(2);
t239 = m(5) * t304;
t238 = Ifges(6,5) * t31 + Ifges(6,6) * t32 + Ifges(6,3) * t63;
t236 = pkin(3) * t244;
t229 = t181 * t342;
t135 = t294 * t187;
t227 = t246 / 0.2e1;
t174 = t188 * qJ(2);
t226 = -pkin(1) * t184 + t174;
t224 = -t242 / 0.2e1;
t223 = (t148 + t243) * qJ(2);
t222 = t294 * t252;
t88 = pkin(4) * t120 + pkin(8) * t119;
t221 = -pkin(4) * t172 - pkin(8) * t171;
t216 = t187 * Ifges(4,1) - t288;
t215 = Ifges(6,1) * t185 - t286;
t214 = -t183 * Ifges(4,2) + t287;
t213 = -Ifges(6,2) * t181 + t285;
t212 = -Ifges(4,5) * t183 - Ifges(4,6) * t187;
t211 = Ifges(6,5) * t185 - Ifges(6,6) * t181;
t209 = -t20 * t181 + t21 * t185;
t208 = -t181 * t66 - t185 * t67;
t89 = pkin(4) * t125 + pkin(8) * t124 + t160;
t134 = t294 * t183;
t94 = -t134 * t186 - t135 * t182;
t44 = -t181 * t94 + t185 * t89;
t45 = t181 * t89 + t185 * t94;
t93 = -t134 * t182 + t186 * t135;
t200 = t124 * t246 - t181 * t91;
t199 = t124 * t247 + t274;
t197 = t183 * (-Ifges(4,2) * t187 - t288);
t114 = Ifges(5,4) * t119;
t40 = t97 * Ifges(6,5) + t96 * Ifges(6,6) + t115 * Ifges(6,3);
t6 = t31 * Ifges(6,4) + t32 * Ifges(6,2) + t63 * Ifges(6,6);
t7 = t31 * Ifges(6,1) + t32 * Ifges(6,4) + t63 * Ifges(6,5);
t80 = -t119 * Ifges(5,2) + t178 * Ifges(5,6) + t283;
t81 = t120 * Ifges(5,1) + t178 * Ifges(5,5) - t114;
t192 = Ifges(5,5) * t64 + Ifges(5,6) * t65 - t17 * mrSges(5,2) + t18 * mrSges(5,1) - t138 * (mrSges(5,1) * t120 - mrSges(5,2) * t119) + (Ifges(6,3) * t120 - t119 * t211) * t310 + (Ifges(6,5) * t120 - t119 * t215) * t313 + (Ifges(6,6) * t120 - t119 * t213) * t314 - t178 * (-Ifges(5,5) * t119 - Ifges(5,6) * t120) / 0.2e1 - t120 * t347 + t69 * t284 + t270 * t341 + t271 * t342 - t68 * t289 + t119 * t338 + ((-t247 - t271) * t21 + (-t246 - t270) * t20 + t343) * mrSges(6,3) + t15 * t344 + (t115 * t211 + t213 * t96 + t215 * t97) * qJD(5) / 0.2e1 - (-Ifges(5,1) * t119 - t283 + t40) * t120 / 0.2e1 + (-Ifges(5,2) * t120 - t114 + t81) * t119 / 0.2e1 + t80 * t307 + t120 * t346 + (Ifges(6,5) * t181 + Ifges(6,6) * t185) * t315 + (t338 + t229) * qJD(5) + Ifges(5,3) * t177 + t181 * t7 / 0.2e1 + t185 * t6 / 0.2e1 + t42 * t227 + (Ifges(6,2) * t185 + t286) * t316 + (Ifges(6,1) * t181 + t285) * t317;
t189 = -pkin(7) - pkin(6);
t169 = -qJDD(1) * pkin(1) + qJDD(2);
t164 = -pkin(4) - t305;
t126 = t219 * qJD(1);
t122 = qJD(3) * t135;
t118 = Ifges(4,5) * qJD(3) + qJD(1) * t216;
t117 = Ifges(4,6) * qJD(3) + qJD(1) * t214;
t113 = t171 * t260 - t265;
t112 = t171 * t264 + t261;
t111 = t171 * t261 + t264;
t110 = -t171 * t265 + t260;
t87 = mrSges(5,1) * t119 + mrSges(5,2) * t120;
t79 = t88 + t236;
t73 = t107 * t186 - t272;
t72 = t107 * t182 + t259;
t55 = -mrSges(5,2) * t177 + mrSges(5,3) * t65;
t54 = mrSges(5,1) * t177 - mrSges(5,3) * t64;
t46 = -qJD(4) * t93 - t186 * t122 + t182 * t222;
t43 = pkin(4) * t92 - pkin(8) * t91 + t149;
t34 = t181 * t88 + t185 * t68;
t33 = -t181 * t68 + t185 * t88;
t28 = t181 * t79 + t185 * t73;
t27 = -t181 * t73 + t185 * t79;
t9 = -qJD(5) * t45 - t181 * t46 + t185 * t43;
t8 = qJD(5) * t44 + t181 * t43 + t185 * t46;
t1 = [t329 * t242 + t115 * (Ifges(6,5) * t199 + Ifges(6,6) * t200 + Ifges(6,3) * t92) / 0.2e1 + t96 * (Ifges(6,4) * t199 + Ifges(6,2) * t200 + Ifges(6,6) * t92) / 0.2e1 + t59 * (-mrSges(6,1) * t200 + mrSges(6,2) * t199) + t94 * t55 + t91 * t81 / 0.2e1 + t92 * t40 / 0.2e1 - t92 * t80 / 0.2e1 + t9 * t67 + t8 * t66 + t44 * t12 + t45 * t13 + (-m(3) * t226 - m(4) * t174 - t113 * mrSges(6,1) + t112 * mrSges(6,2) + t340 * (t188 * t175 + t184 * t189 + t226) + (-m(4) * t190 - t326) * t184 + t320 * t188) * g(1) + (-t111 * mrSges(6,1) - t110 * mrSges(6,2) + t311 * t254 + t340 * (t184 * t175 - t188 * t189 + t254) + (-m(4) * pkin(6) + t326) * t188 + t320 * t184) * g(2) - t183 * (Ifges(4,4) * t129 + Ifges(4,2) * t130) / 0.2e1 - t206 * mrSges(4,3) + (qJD(5) * t42 + t6) * t269 / 0.2e1 + (-m(5) * t18 - t327 - t54) * t93 + t130 * t214 / 0.2e1 + t129 * t216 / 0.2e1 + qJD(3) ^ 2 * t212 / 0.2e1 - t92 * t346 + t274 * t341 + (Ifges(4,1) * t129 + Ifges(4,4) * t130) * t348 + t324 * mrSges(5,3) + t321 * t190 + (Ifges(6,3) * t315 - Ifges(5,4) * t64 - Ifges(5,2) * t65 + t103 * mrSges(5,1) - Ifges(5,6) * t177 + Ifges(6,6) * t316 + Ifges(6,5) * t317 + t238 / 0.2e1 + t322) * t125 + (t219 + 0.2e1 * mrSges(3,3)) * t148 + (-m(5) * t68 - t323) * (qJD(4) * t94 - t122 * t182 - t186 * t222) + m(5) * (t103 * t160 + t138 * t149 + t17 * t94 + t46 * t69) + m(6) * (t2 * t45 + t20 * t9 + t21 * t8 + t3 * t44) + (Ifges(5,1) * t91 - Ifges(5,4) * t92) * t307 + t92 * t347 + t190 * t325 + (-t103 * mrSges(5,2) - Ifges(5,1) * t64 - Ifges(5,4) * t65 - Ifges(5,5) * t177 - t211 * t315 - t213 * t316 - t215 * t317 + t227 * t41) * t124 + (-t199 * t20 + t2 * t269 + t200 * t21 + t268 * t3) * mrSges(6,3) + m(4) * t223 + (Ifges(6,1) * t199 + Ifges(6,4) * t200 + Ifges(6,5) * t92) * t312 + t46 * t101 - t119 * (Ifges(5,4) * t91 - Ifges(5,2) * t92) / 0.2e1 + qJD(2) * t126 + qJ(2) * (-mrSges(4,1) * t130 + mrSges(4,2) * t129) + t138 * (mrSges(5,1) * t92 + mrSges(5,2) * t91) + t149 * t87 + t160 * (-mrSges(5,1) * t65 + mrSges(5,2) * t64) + t169 * mrSges(3,2) + t178 * (Ifges(5,5) * t91 - Ifges(5,6) * t92) / 0.2e1 + qJDD(3) * (Ifges(4,5) * t187 - Ifges(4,6) * t183) + t197 * t224 + t91 * t229 + m(3) * (-pkin(1) * t169 + t223) - pkin(1) * t245 - t117 * t251 / 0.2e1 - t118 * t252 / 0.2e1 - t7 * t268 / 0.2e1 - t217 * t282 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1); t245 + t337 * t91 + (t10 - t54) * t124 + t325 + t328 * t92 + (qJ(2) * t311 - mrSges(3,3)) * t318 + (qJD(5) * t208 + t332 + t55) * t125 + m(3) * t169 - m(5) * t324 + m(6) * (t125 * t193 + t209 * t92 - t59 * t91 + t282) + (-m(5) * t138 - m(6) * t210 - t126 + t208 - t87) * qJD(1) + t330 * (-t311 - t340) + t321; (t197 / 0.2e1 - t329) * t318 + t330 * t220 - t27 * t67 - t28 * t66 + t328 * pkin(3) * t248 + t337 * (-pkin(3) * t249 + t72) + (-(t239 - t291) * t184 - m(6) * (t184 * t304 + t255) + t334) * g(1) + ((t239 - m(6) * (t221 - t304) + t331) * t188 + t335) * g(2) + t319 * (pkin(8) + t306) + t147 * t262 + t192 + (t219 - m(6) * (-t175 + t225) + m(5) * t175 + t333) * g(3) + t54 * t305 + t55 * t306 + (t15 * t164 + (t182 * t59 + t186 * t209) * qJD(4) * pkin(3) - t20 * t27 - t21 * t28 - t59 * t72) * m(6) + t98 * mrSges(4,1) - t99 * mrSges(4,2) - t73 * t101 + ((t17 * t182 + t18 * t186 + (-t182 * t68 + t186 * t69) * qJD(4)) * pkin(3) - t138 * t236 + t68 * t72 - t69 * t73) * m(5) + Ifges(4,5) * t129 + Ifges(4,6) * t130 + t164 * t10 + t212 * t224 + Ifges(4,3) * qJDD(3) - t87 * t236 + t117 * t244 / 0.2e1 + t118 * t253 / 0.2e1 - t136 * t133; -t33 * t67 - t34 * t66 - m(6) * (t20 * t33 + t21 * t34) + t192 - t68 * t101 + t323 * t69 + (t333 - t339) * g(3) + ((-m(6) * t221 + t331) * t188 + t335) * g(2) + (-m(6) * t255 + mrSges(5,2) * t267 + t334) * g(1) + t327 * pkin(4) + t319 * pkin(8); -t59 * (mrSges(6,1) * t97 + mrSges(6,2) * t96) + (Ifges(6,1) * t96 - t296) * t313 + t41 * t312 + (Ifges(6,5) * t96 - Ifges(6,6) * t97) * t310 - t20 * t66 + t21 * t67 - g(1) * (mrSges(6,1) * t110 - mrSges(6,2) * t111) - g(2) * (mrSges(6,1) * t112 + mrSges(6,2) * t113) + g(3) * t217 * t172 + (t20 * t96 + t21 * t97) * mrSges(6,3) + t238 + (-Ifges(6,2) * t97 + t42 + t95) * t314 + t322;];
tau = t1;
