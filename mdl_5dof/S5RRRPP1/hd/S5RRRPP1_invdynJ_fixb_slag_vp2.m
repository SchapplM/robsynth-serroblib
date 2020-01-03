% Calculate vector of inverse dynamics joint torques for
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:24
% EndTime: 2019-12-31 20:49:34
% DurationCPUTime: 5.22s
% Computational Cost: add. (3463->392), mult. (5231->499), div. (0->0), fcn. (3062->12), ass. (0->181)
t323 = mrSges(5,1) + mrSges(6,1);
t322 = mrSges(5,2) - mrSges(6,3);
t318 = Ifges(5,1) + Ifges(6,1);
t308 = Ifges(5,5) + Ifges(6,4);
t190 = sin(qJ(3));
t273 = Ifges(4,4) * t190;
t317 = Ifges(6,5) - Ifges(5,4);
t193 = cos(qJ(3));
t312 = t193 * mrSges(4,1) - t190 * mrSges(4,2);
t321 = -mrSges(3,1) - t312;
t184 = qJ(3) + pkin(8);
t175 = sin(t184);
t176 = cos(t184);
t320 = t322 * t175 - t323 * t176;
t319 = -mrSges(6,2) - mrSges(5,3) - mrSges(4,3) + mrSges(3,2);
t188 = sin(pkin(8));
t259 = cos(pkin(8));
t230 = t259 * t190;
t131 = t188 * t193 + t230;
t288 = t131 / 0.2e1;
t278 = qJD(3) / 0.2e1;
t316 = Ifges(5,6) - Ifges(6,6);
t183 = qJD(1) + qJD(2);
t252 = t188 * t190;
t203 = t259 * t193 - t252;
t111 = t203 * t183;
t104 = Ifges(5,4) * t111;
t112 = t131 * t183;
t270 = Ifges(6,5) * t111;
t315 = t308 * qJD(3) + t318 * t112 + t104 - t270;
t177 = t193 * qJD(4);
t189 = -qJ(4) - pkin(7);
t231 = qJD(3) * t189;
t119 = t190 * t231 + t177;
t202 = -qJD(4) * t190 + t193 * t231;
t194 = cos(qJ(2));
t269 = pkin(1) * qJD(1);
t242 = t194 * t269;
t314 = -t119 * t259 - t188 * t202 + t203 * t242;
t301 = t119 * t188 - t131 * t242 - t202 * t259;
t191 = sin(qJ(2));
t268 = pkin(1) * qJD(2);
t240 = qJD(1) * t268;
t258 = pkin(1) * qJDD(1);
t129 = t191 * t258 + t194 * t240;
t182 = qJDD(1) + qJDD(2);
t117 = pkin(7) * t182 + t129;
t243 = t191 * t269;
t137 = pkin(7) * t183 + t243;
t246 = qJD(3) * t190;
t71 = t193 * t117 - t137 * t246;
t245 = qJD(3) * t193;
t237 = t137 * t245;
t72 = -t117 * t190 - t237;
t313 = -t190 * t72 + t193 * t71;
t187 = qJ(1) + qJ(2);
t178 = sin(t187);
t179 = cos(t187);
t310 = g(1) * t179 + g(2) * t178;
t123 = t182 * t193 - t183 * t246;
t124 = t182 * t190 + t183 * t245;
t66 = -t123 * t259 + t124 * t188;
t52 = -qJDD(3) * mrSges(5,2) - mrSges(5,3) * t66;
t55 = -mrSges(6,2) * t66 + qJDD(3) * mrSges(6,3);
t306 = t52 + t55;
t67 = t188 * t123 + t124 * t259;
t53 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t67;
t54 = -qJDD(3) * mrSges(6,1) + t67 * mrSges(6,2);
t305 = t54 - t53;
t274 = mrSges(5,3) * t112;
t275 = mrSges(6,2) * t112;
t304 = qJD(3) * t323 - t274 - t275;
t266 = t111 * mrSges(5,3);
t267 = t111 * mrSges(6,2);
t94 = qJD(3) * mrSges(6,3) + t267;
t303 = qJD(3) * mrSges(5,2) - t266 - t94;
t128 = -t191 * t240 + t194 * t258;
t116 = -pkin(2) * t182 - t128;
t299 = m(4) * t116 - mrSges(4,1) * t123 + mrSges(4,2) * t124;
t282 = pkin(3) * t193;
t169 = pkin(2) + t282;
t213 = -pkin(4) * t176 - qJ(5) * t175;
t298 = t319 * t179 + (-m(6) * (-t169 + t213) - t320 - t321) * t178;
t256 = t176 * t179;
t257 = t175 * t179;
t297 = t319 * t178 + t321 * t179 - t256 * t323 + t322 * t257;
t254 = t183 * t190;
t135 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t254;
t138 = -pkin(2) * t183 - t242;
t285 = pkin(1) * t194;
t287 = pkin(1) * t191;
t296 = m(4) * pkin(1) * (t138 * t191 + (t190 ^ 2 + t193 ^ 2) * t194 * t137) - t190 * t135 * t285 + (-mrSges(3,1) * t287 - mrSges(3,2) * t285) * t183;
t253 = t183 * t193;
t136 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t253;
t295 = m(4) * t313 - t135 * t245 - t136 * t246 + t193 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t123) - t190 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t124);
t293 = t111 / 0.2e1;
t292 = -t111 / 0.2e1;
t290 = t112 / 0.2e1;
t192 = sin(qJ(1));
t286 = pkin(1) * t192;
t284 = pkin(3) * t188;
t195 = cos(qJ(1));
t181 = t195 * pkin(1);
t25 = -t237 + qJDD(3) * pkin(3) - qJ(4) * t124 + (-qJD(4) * t183 - t117) * t190;
t32 = qJ(4) * t123 + t177 * t183 + t71;
t8 = t188 * t25 + t259 * t32;
t224 = qJ(4) * t183 + t137;
t99 = t224 * t193;
t86 = t259 * t99;
t98 = t224 * t190;
t90 = qJD(3) * pkin(3) - t98;
t42 = t188 * t90 + t86;
t272 = Ifges(4,4) * t193;
t271 = Ifges(5,4) * t112;
t265 = t188 * t99;
t255 = t179 * t189;
t168 = pkin(7) + t287;
t249 = -qJ(4) - t168;
t248 = t179 * pkin(2) + t178 * pkin(7);
t247 = qJD(3) * t183;
t244 = pkin(3) * t254;
t241 = t194 * t268;
t173 = t191 * t268;
t172 = pkin(3) * t246;
t238 = t259 * pkin(3);
t15 = t66 * mrSges(5,1) + t67 * mrSges(5,2);
t14 = t66 * mrSges(6,1) - t67 * mrSges(6,3);
t232 = -pkin(2) * t178 + t179 * pkin(7);
t227 = t249 * t190;
t223 = t179 * t169 - t178 * t189;
t222 = qJD(3) * t249;
t220 = t193 * t241;
t219 = -g(1) * t178 + g(2) * t179;
t217 = mrSges(4,1) * t190 + mrSges(4,2) * t193;
t215 = t193 * Ifges(4,2) + t273;
t214 = Ifges(4,5) * t193 - Ifges(4,6) * t190;
t211 = -t169 * t178 - t255;
t206 = pkin(4) * t256 + qJ(5) * t257 + t223;
t7 = -t188 * t32 + t25 * t259;
t41 = t259 * t90 - t265;
t205 = t138 * t217;
t204 = t190 * (Ifges(4,1) * t193 - t273);
t77 = -pkin(4) * t203 - qJ(5) * t131 - t169;
t121 = t131 * qJD(3);
t122 = t203 * qJD(3);
t44 = pkin(4) * t121 - qJ(5) * t122 - qJD(5) * t131 + t172;
t108 = -t169 * t183 + qJD(4) - t242;
t73 = -pkin(3) * t123 + qJDD(4) + t116;
t199 = (-qJD(4) - t241) * t190 + t193 * t222;
t114 = Ifges(4,6) * qJD(3) + t183 * t215;
t152 = Ifges(4,4) * t253;
t115 = Ifges(4,1) * t254 + Ifges(4,5) * qJD(3) + t152;
t33 = -pkin(4) * t111 - qJ(5) * t112 + t108;
t35 = -qJD(3) * pkin(4) + qJD(5) - t41;
t36 = qJD(3) * qJ(5) + t42;
t4 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t8;
t103 = Ifges(6,5) * t112;
t48 = Ifges(6,6) * qJD(3) - t111 * Ifges(6,3) + t103;
t49 = t111 * Ifges(5,2) + Ifges(5,6) * qJD(3) + t271;
t5 = -qJDD(3) * pkin(4) + qJDD(5) - t7;
t6 = pkin(4) * t66 - qJ(5) * t67 - qJD(5) * t112 + t73;
t196 = (t308 * t131 + t316 * t203) * qJDD(3) / 0.2e1 - t203 * (Ifges(6,5) * t67 + Ifges(6,6) * qJDD(3)) / 0.2e1 + t203 * (Ifges(5,4) * t67 + Ifges(5,6) * qJDD(3)) / 0.2e1 + t73 * (-mrSges(5,1) * t203 + mrSges(5,2) * t131) + t6 * (-mrSges(6,1) * t203 - mrSges(6,3) * t131) + (t183 * (-Ifges(4,2) * t190 + t272) + t115) * t245 / 0.2e1 + t313 * mrSges(4,3) + (t273 + t215) * t123 / 0.2e1 + ((-Ifges(6,3) - Ifges(5,2)) * t203 + 0.2e1 * t317 * t288) * t66 + (-t41 * mrSges(5,3) + t35 * mrSges(6,2) + t315 / 0.2e1 + t108 * mrSges(5,2) - t33 * mrSges(6,3) + Ifges(5,4) * t293 + Ifges(6,5) * t292 + t308 * t278 + t318 * t290) * t122 + (-t42 * t121 - t131 * t7 + t203 * t8) * mrSges(5,3) + (-t121 * t36 + t131 * t5 + t203 * t4) * mrSges(6,2) + t124 * t272 / 0.2e1 + Ifges(3,3) * t182 + t128 * mrSges(3,1) - t129 * mrSges(3,2) + t193 * (Ifges(4,4) * t124 + Ifges(4,2) * t123) / 0.2e1 + t190 * Ifges(4,5) * qJDD(3) + t190 * t124 * Ifges(4,1) - t116 * t312 + (Ifges(6,3) * t292 - Ifges(5,2) * t293 - t49 / 0.2e1 + t108 * mrSges(5,1) + t33 * mrSges(6,1) + t48 / 0.2e1 + t317 * t290 - t316 * t278) * t121 + (t318 * t131 - t317 * t203) * t67 / 0.2e1 + (t308 * qJDD(3) + t318 * t67) * t288 + (t214 * t278 + t205) * qJD(3) + t193 * qJDD(3) * Ifges(4,6) - t114 * t246 / 0.2e1 + t204 * t247 / 0.2e1;
t180 = t193 * qJ(4);
t164 = -t238 - pkin(4);
t155 = qJ(5) + t284;
t149 = pkin(7) * t193 + t180;
t146 = -t169 - t285;
t134 = t173 + t172;
t127 = t168 * t193 + t180;
t120 = t312 * t183;
t89 = t149 * t259 + t189 * t252;
t88 = t149 * t188 - t189 * t230;
t81 = t190 * t222 + t177 + t220;
t68 = t77 - t285;
t57 = -mrSges(5,1) * t111 + mrSges(5,2) * t112;
t56 = -mrSges(6,1) * t111 - mrSges(6,3) * t112;
t47 = pkin(4) * t112 - qJ(5) * t111 + t244;
t46 = -t259 * t98 - t265;
t45 = -t188 * t98 + t86;
t34 = t173 + t44;
t1 = [m(5) * (t108 * t134 + t146 * t73) + m(6) * (t33 * t34 + t6 * t68) + t196 + m(3) * (t128 * t194 + t129 * t191) * pkin(1) + t134 * t57 + t146 * t15 + t68 * t14 + t34 * t56 - t120 * t173 + t136 * t220 + Ifges(2,3) * qJDD(1) + (m(5) * t8 + m(6) * t4 + t306) * (t127 * t259 + t188 * t227) + (-m(5) * t7 + m(6) * t5 + t305) * (t127 * t188 - t227 * t259) + (m(5) * t42 + m(6) * t36 - t303) * (t188 * t199 + t259 * t81) + (-m(5) * t41 + m(6) * t35 - t304) * (t188 * t81 - t199 * t259) + (mrSges(3,1) * t285 - mrSges(3,2) * t287) * t182 + t299 * (-pkin(2) - t285) + t296 * qJD(2) + t295 * t168 + (-mrSges(2,1) * t195 + t192 * mrSges(2,2) - m(5) * (t181 + t223) - m(6) * (t181 + t206) - m(4) * (t181 + t248) - m(3) * t181 + t297) * g(2) + (t192 * mrSges(2,1) + mrSges(2,2) * t195 + m(3) * t286 - m(5) * (t211 - t286) - m(4) * (t232 - t286) - m(6) * (-t255 - t286) + t298) * g(1); -t193 * t136 * t242 + t77 * t14 - t169 * t15 + t57 * t172 + t44 * t56 + t196 + t306 * t89 + t305 * t88 - t299 * pkin(2) + (t301 * t35 - t314 * t36 + t33 * t44 + t4 * t89 + t5 * t88 + t6 * t77) * m(6) + (t108 * t172 - t169 * t73 - t301 * t41 - t314 * t42 - t7 * t88 + t8 * t89) * m(5) + (-m(5) * t108 - m(6) * t33 + t120 - t56 - t57) * t243 - t296 * qJD(1) + (-m(4) * t248 - m(5) * t223 - m(6) * t206 + t297) * g(2) + (-m(4) * t232 - m(5) * t211 + m(6) * t255 + t298) * g(1) + t295 * pkin(7) - t301 * t304 + t314 * t303; (-Ifges(5,2) * t112 + t104 + t315) * t292 - t33 * (mrSges(6,1) * t112 - mrSges(6,3) * t111) - t108 * (mrSges(5,1) * t112 + mrSges(5,2) * t111) + (Ifges(6,3) * t112 + t270) * t293 + t304 * t45 + t308 * t67 + t303 * t46 + t310 * (t217 + (m(6) + m(5)) * pkin(3) * t190 + (-m(6) * qJ(5) + t322) * t176 + (m(6) * pkin(4) + t323) * t175) + t52 * t284 + t49 * t290 + (-t205 - t204 * t183 / 0.2e1) * t183 + t155 * t55 + t164 * t54 + Ifges(4,6) * t123 + Ifges(4,5) * t124 + ((t188 * t8 + t259 * t7) * pkin(3) - t108 * t244 + t41 * t45 - t42 * t46) * m(5) + qJD(5) * t94 - t71 * mrSges(4,2) + t72 * mrSges(4,1) - t47 * t56 + t4 * mrSges(6,3) - t5 * mrSges(6,1) + t7 * mrSges(5,1) - t8 * mrSges(5,2) - t35 * t267 + t42 * t274 + t36 * t275 + t53 * t238 + (-t312 - m(5) * t282 - m(6) * (-t213 + t282) + t320) * g(3) - (t111 * t308 - t112 * t316) * qJD(3) / 0.2e1 - t316 * t66 - (t318 * t111 + t103 - t271 + t48) * t112 / 0.2e1 + (t155 * t4 + t164 * t5 - t33 * t47 - t35 * t45 + (-t46 + qJD(5)) * t36) * m(6) + (Ifges(5,3) + Ifges(4,3) + Ifges(6,2)) * qJDD(3) - (-Ifges(4,2) * t254 + t115 + t152) * t253 / 0.2e1 - t57 * t244 - t214 * t247 / 0.2e1 + t114 * t254 / 0.2e1 + t41 * t266 + (t135 * t193 + t136 * t190) * t137; t304 * t112 + t303 * t111 + t14 + t15 + (-t111 * t36 - t112 * t35 + t219 + t6) * m(6) + (-t111 * t42 + t112 * t41 + t219 + t73) * m(5); -qJD(3) * t94 + t112 * t56 + (g(3) * t176 - t36 * qJD(3) + t33 * t112 - t175 * t310 + t5) * m(6) + t54;];
tau = t1;
