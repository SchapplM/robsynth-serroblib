% Calculate vector of inverse dynamics joint torques for
% S5RRPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:28
% EndTime: 2019-12-31 20:15:34
% DurationCPUTime: 2.82s
% Computational Cost: add. (3272->336), mult. (4252->437), div. (0->0), fcn. (2268->12), ass. (0->178)
t170 = sin(qJ(4));
t297 = -t170 / 0.2e1;
t174 = cos(qJ(4));
t272 = t174 / 0.2e1;
t290 = mrSges(3,1) - mrSges(4,2);
t164 = qJD(1) + qJD(2);
t169 = sin(qJ(5));
t173 = cos(qJ(5));
t233 = t173 * t174;
t282 = -t169 * t170 + t233;
t88 = t282 * t164;
t257 = mrSges(5,2) * t174;
t296 = -t257 + mrSges(3,2) - mrSges(4,3);
t295 = -mrSges(5,3) - mrSges(6,3) - t290;
t171 = sin(qJ(2));
t253 = pkin(1) * qJD(1);
t223 = t171 * t253;
t108 = qJ(3) * t164 + t223;
t200 = mrSges(5,1) * t174 - mrSges(5,2) * t170;
t254 = Ifges(5,4) * t174;
t255 = Ifges(5,4) * t170;
t294 = ((-Ifges(5,1) * t170 - t254) * t272 + (-Ifges(5,2) * t174 - t255) * t297) * t164 + t108 * t200 + qJD(4) * (-Ifges(5,5) * t170 - Ifges(5,6) * t174) / 0.2e1;
t190 = t169 * t174 + t173 * t170;
t236 = t164 * t170;
t178 = -pkin(2) - pkin(7);
t175 = cos(qJ(2));
t222 = t175 * t253;
t202 = qJD(3) - t222;
t90 = t164 * t178 + t202;
t57 = -pkin(8) * t236 + t170 * t90;
t250 = t169 * t57;
t235 = t164 * t174;
t58 = -pkin(8) * t235 + t174 * t90;
t54 = qJD(4) * pkin(4) + t58;
t25 = t173 * t54 - t250;
t247 = t173 * t57;
t26 = t169 * t54 + t247;
t227 = qJD(4) * t170;
t162 = qJDD(1) + qJDD(2);
t252 = pkin(1) * qJD(2);
t221 = t171 * t252;
t267 = pkin(1) * t175;
t98 = -qJD(1) * t221 + qJDD(1) * t267;
t186 = qJDD(3) - t98;
t67 = t162 * t178 + t186;
t39 = t174 * t67 - t227 * t90;
t92 = t162 * t174 - t164 * t227;
t19 = qJDD(4) * pkin(4) - pkin(8) * t92 + t39;
t226 = qJD(4) * t174;
t40 = t170 * t67 + t90 * t226;
t93 = -t162 * t170 - t164 * t226;
t23 = pkin(8) * t93 + t40;
t3 = qJD(5) * t25 + t169 * t19 + t173 * t23;
t286 = qJD(5) * t26;
t4 = -t169 * t23 + t173 * t19 - t286;
t163 = qJD(4) + qJD(5);
t225 = qJD(5) * t169;
t55 = t163 * t233 - t169 * t227 - t170 * t225;
t182 = t190 * qJD(5);
t56 = -qJD(4) * t190 - t182;
t293 = -t190 * t3 - t25 * t56 - t26 * t55 - t282 * t4;
t150 = pkin(8) * t227;
t214 = t178 * t227;
t100 = t150 - t214;
t258 = -pkin(8) + t178;
t117 = t258 * t174;
t101 = qJD(4) * t117;
t116 = t258 * t170;
t62 = t116 * t173 + t117 * t169;
t289 = -qJD(5) * t62 + t100 * t173 - t101 * t169 - t223 * t282;
t61 = -t116 * t169 + t117 * t173;
t288 = qJD(5) * t61 + t100 * t169 + t101 * t173 - t190 * t223;
t199 = mrSges(5,1) * t170 + t257;
t91 = t199 * t164;
t287 = mrSges(4,3) * t164 + t91;
t159 = t170 * pkin(4);
t144 = qJ(3) + t159;
t168 = qJ(1) + qJ(2);
t156 = sin(t168);
t177 = -pkin(8) - pkin(7);
t158 = cos(t168);
t237 = t158 * t170;
t285 = pkin(4) * t237 + t156 * t177;
t149 = -pkin(2) - t267;
t134 = -pkin(7) + t149;
t284 = -t134 * t227 + t174 * t221;
t283 = qJ(3) * t162 + qJD(3) * t164;
t194 = t170 * t40 + t174 * t39;
t280 = -g(1) * t156 + g(2) * t158;
t167 = qJ(4) + qJ(5);
t157 = cos(t167);
t238 = t157 * t158;
t155 = sin(t167);
t241 = t155 * t158;
t279 = -mrSges(5,1) * t237 - mrSges(6,1) * t241 - mrSges(6,2) * t238 + t296 * t158 + (-m(5) * t178 - t295) * t156;
t239 = t156 * t170;
t240 = t156 * t157;
t242 = t155 * t156;
t278 = -mrSges(5,1) * t239 - mrSges(6,1) * t242 - mrSges(6,2) * t240 + t296 * t156 + t295 * t158;
t277 = m(5) * t194 + t174 * (qJDD(4) * mrSges(5,1) - mrSges(5,3) * t92) + t170 * (-qJDD(4) * mrSges(5,2) + mrSges(5,3) * t93);
t107 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t236;
t276 = t107 * t226 + t277;
t274 = t88 / 0.2e1;
t273 = -m(4) - m(5);
t87 = t190 * t164;
t271 = mrSges(6,3) * t87;
t270 = Ifges(6,4) * t88;
t269 = pkin(1) * t171;
t172 = sin(qJ(1));
t268 = pkin(1) * t172;
t176 = cos(qJ(1));
t160 = t176 * pkin(1);
t260 = t88 * mrSges(6,3);
t259 = -pkin(8) + t134;
t251 = t162 * mrSges(4,2);
t243 = t108 * t175;
t231 = t158 * pkin(2) + t156 * qJ(3);
t230 = qJD(2) * t175;
t220 = pkin(1) * t230;
t217 = t160 + t231;
t47 = mrSges(6,1) * t87 + mrSges(6,2) * t88;
t89 = t144 * t164 + t223;
t216 = -m(6) * t89 - t47;
t97 = t259 * t174;
t143 = qJ(3) + t269;
t136 = t158 * qJ(3);
t212 = -pkin(2) * t156 + t136;
t208 = t164 * t222;
t206 = t170 * t221;
t204 = (t170 ^ 2 + t174 ^ 2) * t90 * t171;
t127 = qJD(3) + t220;
t198 = -mrSges(6,1) * t155 - mrSges(6,2) * t157;
t197 = t174 * Ifges(5,1) - t255;
t196 = -Ifges(5,2) * t170 + t254;
t96 = t259 * t170;
t50 = t169 * t97 + t173 * t96;
t49 = -t169 * t96 + t173 * t97;
t99 = (qJD(1) * t230 + qJDD(1) * t171) * pkin(1);
t72 = t99 + t283;
t193 = t108 * t127 + t143 * t72;
t192 = t72 * qJ(3) + t108 * qJD(3);
t109 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t235;
t191 = t107 * t174 - t109 * t170;
t188 = t212 - t268;
t187 = pkin(4) * t239 - t158 * t177 + t231;
t161 = qJDD(4) + qJDD(5);
t37 = -t164 * t182 + t169 * t93 + t173 * t92;
t38 = -qJD(5) * t88 - t169 * t92 + t173 * t93;
t42 = -Ifges(6,2) * t87 + Ifges(6,6) * t163 + t270;
t80 = Ifges(6,4) * t87;
t43 = Ifges(6,1) * t88 + Ifges(6,5) * t163 - t80;
t180 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t25 * t271 + t42 * t274 - t89 * (mrSges(6,1) * t88 - mrSges(6,2) * t87) + Ifges(6,3) * t161 - t88 * (-Ifges(6,1) * t87 - t270) / 0.2e1 + Ifges(6,6) * t38 + Ifges(6,5) * t37 - t163 * (-Ifges(6,5) * t87 - Ifges(6,6) * t88) / 0.2e1 + (-Ifges(6,2) * t88 + t43 - t80) * t87 / 0.2e1;
t48 = -pkin(4) * t93 + t72;
t82 = -pkin(2) * t162 + t186;
t85 = Ifges(5,6) * qJD(4) + t164 * t196;
t86 = Ifges(5,5) * qJD(4) + t164 * t197;
t179 = -t194 * mrSges(5,3) + (Ifges(5,4) * t92 + Ifges(5,2) * t93) * t297 + t163 * (Ifges(6,5) * t56 - Ifges(6,6) * t55) / 0.2e1 + t98 * mrSges(3,1) - t99 * mrSges(3,2) + t93 * t196 / 0.2e1 + t92 * t197 / 0.2e1 + t82 * mrSges(4,2) - t87 * (Ifges(6,4) * t56 - Ifges(6,2) * t55) / 0.2e1 + t89 * (mrSges(6,1) * t55 + mrSges(6,2) * t56) + t56 * t43 / 0.2e1 - t55 * t42 / 0.2e1 + (Ifges(5,1) * t92 + Ifges(5,4) * t93) * t272 + (Ifges(6,1) * t56 - Ifges(6,4) * t55) * t274 - t85 * t226 / 0.2e1 - t86 * t227 / 0.2e1 + (t199 + mrSges(4,3)) * t72 + (Ifges(4,1) + Ifges(3,3)) * t162 + (0.2e1 * Ifges(5,5) * t272 - Ifges(5,6) * t170) * qJDD(4) + t293 * mrSges(6,3) - (-t48 * mrSges(6,2) - Ifges(6,1) * t37 - Ifges(6,4) * t38 - Ifges(6,5) * t161) * t282 + (t48 * mrSges(6,1) - Ifges(6,4) * t37 - Ifges(6,2) * t38 - Ifges(6,6) * t161) * t190 + t294 * qJD(4);
t151 = pkin(4) * t226;
t145 = t158 * pkin(7);
t128 = qJD(3) + t151;
t118 = t143 + t159;
t113 = mrSges(6,2) * t241;
t112 = mrSges(6,1) * t240;
t104 = -pkin(2) * t164 + t202;
t102 = t127 + t151;
t69 = qJD(4) * t97 + t206;
t68 = t150 + t284;
t66 = mrSges(6,1) * t163 - t260;
t65 = -mrSges(6,2) * t163 - t271;
t51 = -mrSges(5,1) * t93 + mrSges(5,2) * t92;
t30 = t173 * t58 - t250;
t29 = -t169 * t58 - t247;
t28 = -mrSges(6,2) * t161 + mrSges(6,3) * t38;
t27 = mrSges(6,1) * t161 - mrSges(6,3) * t37;
t10 = -qJD(5) * t50 - t169 * t69 + t173 * t68;
t9 = qJD(5) * t49 + t169 * t68 + t173 * t69;
t7 = -mrSges(6,1) * t38 + mrSges(6,2) * t37;
t1 = [m(6) * (t10 * t25 + t102 * t89 + t118 * t48 + t26 * t9 + t3 * t50 + t4 * t49) + t143 * t51 + t118 * t7 + t102 * t47 + t10 * t66 + t9 * t65 + t49 * t27 + t50 * t28 + t179 + t149 * t251 + m(5) * (t204 * t252 + t193) + t107 * t206 + m(3) * (t171 * t99 + t175 * t98) * pkin(1) + Ifges(2,3) * qJDD(1) + m(4) * (t104 * t221 + t149 * t82 + t193) + t287 * t127 + t284 * t109 + (mrSges(3,1) * t267 - mrSges(3,2) * t269 + mrSges(4,3) * t143) * t162 + t276 * t134 + (-mrSges(2,1) * t176 + mrSges(2,2) * t172 - m(4) * t217 - m(5) * (t145 + t217) - m(6) * (t160 + t187) - m(3) * t160 + t278) * g(2) + (mrSges(2,1) * t172 + mrSges(2,2) * t176 - m(6) * (t188 + t285) - m(4) * t188 + m(3) * t268 - m(5) * (t136 - t268) + t279) * g(1) + (-mrSges(3,2) * t220 - t290 * t221) * t164; mrSges(3,2) * t208 - pkin(2) * t251 + qJ(3) * t51 + qJD(3) * t91 + t128 * t47 + t144 * t7 + t61 * t27 + t62 * t28 + t179 + t289 * t66 + t288 * t65 + (-t47 - t91) * t222 - t109 * t214 + (t144 * t48 + t3 * t62 + t4 * t61 + (t128 - t222) * t89 + t288 * t26 + t289 * t25) * m(6) + (t192 - (t204 + t243) * t253) * m(5) + (-pkin(2) * t82 + t192 - (t104 * t171 + t243) * t253) * m(4) + t276 * t178 + (-t208 + t283) * mrSges(4,3) + (-m(6) * t187 - m(5) * (t145 + t231) - m(4) * t231 + t278) * g(2) + (-m(5) * t136 - m(6) * (t212 + t285) - m(4) * t212 + t279) * g(1) + (-t107 * t170 - t109 * t174 + t290 * t164) * t223; t251 + t190 * t28 + t282 * t27 + t55 * t65 + t56 * t66 + t191 * qJD(4) + m(4) * t82 - m(6) * t293 + (t108 * t273 + t216 - t287) * t164 + t280 * (m(6) - t273) + t277; -t191 * t90 + (-t112 + (mrSges(6,2) * t155 - t200) * t156) * g(1) + (-t113 + (mrSges(6,1) * t157 + t200) * t158) * g(2) + Ifges(5,5) * t92 + Ifges(5,6) * t93 - t29 * t66 - t30 * t65 + t39 * mrSges(5,1) - t40 * mrSges(5,2) + t26 * t260 + t180 - m(6) * (t25 * t29 + t26 * t30) + (t85 * t272 + t170 * t86 / 0.2e1 + t216 * t174 * pkin(4) - t294) * t164 + ((g(3) * t170 + t280 * t174 - t25 * t225) * m(6) + (m(6) * t3 - qJD(5) * t66 + t28) * t169 + (t27 + qJD(5) * t65 + (t4 + t286) * m(6)) * t173) * pkin(4) + (t199 - t198) * g(3) + Ifges(5,3) * qJDD(4); (t66 + t260) * t26 - g(2) * (-mrSges(6,1) * t238 + t113) - g(1) * (-mrSges(6,2) * t242 + t112) - g(3) * t198 - t25 * t65 + t180;];
tau = t1;
