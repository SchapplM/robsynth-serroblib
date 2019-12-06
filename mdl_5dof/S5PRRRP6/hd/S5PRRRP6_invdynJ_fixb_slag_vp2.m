% Calculate vector of inverse dynamics joint torques for
% S5PRRRP6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:47
% EndTime: 2019-12-05 16:51:05
% DurationCPUTime: 7.53s
% Computational Cost: add. (2547->397), mult. (5521->522), div. (0->0), fcn. (3571->10), ass. (0->183)
t289 = mrSges(5,1) + mrSges(6,1);
t288 = mrSges(5,2) - mrSges(6,3);
t279 = Ifges(5,1) + Ifges(6,1);
t277 = Ifges(6,4) + Ifges(5,5);
t290 = m(6) + m(5);
t153 = sin(qJ(4));
t154 = sin(qJ(3));
t156 = cos(qJ(4));
t157 = cos(qJ(3));
t110 = t153 * t157 + t154 * t156;
t159 = -pkin(7) - pkin(6);
t194 = qJD(3) * t159;
t113 = t154 * t194;
t186 = t157 * t194;
t158 = cos(qJ(2));
t211 = qJD(1) * t158;
t122 = t159 * t154;
t123 = t159 * t157;
t68 = t122 * t153 - t123 * t156;
t273 = -qJD(4) * t68 + t110 * t211 - t113 * t153 + t156 * t186;
t217 = t153 * t154;
t109 = -t156 * t157 + t217;
t170 = t109 * t158;
t175 = t156 * t122 + t123 * t153;
t272 = qJD(1) * t170 + qJD(4) * t175 + t156 * t113 + t153 * t186;
t287 = -m(4) * pkin(6) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t150 = qJ(3) + qJ(4);
t144 = sin(t150);
t145 = cos(t150);
t184 = -mrSges(4,1) * t157 + mrSges(4,2) * t154;
t286 = m(4) * pkin(2) - t288 * t144 + t145 * t289 + mrSges(3,1) - t184;
t278 = -Ifges(5,4) + Ifges(6,5);
t285 = -Ifges(5,6) + Ifges(6,6);
t103 = t110 * qJD(2);
t147 = qJD(3) + qJD(4);
t204 = t157 * qJD(2);
t205 = t154 * qJD(2);
t102 = t153 * t205 - t156 * t204;
t231 = t102 * Ifges(6,5);
t99 = Ifges(5,4) * t102;
t284 = t279 * t103 + t277 * t147 + t231 - t99;
t58 = mrSges(5,1) * t102 + mrSges(5,2) * t103;
t283 = t184 * qJD(2) + t58;
t155 = sin(qJ(2));
t263 = t155 * t147;
t151 = sin(pkin(8));
t152 = cos(pkin(8));
t215 = t154 * t158;
t282 = t151 * t157 - t152 * t215;
t202 = qJD(2) * qJD(3);
t114 = qJDD(2) * t157 - t154 * t202;
t115 = qJDD(2) * t154 + t157 * t202;
t168 = t109 * qJD(4);
t40 = -qJD(2) * t168 + t114 * t153 + t115 * t156;
t260 = t40 / 0.2e1;
t280 = t114 / 0.2e1;
t146 = qJDD(3) + qJDD(4);
t253 = t146 / 0.2e1;
t21 = mrSges(5,1) * t146 - mrSges(5,3) * t40;
t22 = -t146 * mrSges(6,1) + t40 * mrSges(6,2);
t276 = t22 - t21;
t169 = t110 * qJD(4);
t41 = qJD(2) * t169 - t156 * t114 + t115 * t153;
t23 = -mrSges(5,2) * t146 - mrSges(5,3) * t41;
t24 = -mrSges(6,2) * t41 + mrSges(6,3) * t146;
t275 = t23 + t24;
t235 = mrSges(5,3) * t102;
t69 = -mrSges(5,2) * t147 - t235;
t237 = mrSges(6,2) * t102;
t72 = mrSges(6,3) * t147 - t237;
t244 = t69 + t72;
t234 = mrSges(5,3) * t103;
t236 = mrSges(6,2) * t103;
t274 = t289 * t147 - t234 - t236;
t120 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t205;
t121 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t204;
t271 = t157 * t120 + t154 * t121;
t270 = t120 * t154 - t121 * t157;
t269 = t157 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t114) - t154 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t115);
t203 = qJD(1) * qJD(2);
t136 = t158 * t203;
t117 = t155 * qJDD(1) + t136;
t107 = qJDD(2) * pkin(6) + t117;
t212 = qJD(1) * t155;
t127 = qJD(2) * pkin(6) + t212;
t209 = qJD(3) * t154;
t60 = t157 * t107 - t127 * t209;
t208 = qJD(3) * t157;
t61 = -t107 * t154 - t127 * t208;
t176 = -t154 * t61 + t157 * t60;
t218 = t152 * t158;
t89 = t144 * t218 - t151 * t145;
t90 = t144 * t151 + t145 * t218;
t268 = t288 * t90 + t289 * t89;
t219 = t151 * t158;
t87 = t144 * t219 + t145 * t152;
t88 = -t152 * t144 + t145 * t219;
t267 = t288 * t88 + t289 * t87;
t139 = pkin(3) * t157 + pkin(2);
t108 = -qJD(2) * t139 - t211;
t39 = pkin(4) * t102 - qJ(5) * t103 + t108;
t57 = mrSges(6,1) * t102 - mrSges(6,3) * t103;
t264 = -m(6) * t39 - t57;
t38 = qJDD(3) * pkin(3) - pkin(7) * t115 + t61;
t188 = pkin(7) * qJD(2) + t127;
t93 = t188 * t157;
t226 = t156 * t93;
t92 = t188 * t154;
t84 = qJD(3) * pkin(3) - t92;
t43 = t153 * t84 + t226;
t44 = pkin(7) * t114 + t60;
t6 = -qJD(4) * t43 - t153 * t44 + t156 * t38;
t160 = qJD(2) ^ 2;
t258 = -t102 / 0.2e1;
t257 = t102 / 0.2e1;
t255 = t103 / 0.2e1;
t251 = t147 / 0.2e1;
t250 = pkin(3) * t153;
t249 = pkin(3) * t154;
t248 = pkin(3) * t156;
t245 = g(3) * t155;
t238 = mrSges(5,2) * t145;
t233 = Ifges(4,4) * t154;
t232 = Ifges(4,4) * t157;
t230 = t103 * Ifges(5,4);
t229 = t153 * t93;
t213 = t157 * t158;
t210 = qJD(2) * t155;
t207 = qJD(4) * t153;
t206 = qJD(4) * t156;
t201 = pkin(3) * t205;
t200 = pkin(3) * t209;
t198 = pkin(3) * t206;
t135 = t155 * t203;
t192 = (-m(6) * qJ(5) - mrSges(6,3)) * t145 * t155;
t191 = -t87 * pkin(4) + qJ(5) * t88;
t190 = -t89 * pkin(4) + qJ(5) * t90;
t116 = qJDD(1) * t158 - t135;
t183 = mrSges(4,1) * t154 + mrSges(4,2) * t157;
t180 = t157 * Ifges(4,2) + t233;
t179 = Ifges(4,5) * t157 - Ifges(4,6) * t154;
t54 = pkin(4) * t103 + qJ(5) * t102;
t178 = pkin(4) * t145 + qJ(5) * t144;
t42 = t156 * t84 - t229;
t174 = t282 * pkin(3);
t173 = -t151 * t215 - t152 * t157;
t5 = t153 * t38 + t156 * t44 + t84 * t206 - t207 * t93;
t128 = -qJD(2) * pkin(2) - t211;
t172 = t128 * t183;
t171 = t154 * (Ifges(4,1) * t157 - t233);
t106 = -qJDD(2) * pkin(2) - t116;
t166 = t173 * pkin(3);
t165 = t128 * t155 + (t154 ^ 2 + t157 ^ 2) * t158 * t127;
t64 = -pkin(3) * t114 + t106;
t2 = qJ(5) * t146 + qJD(5) * t147 + t5;
t25 = -pkin(4) * t147 + qJD(5) - t42;
t3 = -pkin(4) * t146 + qJDD(5) - t6;
t30 = qJ(5) * t147 + t43;
t98 = Ifges(6,5) * t103;
t50 = t147 * Ifges(6,6) + t102 * Ifges(6,3) + t98;
t51 = -t102 * Ifges(5,2) + t147 * Ifges(5,6) + t230;
t161 = t237 * t25 - t39 * (mrSges(6,1) * t103 + mrSges(6,3) * t102) - t108 * (mrSges(5,1) * t103 - mrSges(5,2) * t102) - t42 * t235 + t30 * t236 + t2 * mrSges(6,3) - t3 * mrSges(6,1) - t5 * mrSges(5,2) + t6 * mrSges(5,1) + t51 * t255 + (Ifges(6,3) * t103 - t231) * t258 + t285 * t41 + t277 * t40 - (-t102 * t277 + t103 * t285) * t147 / 0.2e1 + (Ifges(6,2) + Ifges(5,3)) * t146 + (-Ifges(5,2) * t103 + t284 - t99) * t257 - (-t279 * t102 - t230 + t50 + t98) * t103 / 0.2e1;
t140 = Ifges(4,4) * t204;
t138 = -pkin(4) - t248;
t134 = qJ(5) + t250;
t131 = qJD(5) + t198;
t101 = Ifges(4,1) * t205 + Ifges(4,5) * qJD(3) + t140;
t100 = Ifges(4,6) * qJD(3) + qJD(2) * t180;
t95 = t109 * t155;
t94 = t110 * t155;
t65 = -mrSges(4,1) * t114 + mrSges(4,2) * t115;
t63 = qJD(3) * t110 + t169;
t62 = -qJD(3) * t109 - t168;
t59 = pkin(4) * t109 - qJ(5) * t110 - t139;
t46 = -t156 * t92 - t229;
t45 = -t153 * t92 + t226;
t14 = t153 * t158 * t204 + (t157 * t263 + t158 * t205) * t156 - t217 * t263;
t13 = -qJD(2) * t170 - t110 * t263;
t10 = pkin(4) * t63 - qJ(5) * t62 - qJD(5) * t110 + t200;
t9 = mrSges(5,1) * t41 + mrSges(5,2) * t40;
t8 = mrSges(6,1) * t41 - mrSges(6,3) * t40;
t7 = pkin(4) * t41 - qJ(5) * t40 - qJD(5) * t103 + t64;
t1 = [m(2) * qJDD(1) - t275 * t95 + t276 * t94 - t274 * t14 + t244 * t13 + (-m(2) - m(3) - m(4) - t290) * g(3) + (qJDD(2) * mrSges(3,1) - t160 * mrSges(3,2) - qJD(2) * t270 - t65 - t8 - t9) * t158 + (-t160 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t271 * qJD(3) + (t57 + t283) * qJD(2) + t269) * t155 + m(4) * (qJD(2) * t165 - t106 * t158 + t155 * t176) + m(5) * (t108 * t210 + t13 * t43 - t14 * t42 - t158 * t64 - t5 * t95 - t6 * t94) + m(3) * (t116 * t158 + t117 * t155) + m(6) * (t13 * t30 + t14 * t25 - t158 * t7 - t2 * t95 + t210 * t39 + t3 * t94); (t157 * (-Ifges(4,2) * t154 + t232) + t171) * t202 / 0.2e1 + t115 * t232 / 0.2e1 + (-t117 + t136) * mrSges(3,2) + t180 * t280 + (-pkin(2) * t106 - qJD(1) * t165) * m(4) + qJDD(3) * (Ifges(4,5) * t154 + Ifges(4,6) * t157) - t139 * t9 + (Ifges(4,1) * t115 + Ifges(4,4) * t280) * t154 + (mrSges(5,2) * t64 + mrSges(6,2) * t3 - mrSges(5,3) * t6 - mrSges(6,3) * t7 + (Ifges(6,5) / 0.2e1 - Ifges(5,4) / 0.2e1) * t41 + t279 * t260 + t277 * t253) * t110 + (t146 * t277 + t278 * t41 + t279 * t40) * t110 / 0.2e1 + t275 * t68 + t176 * mrSges(4,3) + (m(4) * t176 - t120 * t208 - t121 * t209 + t269) * pkin(6) + t270 * t211 - pkin(2) * t65 + t10 * t57 + t59 * t8 - t276 * t175 + (t108 * t200 - t139 * t64 + t175 * t6 + t272 * t43 + t273 * t42 + t5 * t68) * m(5) + (t10 * t39 - t175 * t3 + t2 * t68 - t25 * t273 + t272 * t30 + t59 * t7) * m(6) + (-t290 * (t158 * t139 - t155 * t159) + (-m(6) * t178 - t286) * t158 + t287 * t155) * g(3) + (t108 * mrSges(5,1) + t39 * mrSges(6,1) + t50 / 0.2e1 - t51 / 0.2e1 - t30 * mrSges(6,2) - t43 * mrSges(5,3) + Ifges(6,3) * t257 - Ifges(5,2) * t258 + t278 * t255 + t285 * t251) * t63 + (mrSges(5,1) * t64 + mrSges(6,1) * t7 - mrSges(6,2) * t2 - mrSges(5,3) * t5 + (Ifges(6,3) + Ifges(5,2)) * t41 + 0.2e1 * t278 * t260 + 0.2e1 * t285 * t253) * t109 + (-m(5) * t108 + t264 - t283) * t212 + t157 * (Ifges(4,4) * t115 + Ifges(4,2) * t114) / 0.2e1 + (t172 + t179 * qJD(3) / 0.2e1) * qJD(3) + t272 * t244 + t273 * t274 + (t116 + t135) * mrSges(3,1) + (g(1) * t152 + g(2) * t151) * ((t290 * t159 + t287) * t158 + (m(5) * t139 - m(6) * (-t139 - t178) + t286) * t155) + (mrSges(5,2) * t108 + t25 * mrSges(6,2) - t42 * mrSges(5,3) - mrSges(6,3) * t39 + Ifges(5,4) * t258 + Ifges(6,5) * t257 + t251 * t277 + t255 * t279 + t284 / 0.2e1) * t62 + t58 * t200 + t106 * t184 + t101 * t208 / 0.2e1 - t100 * t209 / 0.2e1 + Ifges(3,3) * qJDD(2); -(-Ifges(4,2) * t205 + t101 + t140) * t204 / 0.2e1 + (-t282 * mrSges(4,1) - (-t151 * t154 - t152 * t213) * mrSges(4,2) - m(6) * (t174 + t190) - m(5) * t174 + t268) * g(1) + t43 * t234 + ((t153 * t5 + t156 * t6 + (-t153 * t42 + t156 * t43) * qJD(4)) * pkin(3) - t108 * t201 + t42 * t45 - t43 * t46 + t249 * t245) * m(5) + m(6) * (t131 * t30 + t134 * t2 + t138 * t3) + t131 * t72 + t134 * t24 + t138 * t22 + Ifges(4,6) * t114 + Ifges(4,5) * t115 + (-m(6) * t30 - t244) * t46 + (-t173 * mrSges(4,1) - (-t151 * t213 + t152 * t154) * mrSges(4,2) - m(6) * (t166 + t191) - m(5) * t166 + t267) * g(2) + t271 * t127 + t264 * (t54 + t201) - t60 * mrSges(4,2) + t61 * mrSges(4,1) + (mrSges(5,1) * t144 + t183 + t238) * t245 + (pkin(3) * t207 - t45) * (m(6) * t25 - t274) + t161 + t69 * t198 - qJD(2) * t172 + t21 * t248 + t23 * t250 + Ifges(4,3) * qJDD(3) - t58 * t201 - t179 * t202 / 0.2e1 + t100 * t205 / 0.2e1 - g(3) * ((m(6) * (-pkin(4) * t144 - t249) - t144 * mrSges(6,1)) * t155 - t192) - t160 * t171 / 0.2e1; ((t238 + (m(6) * pkin(4) + t289) * t144) * t155 + t192) * g(3) + qJD(5) * t72 - t54 * t57 + (-m(6) * t191 + t267) * g(2) + (-m(6) * t190 + t268) * g(1) - pkin(4) * t22 + qJ(5) * t24 + t161 - m(6) * (t25 * t43 + t30 * t42 + t39 * t54) + m(6) * (-pkin(4) * t3 + qJ(5) * t2 + qJD(5) * t30) - t244 * t42 + (t234 + t274) * t43; t103 * t57 - t147 * t72 + (-g(1) * t89 - g(2) * t87 + t39 * t103 - t144 * t245 - t30 * t147 + t3) * m(6) + t22;];
tau = t1;
