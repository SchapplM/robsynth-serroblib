% Calculate vector of inverse dynamics joint torques for
% S5RPRRP1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:14
% EndTime: 2019-12-05 17:59:32
% DurationCPUTime: 6.12s
% Computational Cost: add. (2570->342), mult. (4943->432), div. (0->0), fcn. (2873->8), ass. (0->160)
t264 = Ifges(5,4) + Ifges(6,4);
t265 = Ifges(5,1) + Ifges(6,1);
t263 = Ifges(5,5) + Ifges(6,5);
t262 = Ifges(5,2) + Ifges(6,2);
t261 = Ifges(5,6) + Ifges(6,6);
t149 = qJ(3) + qJ(4);
t139 = sin(t149);
t208 = mrSges(5,2) + mrSges(6,2);
t269 = t139 * t208;
t140 = cos(t149);
t276 = t140 * t208;
t209 = mrSges(5,1) + mrSges(6,1);
t268 = t140 * t209;
t275 = m(5) + m(6);
t150 = sin(qJ(3));
t152 = cos(qJ(4));
t153 = cos(qJ(3));
t229 = sin(qJ(4));
t101 = -t152 * t150 - t153 * t229;
t92 = t101 * qJD(1);
t274 = t264 * t92;
t178 = t229 * t150;
t187 = t153 * qJD(1);
t93 = -qJD(1) * t178 + t152 * t187;
t273 = t264 * t93;
t272 = t153 / 0.2e1;
t147 = qJD(3) + qJD(4);
t271 = t261 * t147 + t262 * t92 + t273;
t270 = t263 * t147 + t265 * t93 + t274;
t143 = t150 * pkin(3);
t171 = mrSges(4,1) * t150 + mrSges(4,2) * t153;
t267 = t171 + t275 * t143 + t276 + (m(6) * pkin(4) + t209) * t139;
t235 = -m(3) - m(4);
t251 = t235 - t275;
t145 = qJDD(3) + qJDD(4);
t185 = qJD(1) * qJD(3);
t106 = qJDD(1) * t153 - t150 * t185;
t107 = -qJDD(1) * t150 - t153 * t185;
t196 = t152 * t153;
t165 = t178 - t196;
t37 = qJD(1) * qJD(4) * t165 - t106 * t229 + t152 * t107;
t25 = -mrSges(6,2) * t145 + mrSges(6,3) * t37;
t26 = -mrSges(5,2) * t145 + mrSges(5,3) * t37;
t260 = t25 + t26;
t226 = mrSges(6,3) * t92;
t70 = -mrSges(6,2) * t147 + t226;
t227 = mrSges(5,3) * t92;
t71 = -mrSges(5,2) * t147 + t227;
t259 = t70 + t71;
t212 = t93 * mrSges(6,3);
t72 = mrSges(6,1) * t147 - t212;
t213 = t93 * mrSges(5,3);
t73 = mrSges(5,1) * t147 - t213;
t258 = t72 + t73;
t156 = -pkin(1) - pkin(6);
t122 = qJD(1) * t156 + qJD(2);
t188 = t150 * qJD(1);
t80 = -pkin(7) * t188 + t122 * t150;
t75 = t229 * t80;
t108 = t153 * t122;
t81 = -pkin(7) * t187 + t108;
t78 = qJD(3) * pkin(3) + t81;
t38 = t152 * t78 - t75;
t84 = t93 * qJ(5);
t14 = -t84 + t38;
t151 = sin(qJ(1));
t257 = t151 * t268;
t154 = cos(qJ(1));
t256 = t154 * t269;
t224 = pkin(3) * t153;
t254 = -m(5) * t224 - m(6) * (pkin(4) * t140 + t224);
t186 = qJD(1) * qJD(2);
t123 = qJDD(1) * qJ(2) + t186;
t121 = qJDD(1) * t156 + qJDD(2);
t193 = qJD(3) * t150;
t65 = t153 * t121 - t122 * t193;
t192 = qJD(3) * t153;
t66 = t150 * t121 + t122 * t192;
t167 = t150 * t66 + t153 * t65;
t253 = -g(1) * t151 + g(2) * t154;
t172 = mrSges(4,1) * t153 - mrSges(4,2) * t150;
t205 = Ifges(4,4) * t153;
t252 = (-Ifges(4,1) * t150 - t205) * t272 + qJ(2) * t172;
t113 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t188;
t197 = t150 * (qJD(3) * mrSges(4,1) - mrSges(4,3) * t187);
t250 = (t153 * t113 - t197) * qJD(3);
t76 = t152 * t80;
t166 = -t229 * t78 - t76;
t176 = qJD(4) * t229;
t190 = qJD(4) * t152;
t47 = qJDD(3) * pkin(3) - pkin(7) * t106 + t65;
t55 = pkin(7) * t107 + t66;
t5 = t152 * t55 - t176 * t80 + t78 * t190 + t229 * t47;
t6 = qJD(4) * t166 + t152 * t47 - t229 * t55;
t60 = -qJD(3) * t178 + t147 * t196 - t150 * t176;
t159 = t101 * qJD(4);
t61 = qJD(3) * t101 + t159;
t249 = t101 * t5 + t165 * t6 + t166 * t60 - t38 * t61;
t11 = t147 * pkin(4) + t14;
t200 = t92 * qJ(5);
t15 = -t166 + t200;
t36 = qJD(1) * t159 + t152 * t106 + t107 * t229;
t2 = t145 * pkin(4) - t36 * qJ(5) - t93 * qJD(5) + t6;
t3 = t37 * qJ(5) + t92 * qJD(5) + t5;
t248 = -t101 * t3 + t11 * t61 + t15 * t60 - t165 * t2;
t247 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) - mrSges(3,2);
t246 = m(4) * t167 + t153 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t106) + t150 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t107);
t245 = -mrSges(3,3) + mrSges(2,2) - t267;
t244 = qJD(1) ^ 2;
t236 = t93 / 0.2e1;
t155 = -pkin(7) - pkin(6);
t225 = pkin(3) * t152;
t207 = pkin(7) - t156;
t45 = t152 * t81 - t75;
t111 = t207 * t150;
t112 = t207 * t153;
t63 = -t152 * t111 - t229 * t112;
t206 = Ifges(4,4) * t150;
t128 = qJ(2) + t143;
t115 = pkin(3) * t188 + qJD(1) * qJ(2);
t189 = qJDD(1) * mrSges(3,2);
t124 = pkin(3) * t192 + qJD(2);
t181 = pkin(3) * t187;
t177 = -t37 * mrSges(6,1) + t36 * mrSges(6,2);
t175 = -t185 / 0.2e1;
t174 = (t123 + t186) * qJ(2);
t44 = -t229 * t81 - t76;
t62 = t111 * t229 - t152 * t112;
t170 = t153 * Ifges(4,1) - t206;
t169 = -t150 * Ifges(4,2) + t205;
t168 = -Ifges(4,5) * t150 - Ifges(4,6) * t153;
t163 = t150 * (-Ifges(4,2) * t153 - t206);
t98 = t207 * t193;
t99 = qJD(3) * t112;
t12 = t111 * t176 - t112 * t190 - t152 * t99 + t229 * t98;
t74 = -pkin(3) * t107 + t123;
t13 = -qJD(4) * t63 + t152 * t98 + t229 * t99;
t64 = -pkin(4) * t92 + qJD(5) + t115;
t158 = t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2) + t11 * t226 - t115 * (mrSges(5,1) * t93 + mrSges(5,2) * t92) - t64 * (mrSges(6,1) * t93 + mrSges(6,2) * t92) + t38 * t227 + t261 * t37 + t263 * t36 - (t265 * t92 - t273) * t93 / 0.2e1 + t271 * t236 - (-t261 * t93 + t263 * t92) * t147 / 0.2e1 + (Ifges(6,3) + Ifges(5,3)) * t145 - (-t262 * t93 + t270 + t274) * t92 / 0.2e1;
t146 = -qJ(5) + t155;
t137 = -qJDD(1) * pkin(1) + qJDD(2);
t131 = pkin(4) + t225;
t103 = t171 * qJD(1);
t91 = Ifges(4,5) * qJD(3) + qJD(1) * t170;
t90 = Ifges(4,6) * qJD(3) + qJD(1) * t169;
t79 = -pkin(4) * t101 + t128;
t69 = pkin(4) * t93 + t181;
t59 = -mrSges(5,1) * t92 + mrSges(5,2) * t93;
t58 = -mrSges(6,1) * t92 + mrSges(6,2) * t93;
t54 = pkin(4) * t60 + t124;
t43 = qJ(5) * t101 + t63;
t42 = qJ(5) * t165 + t62;
t24 = mrSges(5,1) * t145 - mrSges(5,3) * t36;
t23 = mrSges(6,1) * t145 - mrSges(6,3) * t36;
t18 = -t84 + t45;
t17 = t44 - t200;
t10 = -pkin(4) * t37 + qJDD(5) + t74;
t8 = -t61 * qJ(5) + qJD(5) * t165 + t13;
t7 = -t60 * qJ(5) + t101 * qJD(5) + t12;
t1 = [-t90 * t192 / 0.2e1 - t91 * t193 / 0.2e1 - pkin(1) * t189 + (-t264 * t60 + t265 * t61) * t236 + t79 * t177 + m(3) * (-pkin(1) * t137 + t174) + (-mrSges(5,2) * t74 - mrSges(6,2) * t10 - t145 * t263 - t264 * t37 - t265 * t36) * t165 + ((-m(5) * (-pkin(1) + t155) + m(3) * pkin(1) - m(4) * t156 - m(6) * (-pkin(1) + t146) + t247) * t151 + (t251 * qJ(2) + t245) * t154) * g(1) + qJDD(3) * (Ifges(4,5) * t153 - Ifges(4,6) * t150) + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + t128 * (-mrSges(5,1) * t37 + mrSges(5,2) * t36) + t137 * mrSges(3,2) + t115 * (mrSges(5,1) * t60 + mrSges(5,2) * t61) + t124 * t59 + qJD(2) * t103 + qJ(2) * (-mrSges(4,1) * t107 + mrSges(4,2) * t106) + t7 * t70 + t12 * t71 + t8 * t72 + t13 * t73 + t62 * t24 + t63 * t26 + t64 * (mrSges(6,1) * t60 + mrSges(6,2) * t61) + t54 * t58 + t42 * t23 + t43 * t25 + (-t261 * t60 + t263 * t61) * t147 / 0.2e1 + (-t262 * t60 + t264 * t61) * t92 / 0.2e1 - t150 * (Ifges(4,4) * t106 + Ifges(4,2) * t107) / 0.2e1 + m(4) * t174 + t252 * t185 - t167 * mrSges(4,3) + m(5) * (t115 * t124 - t12 * t166 + t128 * t74 + t13 * t38 + t5 * t63 + t6 * t62) + t270 * t61 / 0.2e1 - t271 * t60 / 0.2e1 + t156 * t250 + (-mrSges(5,1) * t74 - mrSges(6,1) * t10 + t145 * t261 + t262 * t37 + t264 * t36) * t101 + (t251 * (t154 * pkin(1) + t151 * qJ(2)) + (-m(4) * pkin(6) + m(5) * t155 + m(6) * t146 - t247) * t154 + t245 * t151) * g(2) - t248 * mrSges(6,3) + t249 * mrSges(5,3) + t246 * t156 + m(6) * (t10 * t79 + t11 * t8 + t15 * t7 + t2 * t42 + t3 * t43 + t54 * t64) + qJD(3) ^ 2 * t168 / 0.2e1 + t107 * t169 / 0.2e1 + t106 * t170 / 0.2e1 + (0.2e1 * mrSges(3,3) + t171) * t123 + t163 * t175 + (Ifges(4,1) * t106 + Ifges(4,4) * t107) * t272; t189 + t258 * t61 + t259 * t60 - (t23 + t24) * t165 - t260 * t101 + t250 + (qJ(2) * t235 - mrSges(3,3)) * t244 + (-m(5) * t115 - m(6) * t64 - t103 - t58 - t59) * qJD(1) + m(3) * t137 - m(5) * t249 + m(6) * t248 - t253 * t251 + t246; -t113 * t108 + t90 * t187 / 0.2e1 + t91 * t188 / 0.2e1 + t267 * g(3) - t59 * t181 + ((-t254 + t268) * t154 - t256) * g(2) + ((t254 + t269) * t151 - t257) * g(1) + t131 * t23 + Ifges(4,5) * t106 + Ifges(4,6) * t107 + t65 * mrSges(4,1) - t66 * mrSges(4,2) - t69 * t58 - t18 * t70 - t45 * t71 - t17 * t72 - t44 * t73 + (t163 / 0.2e1 - t252) * t244 + t253 * t172 + (-t11 * t17 + t2 * t131 - t15 * t18 - t64 * t69) * m(6) + (-t115 * t181 + t166 * t45 - t38 * t44) * m(5) + (t259 * t190 + t260 * t229 - t258 * t176 + (t229 * t5 + t152 * t6 + (-t152 * t166 - t229 * t38) * qJD(4)) * m(5) + (t229 * t3 + (-t11 * t229 + t15 * t152) * qJD(4)) * m(6)) * pkin(3) - t166 * t213 + t158 + Ifges(4,3) * qJDD(3) + t168 * t175 + t122 * t197 + t15 * t212 + t24 * t225; (t151 * t269 - t257) * g(1) + (t154 * t268 - t256) * g(2) + (t139 * t209 + t276) * g(3) - t14 * t70 - t38 * t71 + t15 * t72 - t166 * t73 + pkin(4) * t23 + (-(-t11 + t14) * t15 + (g(3) * t139 + t140 * t253 - t64 * t93 + t2) * pkin(4)) * m(6) + t158 + (-mrSges(5,3) * t166 + t15 * mrSges(6,3) - pkin(4) * t58) * t93; -t92 * t70 + t93 * t72 + (-g(1) * t154 - g(2) * t151 + t11 * t93 - t15 * t92 + t10) * m(6) + t177;];
tau = t1;
