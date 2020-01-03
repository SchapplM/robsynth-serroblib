% Calculate vector of inverse dynamics joint torques for
% S5RPRRR8
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:42
% EndTime: 2019-12-31 19:05:51
% DurationCPUTime: 4.04s
% Computational Cost: add. (3811->331), mult. (5583->443), div. (0->0), fcn. (3013->10), ass. (0->174)
t151 = -pkin(8) - pkin(7);
t276 = m(5) * pkin(7) - m(6) * t151 - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t149 = cos(qJ(4));
t229 = pkin(4) * t149;
t126 = pkin(3) + t229;
t146 = sin(qJ(4));
t107 = -mrSges(5,1) * t149 + mrSges(5,2) * t146;
t144 = qJ(4) + qJ(5);
t130 = sin(t144);
t131 = cos(t144);
t182 = -mrSges(6,1) * t131 + t130 * mrSges(6,2);
t270 = -t107 - t182;
t275 = m(5) * pkin(3) + m(6) * t126 + mrSges(4,1) + t270;
t274 = qJD(1) - qJD(3);
t147 = sin(qJ(3));
t150 = cos(qJ(3));
t231 = sin(qJ(1));
t232 = cos(qJ(1));
t89 = -t147 * t231 - t150 * t232;
t91 = t147 * t232 - t150 * t231;
t273 = t275 * t91 + t276 * t89;
t272 = -t275 * t89 + t276 * t91;
t271 = t150 * t274;
t235 = t146 / 0.2e1;
t152 = -pkin(1) - pkin(2);
t201 = qJ(2) * qJD(1);
t269 = -qJD(3) * t201 + qJDD(1) * t152 + qJDD(2);
t114 = qJD(1) * t152 + qJD(2);
t79 = t114 * t147 + t150 * t201;
t64 = -pkin(7) * t274 + t79;
t183 = (t146 ^ 2 + t149 ^ 2) * t64;
t209 = t114 * t150;
t78 = -t147 * t201 + t209;
t63 = pkin(3) * t274 - t78;
t268 = t147 * t63 + t150 * t183;
t248 = -g(1) * t89 - g(2) * t91;
t145 = sin(qJ(5));
t148 = cos(qJ(5));
t90 = t145 * t149 + t146 * t148;
t75 = t90 * t274;
t265 = t75 / 0.2e1;
t55 = t126 * t274 - t78;
t264 = m(6) * t55;
t140 = qJD(4) + qJD(5);
t263 = -t140 / 0.2e1;
t261 = mrSges(3,1) + mrSges(2,1);
t260 = -mrSges(3,3) + mrSges(2,2);
t111 = t151 * t146;
t112 = t151 * t149;
t61 = t111 * t145 - t112 * t148;
t190 = qJD(4) * t151;
t98 = t146 * t190;
t99 = t149 * t190;
t259 = -qJD(5) * t61 - t145 * t98 + t148 * t99 + t90 * t78;
t60 = t111 * t148 + t112 * t145;
t206 = t145 * t146;
t88 = -t148 * t149 + t206;
t258 = qJD(5) * t60 + t145 * t99 + t148 * t98 + t88 * t78;
t81 = t88 * t147;
t257 = t140 * t81 + t90 * t271;
t161 = t90 * qJD(5);
t256 = (-qJD(4) * t90 - t161) * t147 + t88 * t271;
t169 = mrSges(5,1) * t146 + mrSges(5,2) * t149;
t255 = t63 * t169;
t254 = (mrSges(4,1) - t107) * t274;
t186 = -pkin(8) * t274 + t64;
t54 = t186 * t149;
t215 = t145 * t54;
t53 = t186 * t146;
t51 = qJD(4) * pkin(4) - t53;
t14 = t148 * t51 - t215;
t253 = qJD(5) * t14;
t103 = -t147 * qJ(2) + t150 * t152;
t104 = t150 * qJ(2) + t147 * t152;
t195 = qJD(1) * qJD(2);
t116 = qJDD(1) * qJ(2) + t195;
t139 = -qJDD(1) + qJDD(3);
t199 = qJD(4) * t146;
t84 = t139 * t149 + t199 * t274;
t198 = qJD(4) * t149;
t85 = t139 * t146 - t198 * t274;
t249 = t149 * (-qJDD(4) * mrSges(5,2) + mrSges(5,3) * t84) - t146 * (qJDD(4) * mrSges(5,1) - mrSges(5,3) * t85);
t45 = qJD(3) * t209 + t150 * t116 + t269 * t147;
t34 = pkin(7) * t139 + t45;
t12 = t149 * t34 - t199 * t64;
t13 = -t146 * t34 - t198 * t64;
t166 = t12 * t149 - t13 * t146;
t247 = -m(6) - m(4) - m(5);
t208 = t274 * t146;
t101 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t208;
t207 = t274 * t149;
t102 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t207;
t158 = -mrSges(4,2) * t274 + t146 * t101 - t149 * t102;
t200 = qJD(3) * t147;
t46 = -t114 * t200 - t147 * t116 + t269 * t150;
t35 = -pkin(3) * t139 - t46;
t52 = -mrSges(5,1) * t84 + mrSges(5,2) * t85;
t245 = m(5) * t35 + t52;
t244 = m(5) * t63 + t254;
t243 = -m(5) * t183 + t158;
t241 = m(5) * t166 - t101 * t198 - t102 * t199 + t249;
t221 = Ifges(5,2) * t149;
t223 = Ifges(5,4) * t146;
t109 = t221 + t223;
t138 = qJDD(4) + qJDD(5);
t167 = Ifges(5,1) * t149 - t223;
t11 = pkin(8) * t84 + t12;
t9 = qJDD(4) * pkin(4) - pkin(8) * t85 + t13;
t2 = t11 * t148 + t145 * t9 + t253;
t22 = -pkin(4) * t84 + t35;
t222 = Ifges(5,4) * t149;
t212 = t148 * t54;
t15 = t145 * t51 + t212;
t227 = t15 * mrSges(6,3);
t228 = t14 * mrSges(6,3);
t74 = t88 * t274;
t27 = qJD(5) * t74 + t145 * t84 + t148 * t85;
t28 = -t145 * t85 + t148 * t84 + t161 * t274;
t3 = -qJD(5) * t15 - t11 * t145 + t148 * t9;
t226 = t75 * Ifges(6,4);
t37 = t74 * Ifges(6,2) + t140 * Ifges(6,6) - t226;
t67 = Ifges(6,4) * t74;
t38 = -t75 * Ifges(6,1) + t140 * Ifges(6,5) + t67;
t197 = qJD(5) * t148;
t56 = t140 * t206 - t148 * t198 - t149 * t197;
t57 = t140 * t90;
t72 = Ifges(5,6) * qJD(4) - t109 * t274;
t115 = Ifges(5,4) * t207;
t73 = -Ifges(5,1) * t208 + Ifges(5,5) * qJD(4) - t115;
t93 = qJD(4) * (Ifges(5,5) * t149 - Ifges(5,6) * t146);
t240 = t85 * (Ifges(5,1) * t146 + t222) / 0.2e1 + (Ifges(6,5) * t56 + Ifges(6,6) * t57) * t263 + (t93 - (-Ifges(5,2) * t146 + t222) * t207 - t167 * t208) * qJD(4) / 0.2e1 + (0.2e1 * Ifges(5,5) * t235 + Ifges(5,6) * t149) * qJDD(4) + t149 * (Ifges(5,4) * t85 + Ifges(5,2) * t84) / 0.2e1 + Ifges(4,3) * t139 + t35 * t107 + t84 * t109 / 0.2e1 + qJD(4) * t255 + (-t37 / 0.2e1 - t227) * t57 + (Ifges(5,1) * t85 + Ifges(5,4) * t84) * t235 + t55 * (mrSges(6,1) * t57 - mrSges(6,2) * t56) + t46 * mrSges(4,1) - t45 * mrSges(4,2) + (-t38 / 0.2e1 + t228) * t56 + t166 * mrSges(5,3) + t73 * t198 / 0.2e1 - t72 * t199 / 0.2e1 + (t22 * mrSges(6,2) - t3 * mrSges(6,3) + t27 * Ifges(6,1) + t28 * Ifges(6,4) + t138 * Ifges(6,5)) * t90 + (t22 * mrSges(6,1) - t2 * mrSges(6,3) - t27 * Ifges(6,4) - t28 * Ifges(6,2) - t138 * Ifges(6,6)) * t88;
t239 = -t74 / 0.2e1;
t238 = -t75 / 0.2e1;
t97 = -pkin(7) + t104;
t233 = pkin(8) - t97;
t230 = mrSges(6,3) * t75;
t218 = t139 * mrSges(4,1);
t217 = t139 * mrSges(4,2);
t202 = t232 * pkin(1) + t231 * qJ(2);
t196 = qJDD(1) * mrSges(3,1);
t193 = pkin(4) * t199;
t185 = qJD(4) * t233;
t48 = -mrSges(6,1) * t74 - mrSges(6,2) * t75;
t181 = -t48 - t254;
t96 = pkin(3) - t103;
t180 = -pkin(1) * t231 + t232 * qJ(2);
t175 = Ifges(6,1) * t56 + Ifges(6,4) * t57;
t173 = Ifges(6,4) * t56 + Ifges(6,2) * t57;
t168 = mrSges(6,1) * t130 + mrSges(6,2) * t131;
t70 = t233 * t146;
t71 = t233 * t149;
t41 = t145 * t71 + t148 * t70;
t42 = t145 * t70 - t148 * t71;
t165 = -t147 * t78 + t150 * t79;
t164 = t101 * t149 + t102 * t146;
t77 = qJD(2) * t147 + t104 * qJD(3);
t155 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + t74 * t228 + t37 * t238 - t55 * (-mrSges(6,1) * t75 + mrSges(6,2) * t74) + Ifges(6,3) * t138 + (Ifges(6,1) * t74 + t226) * t265 + Ifges(6,6) * t28 + Ifges(6,5) * t27 + (Ifges(6,5) * t74 + Ifges(6,6) * t75) * t263 + (Ifges(6,2) * t75 + t38 + t67) * t239;
t129 = -qJDD(1) * pkin(1) + qJDD(2);
t83 = t96 + t229;
t80 = t90 * t147;
t76 = qJD(2) * t150 + t103 * qJD(3);
t62 = t77 - t193;
t59 = mrSges(6,1) * t140 + t230;
t58 = -mrSges(6,2) * t140 + mrSges(6,3) * t74;
t44 = -t146 * t76 + t149 * t185;
t43 = t146 * t185 + t149 * t76;
t21 = -mrSges(6,2) * t138 + mrSges(6,3) * t28;
t20 = mrSges(6,1) * t138 - mrSges(6,3) * t27;
t17 = -t148 * t53 - t215;
t16 = t145 * t53 - t212;
t8 = -mrSges(6,1) * t28 + mrSges(6,2) * t27;
t5 = -qJD(5) * t42 - t145 * t43 + t148 * t44;
t4 = qJD(5) * t41 + t145 * t44 + t148 * t43;
t1 = [(-m(3) * t180 + t247 * (-pkin(2) * t231 + t180) + t260 * t232 + t261 * t231 - t273) * g(1) + (-m(3) * t202 + t247 * (t232 * pkin(2) + t202) - t261 * t232 + t260 * t231 - t272) * g(2) + m(4) * (t103 * t46 + t104 * t45) + 0.2e1 * t116 * mrSges(3,3) + t241 * t97 + t245 * t96 + (-m(4) * t78 + t244) * t77 + (m(4) * t79 - t243) * t76 + (Ifges(3,2) + Ifges(2,3)) * qJDD(1) + m(6) * (t14 * t5 + t15 * t4 + t2 * t42 + t22 * t83 + t3 * t41 + t55 * t62) - t129 * mrSges(3,1) + t83 * t8 + t4 * t58 + t5 * t59 + t62 * t48 + t41 * t20 + t42 * t21 - t240 + t74 * t173 / 0.2e1 + m(3) * (-pkin(1) * t129 + (t116 + t195) * qJ(2)) - t104 * t217 + pkin(1) * t196 + t103 * t218 + t175 * t238; -t196 - t80 * t20 - t81 * t21 + t257 * t59 + t256 * t58 + (-qJD(3) * t158 + t218 - t52 - t8) * t150 + (-qJD(3) * t181 - qJD(4) * t164 - t217 + t249) * t147 + m(5) * (t268 * qJD(3) + t166 * t147 - t150 * t35) + m(4) * (qJD(3) * t165 + t147 * t45 + t150 * t46) + m(3) * t129 + (t257 * t14 + t256 * t15 - t150 * t22 - t2 * t81 + t55 * t200 - t3 * t80) * m(6) + (t158 * t150 - m(5) * t268 - m(4) * t165 + (-m(3) * qJ(2) - mrSges(3,3)) * qJD(1) + (t181 - t264) * t147) * qJD(1) + (-g(1) * t231 + g(2) * t232) * (m(3) - t247); t259 * t59 + t258 * t58 + (-t126 * t22 + t259 * t14 + t258 * t15 + t193 * t55 + t2 * t61 + t3 * t60) * m(6) + t273 * g(1) + t272 * g(2) - t126 * t8 + t60 * t20 + t61 * t21 + t48 * t193 + t241 * pkin(7) + t240 - t245 * pkin(3) + t175 * t265 + t173 * t239 + t243 * t78 + (-t244 - t48 - t264) * t79; t164 * t64 + (t48 * t208 + (g(3) * t149 + t15 * t197 + (t274 * t55 + t248) * t146) * m(6) + (m(6) * t3 + qJD(5) * t58 + t20) * t148 + (t21 - t59 * qJD(5) + (t2 - t253) * m(6)) * t145) * pkin(4) - m(6) * (t14 * t16 + t15 * t17) + Ifges(5,5) * t85 + Ifges(5,6) * t84 - t17 * t58 - t16 * t59 - t12 * mrSges(5,2) + t13 * mrSges(5,1) - t75 * t227 + t155 + t270 * g(3) + Ifges(5,3) * qJDD(4) - (-t255 - t93 / 0.2e1 + t72 * t235 - (t221 * t235 - t146 * t167 / 0.2e1) * t274 - (-t115 + t73) * t149 / 0.2e1) * t274 - t248 * (-t168 - t169); (t59 - t230) * t15 - g(3) * t182 - t14 * t58 + t155 + t248 * t168;];
tau = t1;
