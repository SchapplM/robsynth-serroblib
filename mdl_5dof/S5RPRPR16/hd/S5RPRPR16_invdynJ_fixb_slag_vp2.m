% Calculate vector of inverse dynamics joint torques for
% S5RPRPR16
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR16_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR16_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR16_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:45
% EndTime: 2019-12-31 18:38:59
% DurationCPUTime: 9.06s
% Computational Cost: add. (2014->426), mult. (3788->567), div. (0->0), fcn. (1848->6), ass. (0->202)
t110 = sin(qJ(1));
t241 = g(1) * t110;
t275 = Ifges(5,4) - Ifges(4,5);
t274 = Ifges(5,5) - Ifges(4,6);
t108 = sin(qJ(5));
t111 = cos(qJ(5));
t112 = cos(qJ(3));
t109 = sin(qJ(3));
t195 = qJD(3) * t109;
t115 = -pkin(1) - pkin(6);
t83 = qJDD(1) * t115 + qJDD(2);
t84 = qJD(1) * t115 + qJD(2);
t30 = t112 * t83 - t84 * t195;
t155 = qJDD(4) - t30;
t245 = pkin(3) + pkin(7);
t184 = qJD(1) * qJD(3);
t72 = -t112 * qJDD(1) + t109 * t184;
t13 = -pkin(4) * t72 - qJDD(3) * t245 + t155;
t159 = -qJD(4) * t112 + qJD(2);
t183 = qJDD(1) * qJ(2);
t118 = qJ(4) * t72 + qJD(1) * t159 + t183;
t73 = qJDD(1) * t109 + t112 * t184;
t8 = t245 * t73 + t118;
t41 = (pkin(4) * qJD(1) - t84) * t112;
t264 = t41 + qJD(4);
t32 = -qJD(3) * t245 + t264;
t101 = t112 * qJ(4);
t242 = pkin(7) * t109;
t136 = -t101 + t242;
t198 = qJD(1) * t109;
t211 = pkin(3) * t198 + qJD(1) * qJ(2);
t36 = qJD(1) * t136 + t211;
t9 = -t108 * t36 + t111 * t32;
t1 = qJD(5) * t9 + t108 * t13 + t111 * t8;
t279 = t1 * mrSges(6,2);
t10 = t108 * t32 + t111 * t36;
t2 = -qJD(5) * t10 - t108 * t8 + t111 * t13;
t278 = t2 * mrSges(6,1);
t152 = mrSges(4,1) * t109 + mrSges(4,2) * t112;
t277 = -t109 * mrSges(6,3) - t152;
t196 = qJD(3) * t108;
t64 = t111 * t198 - t196;
t23 = qJD(5) * t64 + qJDD(3) * t111 + t108 * t73;
t63 = qJDD(5) - t72;
t11 = mrSges(6,1) * t63 - mrSges(6,3) * t23;
t194 = qJD(3) * t111;
t65 = t108 * t198 + t194;
t24 = -qJD(5) * t65 - qJDD(3) * t108 + t111 * t73;
t12 = -mrSges(6,2) * t63 + mrSges(6,3) * t24;
t186 = t112 * qJD(1);
t87 = qJD(5) + t186;
t34 = -mrSges(6,2) * t87 + mrSges(6,3) * t64;
t35 = mrSges(6,1) * t87 - mrSges(6,3) * t65;
t254 = -t108 * t35 + t111 * t34;
t276 = -qJD(5) * t254 - t108 * t12 - t111 * t11;
t260 = m(5) + m(6);
t232 = Ifges(4,4) * t109;
t147 = t112 * Ifges(4,1) - t232;
t273 = Ifges(4,5) * qJD(3) + t65 * Ifges(6,5) + t64 * Ifges(6,6) + t87 * Ifges(6,3) + qJD(1) * t147;
t78 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t198;
t178 = mrSges(5,1) * t198;
t80 = -qJD(3) * mrSges(5,3) + t178;
t272 = t80 - t78;
t177 = mrSges(5,1) * t186;
t271 = mrSges(4,3) * t186 + t177 + (-mrSges(4,1) + mrSges(5,2)) * qJD(3);
t28 = -qJDD(3) * pkin(3) + t155;
t197 = qJD(3) * qJ(4);
t74 = t109 * t84;
t52 = -t74 - t197;
t270 = -qJD(3) * t52 - t28;
t113 = cos(qJ(1));
t240 = g(2) * t113;
t157 = -t240 + t241;
t269 = t260 * t157;
t227 = Ifges(5,6) * t112;
t268 = t109 * (-Ifges(4,2) * t112 - t232) + t112 * (Ifges(5,2) * t109 + t227);
t267 = t275 * t109 + t274 * t112;
t185 = qJD(1) * qJD(2);
t85 = t183 + t185;
t193 = qJD(3) * t112;
t31 = t109 * t83 + t84 * t193;
t266 = t109 * t31 + t112 * t30;
t27 = -qJDD(3) * qJ(4) - qJD(3) * qJD(4) - t31;
t265 = -t109 * t27 - t112 * t28;
t137 = t109 * Ifges(5,3) - t227;
t243 = Ifges(6,4) * t65;
t19 = Ifges(6,2) * t64 + Ifges(6,6) * t87 + t243;
t60 = Ifges(6,4) * t64;
t20 = Ifges(6,1) * t65 + Ifges(6,5) * t87 + t60;
t263 = Ifges(5,5) * qJD(3) + qJD(1) * t137 + t108 * t20 + t111 * t19;
t262 = mrSges(2,1) + mrSges(5,1) + mrSges(4,3) - mrSges(3,2);
t223 = t109 * mrSges(5,2);
t261 = -m(6) * t136 + mrSges(2,2) - mrSges(3,3) + t223 - (-m(5) * qJ(4) - mrSges(5,3)) * t112 + t277;
t250 = t23 / 0.2e1;
t249 = t24 / 0.2e1;
t248 = t63 / 0.2e1;
t44 = mrSges(5,1) * t73 - qJDD(3) * mrSges(5,3);
t7 = -mrSges(6,1) * t24 + mrSges(6,2) * t23;
t259 = -t44 + t7;
t236 = pkin(4) - t115;
t29 = -mrSges(6,1) * t64 + mrSges(6,2) * t65;
t257 = -t80 + t29;
t256 = pkin(3) * t193 + qJ(4) * t195;
t189 = qJD(5) * t111;
t191 = qJD(5) * t108;
t255 = -t10 * t189 + t9 * t191;
t45 = -t72 * mrSges(5,1) + qJDD(3) * mrSges(5,2);
t252 = t276 - t45;
t116 = qJD(1) ^ 2;
t251 = Ifges(6,1) * t250 + Ifges(6,4) * t249 + Ifges(6,5) * t248;
t246 = t65 / 0.2e1;
t103 = t109 * pkin(3);
t239 = t111 * t2;
t204 = t110 * t112;
t207 = t109 * t110;
t234 = pkin(3) * t204 + qJ(4) * t207;
t68 = pkin(3) * t186 + qJ(4) * t198;
t233 = mrSges(6,3) * t108;
t231 = Ifges(4,4) * t112;
t230 = Ifges(6,4) * t108;
t229 = Ifges(6,4) * t111;
t228 = Ifges(5,6) * t109;
t214 = t112 * t84;
t40 = -pkin(4) * t198 + t74;
t37 = t40 + t197;
t210 = qJD(3) * t37;
t208 = t108 * t109;
t206 = t109 * t111;
t203 = t111 * t112;
t202 = t112 * t113;
t200 = t103 - t101;
t199 = t113 * pkin(1) + t110 * qJ(2);
t190 = qJD(5) * t109;
t187 = qJDD(1) * mrSges(3,2);
t180 = Ifges(6,5) * t23 + Ifges(6,6) * t24 + Ifges(6,3) * t63;
t175 = t113 * pkin(6) + t199;
t167 = -t191 / 0.2e1;
t164 = (t185 + t85) * qJ(2);
t163 = qJD(3) * t236;
t162 = -t184 / 0.2e1;
t158 = -t200 - t242;
t156 = t1 * t108 + t239;
t154 = t10 * t111 - t9 * t108;
t153 = mrSges(4,1) * t112 - mrSges(4,2) * t109;
t151 = mrSges(6,1) * t111 - mrSges(6,2) * t108;
t150 = mrSges(6,1) * t108 + mrSges(6,2) * t111;
t149 = -t112 * mrSges(5,2) + t109 * mrSges(5,3);
t148 = -t112 * mrSges(5,3) - t223;
t146 = Ifges(6,1) * t111 - t230;
t145 = Ifges(6,1) * t108 + t229;
t144 = -t109 * Ifges(4,2) + t231;
t142 = -Ifges(6,2) * t108 + t229;
t141 = Ifges(6,2) * t111 + t230;
t139 = Ifges(6,5) * t111 - Ifges(6,6) * t108;
t138 = Ifges(6,5) * t108 + Ifges(6,6) * t111;
t58 = qJ(2) - t158;
t77 = t236 * t112;
t26 = t108 * t77 + t111 * t58;
t25 = -t108 * t58 + t111 * t77;
t46 = -qJD(3) * pkin(3) + qJD(4) - t214;
t133 = t46 * t109 - t52 * t112;
t69 = t148 * qJD(1);
t132 = -t254 - t69;
t47 = -qJ(4) * t186 + t211;
t131 = t47 * t149;
t130 = qJ(2) * t153;
t128 = t112 * (-Ifges(4,1) * t109 - t231);
t126 = -t116 * qJ(2) - t157;
t125 = -t108 * t190 + t111 * t193;
t124 = t108 * t193 + t109 * t189;
t121 = -Ifges(6,5) * t109 + t112 * t145;
t120 = -Ifges(6,6) * t109 + t112 * t141;
t119 = -Ifges(6,3) * t109 + t112 * t138;
t117 = qJD(5) * t154 + t156;
t102 = t113 * qJ(2);
t97 = -qJDD(1) * pkin(1) + qJDD(2);
t93 = Ifges(5,6) * t198;
t76 = t236 * t109;
t75 = qJ(2) + t200;
t70 = t152 * qJD(1);
t62 = t112 * t163;
t61 = t109 * t163;
t59 = t151 * t109;
t56 = -t108 * t204 + t111 * t113;
t55 = -t108 * t113 - t110 * t203;
t54 = -t108 * t202 - t110 * t111;
t53 = t108 * t110 - t111 * t202;
t51 = Ifges(5,4) * qJD(3) - Ifges(5,2) * t186 + t93;
t48 = Ifges(4,6) * qJD(3) + qJD(1) * t144;
t43 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t73;
t42 = qJDD(3) * mrSges(4,1) + mrSges(4,3) * t72;
t39 = pkin(7) * t186 + t68;
t38 = t159 + t256;
t33 = qJD(2) + (qJD(3) * pkin(7) - qJD(4)) * t112 + t256;
t17 = pkin(3) * t73 + t118;
t16 = t108 * t40 + t111 * t39;
t15 = -t108 * t39 + t111 * t40;
t14 = -pkin(4) * t73 - t27;
t6 = -qJD(5) * t26 - t108 * t33 - t111 * t61;
t5 = qJD(5) * t25 - t108 * t61 + t111 * t33;
t3 = t23 * Ifges(6,4) + t24 * Ifges(6,2) + t63 * Ifges(6,6);
t4 = [t112 * t278 + (-m(3) * t199 - m(4) * t175 - t56 * mrSges(6,1) - t55 * mrSges(6,2) - t260 * (pkin(3) * t207 + t175) + (-m(6) * pkin(4) - t262) * t113 + t261 * t110) * g(2) + (t274 * t109 - t275 * t112) * qJDD(3) / 0.2e1 + (t263 / 0.2e1 + t52 * mrSges(5,1) - t48 / 0.2e1) * t193 + (-t54 * mrSges(6,1) - t53 * mrSges(6,2) - t260 * (t113 * t103 + t102) + (-m(3) - m(4)) * t102 + (m(3) * pkin(1) + t236 * m(6) + t262) * t110 + t261 * t113) * g(1) - t112 * t279 - t265 * mrSges(5,1) + t73 * t137 / 0.2e1 - t73 * t144 / 0.2e1 - t72 * t147 / 0.2e1 + t17 * t148 - t266 * mrSges(4,3) + t267 * qJD(3) ^ 2 / 0.2e1 + t268 * t162 + t37 * (-mrSges(6,1) * t125 + mrSges(6,2) * t124) + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + t109 * t19 * t167 - t112 * (Ifges(5,4) * qJDD(3) + Ifges(5,6) * t73) / 0.2e1 + t72 * t228 / 0.2e1 + (0.2e1 * mrSges(3,3) + t152) * t85 + m(5) * (t17 * t75 + t38 * t47) + t109 * (Ifges(5,5) * qJDD(3) + Ifges(5,6) * t72 + Ifges(5,3) * t73) / 0.2e1 - t109 * (-Ifges(4,4) * t72 - Ifges(4,2) * t73 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t97 * mrSges(3,2) + t75 * (-mrSges(5,2) * t73 + mrSges(5,3) * t72) - t76 * t7 + t38 * t69 + qJD(2) * t70 + qJ(2) * (mrSges(4,1) * t73 - mrSges(4,2) * t72) - t14 * t59 - t62 * t29 + t5 * t34 + t6 * t35 + t25 * t11 + t26 * t12 + (t1 * t206 + t10 * t125 - t124 * t9 - t2 * t208) * mrSges(6,3) + (qJD(3) * t121 + t146 * t190) * t246 + (Ifges(6,3) * t112 + t109 * t138) * t248 + (Ifges(6,6) * t112 + t109 * t141) * t249 + (Ifges(6,5) * t112 + t109 * t145) * t250 + t208 * t251 + (-t273 / 0.2e1 - t46 * mrSges(5,1) + t51 / 0.2e1 - t9 * mrSges(6,1) + t10 * mrSges(6,2)) * t195 + ((t43 - t44) * t109 + (t42 - t45) * t112 + (-m(4) - m(5)) * t241 + m(5) * (qJD(3) * t133 + t265) + m(4) * t266 + (t109 * t271 - t112 * t272) * qJD(3)) * t115 - t72 * t112 * Ifges(5,2) + qJD(3) * t131 + t130 * t184 + m(3) * (-pkin(1) * t97 + t164) + (t109 * (Ifges(5,3) * t112 + t228) + t128) * t184 / 0.2e1 + (qJD(5) * t20 + t3) * t206 / 0.2e1 + (-Ifges(4,1) * t72 - Ifges(4,4) * t73 + Ifges(4,5) * qJDD(3) + t180) * t112 / 0.2e1 + m(6) * (t1 * t26 + t10 * t5 - t14 * t76 + t2 * t25 - t37 * t62 + t6 * t9) - pkin(1) * t187 + t64 * (qJD(3) * t120 + t142 * t190) / 0.2e1 + t87 * (qJD(3) * t119 + t139 * t190) / 0.2e1 + m(4) * t164; t187 - t116 * mrSges(3,3) + t126 * m(4) + (t126 + t97) * m(3) + (t43 + (t108 * t34 + t111 * t35 + t271) * qJD(3) + m(6) * (t10 * t196 + t194 * t9 + t14) + m(5) * (qJD(3) * t46 - t27) + m(4) * t31 + t259) * t109 + (t42 + (t78 + t257) * qJD(3) + m(6) * (-t156 + t210 + t255) + m(5) * t270 + m(4) * t30 + t252) * t112 + (-m(5) * t47 - m(6) * t154 + t132 - t70) * qJD(1) - t269; -t263 * t186 / 0.2e1 + (t268 / 0.2e1 - t130 - t128 / 0.2e1) * t116 + t274 * t73 + t275 * t72 + t272 * t214 + t273 * t198 / 0.2e1 + (m(5) * t200 - m(6) * t158 - t112 * t150 + t148 - t277) * g(3) + (-m(6) * t117 + t276) * t245 + t14 * t150 + (-m(6) * t234 + (-t149 - (m(6) * pkin(7) + mrSges(6,3)) * t112 - t150 * t109) * t110) * g(1) + (qJ(4) * t14 - t10 * t16 - t15 * t9 + t264 * t37) * m(6) + t267 * t162 - t271 * t74 + t46 * t178 - (Ifges(5,3) * t186 + t51 + t93) * t198 / 0.2e1 - (t138 * t87 + t141 * t64 + t145 * t65) * qJD(5) / 0.2e1 - (t119 * t87 + t120 * t64 + t121 * t65) * qJD(1) / 0.2e1 - t108 * t3 / 0.2e1 - t68 * t69 + t41 * t29 - pkin(3) * t45 + t30 * mrSges(4,1) - t31 * mrSges(4,2) - t16 * t34 - t15 * t35 - t27 * mrSges(5,3) + t28 * mrSges(5,2) + (t153 - m(5) * (-pkin(3) * t112 - qJ(4) * t109) + t149 - (-m(6) * t245 - mrSges(6,3)) * t112 - (-m(6) * qJ(4) - t150) * t109) * t240 + (-t131 - t10 * (mrSges(6,2) * t109 + mrSges(6,3) * t203) - t9 * (-mrSges(6,1) * t109 - t112 * t233)) * qJD(1) + (Ifges(5,1) + Ifges(4,3)) * qJDD(3) + t111 * t251 + t139 * t248 + t142 * t249 + t146 * t250 + (-pkin(3) * t28 - g(1) * t234 - qJ(4) * t27 - qJD(4) * t52 - t133 * t84 - t47 * t68) * m(5) + (-t239 + t255) * mrSges(6,3) + t257 * qJD(4) + t259 * qJ(4) + t87 * t37 * t151 + t20 * t167 - t52 * t177 + t48 * t186 / 0.2e1 - t19 * t189 / 0.2e1 - t1 * t233 - t153 * t241; -t257 * qJD(3) - t260 * t109 * g(3) + (-qJD(1) * t132 + t269) * t112 - t252 + (t154 * t186 + t117 - t210) * m(6) + (t186 * t47 - t270) * m(5); -t279 + t278 - t37 * (mrSges(6,1) * t65 + mrSges(6,2) * t64) - t65 * (Ifges(6,1) * t64 - t243) / 0.2e1 + t19 * t246 - t87 * (Ifges(6,5) * t64 - Ifges(6,6) * t65) / 0.2e1 - t9 * t34 + t10 * t35 - g(1) * (mrSges(6,1) * t55 - mrSges(6,2) * t56) - g(2) * (-mrSges(6,1) * t53 + mrSges(6,2) * t54) - g(3) * t59 + (t10 * t65 + t64 * t9) * mrSges(6,3) + t180 - (-Ifges(6,2) * t65 + t20 + t60) * t64 / 0.2e1;];
tau = t4;
