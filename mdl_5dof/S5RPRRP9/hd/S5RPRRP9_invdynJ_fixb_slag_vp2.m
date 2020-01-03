% Calculate vector of inverse dynamics joint torques for
% S5RPRRP9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:37
% EndTime: 2019-12-31 18:48:54
% DurationCPUTime: 9.91s
% Computational Cost: add. (4262->412), mult. (10546->512), div. (0->0), fcn. (7770->12), ass. (0->178)
t277 = mrSges(5,1) + mrSges(6,1);
t260 = Ifges(5,1) + Ifges(6,1);
t259 = Ifges(6,4) + Ifges(5,5);
t149 = sin(pkin(8));
t145 = t149 ^ 2;
t150 = cos(pkin(8));
t146 = t150 ^ 2;
t204 = t145 + t146;
t201 = qJD(1) * qJD(2);
t126 = qJDD(1) * qJ(2) + t201;
t148 = qJD(3) + qJD(4);
t152 = sin(qJ(4));
t153 = sin(qJ(3));
t155 = cos(qJ(3));
t207 = t150 * t155;
t173 = t149 * t153 - t207;
t169 = t173 * qJD(1);
t218 = pkin(6) + qJ(2);
t120 = t218 * t149;
t111 = qJD(1) * t120;
t121 = t218 * t150;
t112 = qJD(1) * t121;
t82 = -t111 * t153 + t112 * t155;
t56 = -pkin(7) * t169 + t82;
t211 = t152 * t56;
t236 = cos(qJ(4));
t110 = t149 * t155 + t150 * t153;
t101 = t110 * qJD(1);
t210 = t112 * t153;
t81 = -t111 * t155 - t210;
t55 = -pkin(7) * t101 + t81;
t54 = qJD(3) * pkin(3) + t55;
t30 = t236 * t54 - t211;
t254 = qJD(5) - t30;
t23 = -pkin(4) * t148 + t254;
t165 = t152 * t169;
t159 = t101 * t236 - t165;
t92 = t236 * t169;
t72 = t101 * t152 + t92;
t231 = Ifges(6,5) * t72;
t68 = Ifges(5,4) * t72;
t257 = t148 * t259 + t159 * t260 + t231 - t68;
t135 = pkin(2) * t150 + pkin(1);
t91 = pkin(3) * t173 - t135;
t83 = qJD(1) * t91 + qJD(2);
t32 = pkin(4) * t72 - qJ(5) * t159 + t83;
t276 = t83 * mrSges(5,2) + t23 * mrSges(6,2) - t30 * mrSges(5,3) - t32 * mrSges(6,3) + t257 / 0.2e1;
t272 = Ifges(5,4) - Ifges(6,5);
t147 = pkin(8) + qJ(3);
t140 = sin(t147);
t141 = cos(t147);
t142 = qJ(4) + t147;
t132 = sin(t142);
t133 = cos(t142);
t251 = -t277 * t133 + (mrSges(5,2) - mrSges(6,3)) * t132;
t275 = -mrSges(4,1) * t141 + mrSges(4,2) * t140 + t251;
t80 = t110 * t236 - t152 * t173;
t241 = t80 / 0.2e1;
t274 = -m(6) - m(5);
t271 = Ifges(5,6) - Ifges(6,6);
t234 = mrSges(5,3) * t72;
t58 = -mrSges(5,2) * t148 - t234;
t61 = -mrSges(6,2) * t72 + mrSges(6,3) * t148;
t256 = -t58 - t61;
t233 = mrSges(5,3) * t159;
t255 = -mrSges(6,2) * t159 + t148 * t277 - t233;
t202 = qJD(4) * t152;
t195 = t236 * t56;
t33 = t152 * t55 + t195;
t270 = -pkin(3) * t202 + t33;
t267 = t133 * (-m(6) * qJ(5) - mrSges(6,3));
t43 = pkin(4) * t159 + qJ(5) * t72;
t235 = m(3) * qJ(2);
t266 = -m(4) * t218 + mrSges(2,2) - mrSges(6,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - t235;
t179 = -mrSges(3,1) * t150 + mrSges(3,2) * t149;
t265 = m(3) * pkin(1) + m(4) * t135 + mrSges(2,1) - t179 - t275;
t31 = t152 * t54 + t195;
t28 = qJ(5) * t148 + t31;
t67 = Ifges(6,5) * t159;
t39 = Ifges(6,6) * t148 + Ifges(6,3) * t72 + t67;
t232 = Ifges(5,4) * t159;
t40 = -Ifges(5,2) * t72 + Ifges(5,6) * t148 + t232;
t263 = t31 * mrSges(5,3) + t28 * mrSges(6,2) - t39 / 0.2e1 + t40 / 0.2e1 - t32 * mrSges(6,1) - t83 * mrSges(5,1);
t85 = -t120 * t153 + t121 * t155;
t252 = pkin(4) * t133 + qJ(5) * t132;
t154 = sin(qJ(1));
t156 = cos(qJ(1));
t250 = -g(1) * t154 + g(2) * t156;
t249 = g(1) * t156 + g(2) * t154;
t183 = pkin(6) * qJDD(1) + t126;
t95 = t183 * t149;
t96 = t183 * t150;
t47 = -qJD(3) * t82 - t153 * t96 - t155 * t95;
t102 = t173 * qJD(3);
t77 = -qJD(1) * t102 + qJDD(1) * t110;
t25 = qJDD(3) * pkin(3) - pkin(7) * t77 + t47;
t203 = qJD(3) * t155;
t46 = -qJD(3) * t210 - t111 * t203 - t153 * t95 + t155 * t96;
t103 = t110 * qJD(3);
t78 = -qJD(1) * t103 - qJDD(1) * t173;
t29 = pkin(7) * t78 + t46;
t6 = -qJD(4) * t31 - t152 * t29 + t236 * t25;
t246 = -t72 / 0.2e1;
t245 = t72 / 0.2e1;
t243 = -t159 / 0.2e1;
t242 = t159 / 0.2e1;
t240 = t101 / 0.2e1;
t238 = -t148 / 0.2e1;
t237 = t148 / 0.2e1;
t230 = pkin(3) * t101;
t229 = pkin(3) * t103;
t228 = pkin(3) * t140;
t130 = pkin(3) * t141;
t227 = pkin(3) * t152;
t217 = mrSges(5,2) * t133;
t216 = mrSges(4,3) * t101;
t215 = Ifges(4,4) * t101;
t213 = t145 * mrSges(3,3);
t212 = t146 * mrSges(3,3);
t199 = qJDD(1) * t149;
t198 = qJDD(1) * t150;
t197 = t236 * pkin(3);
t190 = -t78 * mrSges(4,1) + mrSges(4,2) * t77;
t26 = t236 * t77 + (-qJD(4) * t101 + t78) * t152 - qJD(4) * t92;
t187 = qJD(4) * t236;
t27 = -qJD(4) * t165 + t101 * t187 + t152 * t77 - t236 * t78;
t189 = t27 * mrSges(5,1) + mrSges(5,2) * t26;
t188 = mrSges(6,1) * t27 - t26 * mrSges(6,3);
t186 = t154 * t267;
t185 = t156 * t267;
t143 = qJDD(3) + qJDD(4);
t12 = -t143 * mrSges(6,1) + mrSges(6,2) * t26;
t84 = -t120 * t155 - t121 * t153;
t181 = pkin(3) * t187;
t180 = -mrSges(3,1) * t198 + mrSges(3,2) * t199;
t65 = -pkin(7) * t110 + t84;
t66 = -pkin(7) * t173 + t85;
t171 = -t152 * t66 + t236 * t65;
t38 = t152 * t65 + t236 * t66;
t5 = t152 * t25 + t187 * t54 - t202 * t56 + t236 * t29;
t170 = -t110 * t152 - t173 * t236;
t114 = -qJDD(1) * t135 + qJDD(2);
t166 = mrSges(4,3) * t169;
t162 = t217 + (m(6) * pkin(4) + t277) * t132;
t161 = m(6) * (-pkin(4) * t132 - t228) - t132 * mrSges(6,1);
t57 = -pkin(3) * t78 + t114;
t62 = -t120 * t203 + qJD(2) * t207 + (-qJD(2) * t149 - qJD(3) * t121) * t153;
t2 = qJ(5) * t143 + qJD(5) * t148 + t5;
t3 = -pkin(4) * t143 + qJDD(5) - t6;
t160 = t6 * mrSges(5,1) - t3 * mrSges(6,1) - t5 * mrSges(5,2) + t2 * mrSges(6,3) - t271 * t27 + t259 * t26 + (Ifges(6,2) + Ifges(5,3)) * t143;
t63 = -qJD(2) * t110 - qJD(3) * t85;
t158 = pkin(7) * t102 + t63;
t144 = -pkin(7) - t218;
t139 = -qJDD(1) * pkin(1) + qJDD(2);
t136 = -t197 - pkin(4);
t134 = qJ(5) + t227;
t127 = t181 + qJD(5);
t115 = -qJD(1) * t135 + qJD(2);
t113 = t130 + t135;
t97 = Ifges(4,4) * t169;
t87 = qJD(3) * mrSges(4,1) - t216;
t86 = -qJD(3) * mrSges(4,2) - t166;
t70 = t101 * Ifges(4,1) + Ifges(4,5) * qJD(3) - t97;
t69 = -Ifges(4,2) * t169 + Ifges(4,6) * qJD(3) + t215;
t51 = -pkin(7) * t103 + t62;
t49 = qJD(4) * t80 - t102 * t152 + t103 * t236;
t48 = qJD(4) * t170 - t102 * t236 - t103 * t152;
t45 = mrSges(5,1) * t72 + mrSges(5,2) * t159;
t44 = mrSges(6,1) * t72 - mrSges(6,3) * t159;
t36 = -pkin(4) * t170 - qJ(5) * t80 + t91;
t35 = t230 + t43;
t34 = t236 * t55 - t211;
t14 = -mrSges(6,2) * t27 + mrSges(6,3) * t143;
t13 = -mrSges(5,2) * t143 - mrSges(5,3) * t27;
t11 = mrSges(5,1) * t143 - mrSges(5,3) * t26;
t10 = pkin(4) * t49 - qJ(5) * t48 - qJD(5) * t80 + t229;
t7 = pkin(4) * t27 - qJ(5) * t26 - qJD(5) * t159 + t57;
t1 = [((-Ifges(6,3) - Ifges(5,2)) * t170 - 0.2e1 * t272 * t241) * t27 + t10 * t44 + (-m(5) * t30 + m(6) * t23 - t255) * (qJD(4) * t38 + t152 * t51 - t158 * t236) + (m(5) * t31 + m(6) * t28 - t256) * (qJD(4) * t171 + t152 * t158 + t236 * t51) + (t271 * t170 + t259 * t80) * t143 / 0.2e1 + (-Ifges(5,2) * t246 + Ifges(6,3) * t245 - t237 * t271 - t242 * t272 - t263) * t49 + (t274 * (t113 * t156 - t144 * t154) + (-m(6) * t252 - t265) * t156 + t266 * t154) * g(2) + ((-t144 * t274 + t266) * t156 + (-m(6) * (-t113 - t252) + m(5) * t113 + t265) * t154) * g(1) + (t143 * t259 + t26 * t260) * t241 + (t170 * t272 + t260 * t80) * t26 / 0.2e1 + m(5) * (t229 * t83 + t57 * t91) + m(6) * (t10 * t32 + t36 * t7) + (m(5) * t6 - m(6) * t3 + t11 - t12) * t171 - t170 * (Ifges(6,5) * t26 + Ifges(6,6) * t143) / 0.2e1 + (t170 * t2 + t3 * t80) * mrSges(6,2) + t170 * (Ifges(5,4) * t26 + Ifges(5,6) * t143) / 0.2e1 + t7 * (-mrSges(6,1) * t170 - mrSges(6,3) * t80) + t57 * (-mrSges(5,1) * t170 + mrSges(5,2) * t80) + (t170 * t5 - t6 * t80) * mrSges(5,3) + (mrSges(4,1) * t114 - mrSges(4,3) * t46 - Ifges(4,4) * t77 - Ifges(4,2) * t78 - Ifges(4,6) * qJDD(3)) * t173 + t126 * t212 + t126 * t213 + (Ifges(3,4) * t149 + Ifges(3,2) * t150) * t198 + (Ifges(3,1) * t149 + Ifges(3,4) * t150) * t199 + (mrSges(4,2) * t114 - mrSges(4,3) * t47 + Ifges(4,1) * t77 + Ifges(4,4) * t78 + Ifges(4,5) * qJDD(3)) * t110 + (-Ifges(4,1) * t102 - Ifges(4,4) * t103) * t240 + qJD(3) * (-Ifges(4,5) * t102 - Ifges(4,6) * t103) / 0.2e1 + t115 * (mrSges(4,1) * t103 - mrSges(4,2) * t102) - (-Ifges(4,4) * t102 - Ifges(4,2) * t103) * t169 / 0.2e1 + t204 * t126 * mrSges(3,3) + m(4) * (-t114 * t135 + t46 * t85 + t47 * t84 + t62 * t82 + t63 * t81) + (m(5) * t5 + m(6) * t2 + t13 + t14) * t38 + (t102 * t81 - t103 * t82 - t77 * t84 + t78 * t85) * mrSges(4,3) + t62 * t86 + t63 * t87 - t102 * t70 / 0.2e1 - t103 * t69 / 0.2e1 + t139 * t179 - pkin(1) * t180 + t36 * t188 + t91 * t189 - t135 * t190 + t45 * t229 + Ifges(2,3) * qJDD(1) + m(3) * (-pkin(1) * t139 + (t126 + t201) * qJ(2) * t204) + (mrSges(4,1) * t84 - mrSges(4,2) * t85) * qJDD(3) + (Ifges(5,4) * t246 + Ifges(6,5) * t245 + t237 * t259 + t242 * t260 + t276) * t48; t101 * t87 + t86 * t169 + t180 + t188 + t189 + t190 - t256 * t72 + t255 * t159 + (-t204 * t235 - t212 - t213) * qJD(1) ^ 2 + (t139 + t250) * m(3) + (-t159 * t23 + t28 * t72 + t250 + t7) * m(6) + (t159 * t30 + t31 * t72 + t250 + t57) * m(5) + (t81 * t101 + t169 * t82 + t114 + t250) * m(4); -t35 * t44 - t46 * mrSges(4,2) + t47 * mrSges(4,1) + (-m(6) * (t130 + t252) - m(5) * t130 + t275) * g(3) + ((t236 * t6 + t152 * t5 + (-t152 * t30 + t236 * t31) * qJD(4)) * pkin(3) - t230 * t83 + t30 * t33 - t31 * t34) * m(5) + t270 * t255 + (t134 * t2 + t136 * t3 - t32 * t35 + (-t34 + t127) * t28 - t270 * t23) * m(6) + t58 * t181 + (t216 + t87) * t82 + (-Ifges(4,2) * t101 + t70 - t97) * t169 / 0.2e1 + (-t166 - t86) * t81 + t11 * t197 + (m(5) * t228 + mrSges(4,1) * t140 + mrSges(5,1) * t132 + mrSges(4,2) * t141 + t217) * t249 - t115 * (t101 * mrSges(4,1) - mrSges(4,2) * t169) - qJD(3) * (-Ifges(4,5) * t169 - Ifges(4,6) * t101) / 0.2e1 - t101 * (-Ifges(4,1) * t169 - t215) / 0.2e1 + t69 * t240 + t13 * t227 + t160 + (-Ifges(5,2) * t245 + Ifges(6,3) * t246 - t238 * t271 - t243 * t272 + t263) * t159 + Ifges(4,5) * t77 + Ifges(4,6) * t78 + t127 * t61 + t134 * t14 + t136 * t12 + Ifges(4,3) * qJDD(3) - g(1) * (t156 * t161 - t185) - g(2) * (t154 * t161 - t186) - t45 * t230 + t256 * t34 + (-Ifges(5,4) * t245 - Ifges(6,5) * t246 - t238 * t259 - t243 * t260 + t276) * t72; qJD(5) * t61 - t43 * t44 - pkin(4) * t12 + qJ(5) * t14 + (t159 * t28 + t23 * t72) * mrSges(6,2) + t40 * t242 + t160 - t32 * (mrSges(6,1) * t159 + mrSges(6,3) * t72) - t83 * (mrSges(5,1) * t159 - mrSges(5,2) * t72) + (Ifges(6,3) * t159 - t231) * t246 + (t156 * t162 + t185) * g(1) + (t154 * t162 + t186) * g(2) + t251 * g(3) + (t233 + t255) * t31 + (-t234 + t256) * t30 + (-t159 * t271 - t259 * t72) * t238 + (-pkin(4) * t3 - g(3) * t252 + qJ(5) * t2 - t23 * t31 + t254 * t28 - t32 * t43) * m(6) + (-Ifges(5,2) * t159 + t257 - t68) * t245 + (-t260 * t72 - t232 + t39 + t67) * t243; -t148 * t61 + t159 * t44 + (g(3) * t133 - t132 * t249 - t28 * t148 + t159 * t32 + t3) * m(6) + t12;];
tau = t1;
