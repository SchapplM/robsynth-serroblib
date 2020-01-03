% Calculate vector of inverse dynamics joint torques for
% S5RPRRP5
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:42
% EndTime: 2019-12-31 18:40:47
% DurationCPUTime: 2.58s
% Computational Cost: add. (1904->271), mult. (3318->338), div. (0->0), fcn. (1677->12), ass. (0->128)
t253 = mrSges(5,1) + mrSges(6,1);
t126 = sin(qJ(4));
t129 = cos(qJ(4));
t122 = qJD(1) + qJD(3);
t186 = t122 * t129;
t173 = mrSges(5,3) * t186;
t174 = mrSges(6,2) * t186;
t80 = qJD(4) * mrSges(6,3) + t174;
t234 = -qJD(4) * mrSges(5,2) + t173 + t80;
t187 = t122 * t126;
t175 = mrSges(5,3) * t187;
t176 = mrSges(6,2) * t187;
t250 = t253 * qJD(4) - t175 - t176;
t252 = -t126 * t250 + t234 * t129;
t247 = Ifges(6,4) + Ifges(5,5);
t246 = Ifges(6,6) - Ifges(5,6);
t155 = t129 * mrSges(6,1) + t126 * mrSges(6,3);
t91 = -t129 * mrSges(5,1) + mrSges(5,2) * t126;
t249 = t91 - t155;
t248 = -mrSges(6,2) - mrSges(5,3) + mrSges(4,2);
t204 = Ifges(6,5) * t129;
t153 = Ifges(6,1) * t126 - t204;
t97 = Ifges(5,4) * t186;
t245 = Ifges(5,1) * t187 + t247 * qJD(4) + t122 * t153 + t97;
t125 = cos(pkin(8));
t112 = pkin(1) * t125 + pkin(2);
t124 = sin(pkin(8));
t219 = pkin(1) * t124;
t177 = qJD(1) * t219;
t244 = -qJD(3) * t177 + t112 * qJDD(1);
t154 = t126 * mrSges(6,1) - t129 * mrSges(6,3);
t156 = mrSges(5,1) * t126 + mrSges(5,2) * t129;
t127 = sin(qJ(3));
t130 = cos(qJ(3));
t90 = t112 * qJD(1);
t43 = -t127 * t177 + t130 * t90;
t149 = pkin(4) * t129 + qJ(5) * t126;
t86 = -pkin(3) - t149;
t23 = t122 * t86 - t43;
t36 = -pkin(3) * t122 - t43;
t243 = t23 * t154 + t36 * t156;
t242 = t246 * t126 + t247 * t129;
t192 = pkin(1) * qJDD(1);
t171 = t124 * t192;
t241 = qJD(3) * t90 + t171;
t121 = qJDD(1) + qJDD(3);
t16 = t244 * t127 + t241 * t130;
t13 = pkin(7) * t121 + t16;
t181 = qJD(4) * t129;
t172 = qJD(2) * t181 + t126 * qJDD(2) + t129 * t13;
t182 = qJD(4) * t126;
t44 = t127 * t90 + t130 * t177;
t37 = pkin(7) * t122 + t44;
t7 = -t182 * t37 + t172;
t31 = qJD(2) * t126 + t129 * t37;
t8 = -qJD(4) * t31 + qJDD(2) * t129 - t126 * t13;
t240 = -t8 * t126 + t129 * t7;
t201 = t126 * t37;
t3 = qJDD(4) * qJ(5) + (qJD(5) - t201) * qJD(4) + t172;
t5 = -qJDD(4) * pkin(4) + qJDD(5) - t8;
t239 = t126 * t5 + t129 * t3;
t69 = t121 * t126 + t122 * t181;
t47 = -qJDD(4) * mrSges(6,1) + t69 * mrSges(6,2);
t236 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t69 - t47;
t68 = -t129 * t121 + t122 * t182;
t48 = -mrSges(6,2) * t68 + qJDD(4) * mrSges(6,3);
t237 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t68 + t48;
t30 = qJD(2) * t129 - t201;
t24 = -qJD(4) * pkin(4) + qJD(5) - t30;
t26 = qJD(4) * qJ(5) + t31;
t238 = -t126 * t236 + t237 * t129 + m(5) * ((-t126 * t31 - t129 * t30) * qJD(4) + t240) + m(6) * ((-t126 * t26 + t129 * t24) * qJD(4) + t239) - t234 * t182 - t250 * t181;
t123 = qJ(1) + pkin(8);
t119 = qJ(3) + t123;
t111 = cos(t119);
t189 = t111 * t129;
t190 = t111 * t126;
t233 = pkin(4) * t189 + qJ(5) * t190;
t62 = t112 * t130 - t127 * t219;
t110 = sin(t119);
t232 = g(1) * t111 + g(2) * t110;
t231 = t248 * t111 + (-m(6) * t86 + mrSges(4,1) - t249) * t110;
t230 = -t111 * mrSges(4,1) + (mrSges(5,2) - mrSges(6,3)) * t190 - t253 * t189 + t248 * t110;
t17 = -t241 * t127 + t244 * t130;
t14 = -pkin(3) * t121 - t17;
t229 = m(5) * t14 + mrSges(5,1) * t68 + mrSges(5,2) * t69;
t226 = m(5) * t36 + (-mrSges(4,1) + t91) * t122;
t146 = -t126 * t30 + t129 * t31;
t147 = t126 * t24 + t129 * t26;
t224 = -m(5) * t146 - m(6) * t147 + t122 * mrSges(4,2) - t252;
t221 = m(3) + m(4);
t220 = t126 / 0.2e1;
t128 = sin(qJ(1));
t218 = pkin(1) * t128;
t217 = pkin(3) * t110;
t131 = cos(qJ(1));
t120 = t131 * pkin(1);
t207 = Ifges(5,4) * t126;
t206 = Ifges(5,4) * t129;
t205 = Ifges(6,5) * t126;
t63 = t127 * t112 + t130 * t219;
t185 = t111 * pkin(3) + t110 * pkin(7);
t117 = cos(t123);
t184 = pkin(2) * t117 + t120;
t183 = qJD(4) * t122;
t180 = qJD(5) * t126;
t164 = -t183 / 0.2e1;
t159 = t184 + t185;
t116 = sin(t123);
t158 = -pkin(2) * t116 - t218;
t152 = Ifges(5,2) * t129 + t207;
t148 = pkin(4) * t126 - qJ(5) * t129;
t104 = t111 * pkin(7);
t142 = t104 + t158;
t139 = t126 * (Ifges(5,1) * t129 - t207);
t138 = t129 * (Ifges(6,3) * t126 + t204);
t64 = pkin(4) * t182 - qJ(5) * t181 - t180;
t54 = t63 * qJD(3);
t96 = Ifges(6,5) * t187;
t55 = Ifges(6,6) * qJD(4) - Ifges(6,3) * t186 + t96;
t56 = Ifges(5,6) * qJD(4) + t122 * t152;
t9 = pkin(4) * t68 - qJ(5) * t69 - t122 * t180 + t14;
t133 = -t9 * t155 + t129 * (Ifges(5,4) * t69 + Ifges(5,6) * qJDD(4)) / 0.2e1 - t129 * (Ifges(6,5) * t69 + Ifges(6,6) * qJDD(4)) / 0.2e1 + Ifges(4,3) * t121 + t14 * t91 - t16 * mrSges(4,2) + t17 * mrSges(4,1) + t138 * t164 - t56 * t182 / 0.2e1 + t139 * t183 / 0.2e1 + (Ifges(5,1) * t126 + t153 + t206) * t69 / 0.2e1 + ((Ifges(5,1) + Ifges(6,1)) * t69 + t247 * qJDD(4)) * t220 + (t247 * t126 - t246 * t129) * qJDD(4) / 0.2e1 + (t122 * (Ifges(6,1) * t129 + t205) + t55) * t182 / 0.2e1 + (-t181 * t30 - t182 * t31 + t240) * mrSges(5,3) + (-t152 / 0.2e1 + t205 / 0.2e1 + (Ifges(6,5) - Ifges(5,4)) * t220 + (-Ifges(5,2) / 0.2e1 - Ifges(6,3)) * t129) * t68 + (t122 * (-Ifges(5,2) * t126 + t206) + t245) * t181 / 0.2e1 + (t24 * t181 - t26 * t182 + t239) * mrSges(6,2) + (t243 + t242 * qJD(4) / 0.2e1) * qJD(4);
t67 = t148 * t122;
t65 = t155 * t122;
t38 = -t62 + t86;
t28 = mrSges(6,1) * t68 - mrSges(6,3) * t69;
t27 = t54 + t64;
t1 = [m(6) * (t23 * t27 + t38 * t9) + t133 + m(4) * (t16 * t63 + t17 * t62) - t27 * t65 + t38 * t28 + 0.2e1 * t125 * mrSges(3,1) * t192 - 0.2e1 * mrSges(3,2) * t171 + t229 * (-pkin(3) - t62) + (-m(4) * t43 + t226) * t54 + (t62 * mrSges(4,1) - t63 * mrSges(4,2)) * t121 + (m(3) * (t124 ^ 2 + t125 ^ 2) * pkin(1) ^ 2 + Ifges(3,3) + Ifges(2,3)) * qJDD(1) + (m(4) * t44 - t224) * t62 * qJD(3) + (-mrSges(2,1) * t131 + mrSges(2,2) * t128 - m(5) * t159 - m(6) * (t159 + t233) - m(4) * t184 - m(3) * t120 - mrSges(3,1) * t117 + mrSges(3,2) * t116 + t230) * g(2) + (-m(4) * t158 - m(6) * t142 + mrSges(2,1) * t128 + mrSges(2,2) * t131 - m(5) * (t142 - t217) + m(3) * t218 + mrSges(3,1) * t116 + mrSges(3,2) * t117 + t231) * g(1) + t238 * (pkin(7) + t63); t236 * t129 + t237 * t126 + t221 * qJDD(2) + t252 * qJD(4) + m(5) * (qJD(4) * t146 + t126 * t7 + t129 * t8) + m(6) * (qJD(4) * t147 + t126 * t3 - t129 * t5) + (-m(5) - m(6) - t221) * g(3); t133 + t86 * t28 - t64 * t65 + m(6) * (t23 * t64 + t86 * t9) - t229 * pkin(3) + (-m(6) * (t185 + t233) - m(5) * t185 + t230) * g(2) + (-m(6) * t104 - m(5) * (t104 - t217) + t231) * g(1) + (-m(6) * t23 - t226 + t65) * t44 + t224 * t43 + t238 * pkin(7); t246 * t68 + t247 * t69 - (-Ifges(5,2) * t187 + t245 + t97) * t186 / 0.2e1 - (Ifges(6,1) * t186 + t55 + t96) * t187 / 0.2e1 + t242 * t164 + qJD(5) * t80 + t67 * t65 + (-pkin(4) * t5 - g(3) * t149 + qJ(5) * t3 + qJD(5) * t26 - t23 * t67) * m(6) + t249 * g(3) - pkin(4) * t47 + qJ(5) * t48 - t5 * mrSges(6,1) - t7 * mrSges(5,2) + t8 * mrSges(5,1) + t3 * mrSges(6,3) + (Ifges(6,2) + Ifges(5,3)) * qJDD(4) + (-m(6) * t26 + t173 - t234) * t30 + (-m(6) * t24 + t175 + t250) * t31 + t26 * t176 - t24 * t174 + t56 * t187 / 0.2e1 + (m(6) * t148 + t154 + t156) * t232 + ((-t139 / 0.2e1 + t138 / 0.2e1) * t122 - t243) * t122; -t65 * t187 - qJD(4) * t80 + (g(3) * t129 - t26 * qJD(4) - t126 * t232 + t23 * t187 + t5) * m(6) + t47;];
tau = t1;
