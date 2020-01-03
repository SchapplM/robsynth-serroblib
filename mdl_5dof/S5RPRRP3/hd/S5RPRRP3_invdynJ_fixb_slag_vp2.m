% Calculate vector of inverse dynamics joint torques for
% S5RPRRP3
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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:46:55
% EndTime: 2020-01-03 11:47:09
% DurationCPUTime: 6.50s
% Computational Cost: add. (2644->349), mult. (5640->441), div. (0->0), fcn. (3499->12), ass. (0->157)
t237 = Ifges(5,4) + Ifges(6,4);
t153 = sin(qJ(3));
t156 = cos(qJ(3));
t123 = -mrSges(4,1) * t156 + mrSges(4,2) * t153;
t149 = qJ(3) + qJ(4);
t140 = sin(t149);
t141 = cos(t149);
t206 = mrSges(5,2) + mrSges(6,2);
t207 = mrSges(5,1) + mrSges(6,1);
t169 = t140 * t206 - t207 * t141;
t247 = t123 + t169;
t238 = Ifges(5,1) + Ifges(6,1);
t236 = Ifges(5,5) + Ifges(6,5);
t235 = Ifges(5,2) + Ifges(6,2);
t234 = Ifges(5,6) + Ifges(6,6);
t152 = sin(qJ(4));
t155 = cos(qJ(4));
t106 = -t152 * t153 + t155 * t156;
t100 = t106 * qJD(1);
t246 = t100 / 0.2e1;
t147 = qJD(3) + qJD(4);
t245 = t147 / 0.2e1;
t244 = t237 * t100;
t179 = t206 * t141;
t107 = t152 * t156 + t153 * t155;
t101 = t107 * qJD(1);
t243 = t237 * t101;
t242 = t235 * t100 + t234 * t147 + t243;
t241 = t238 * t101 + t236 * t147 + t244;
t230 = -m(3) - m(4) - m(5) - m(6);
t240 = pkin(1) * t230 - mrSges(2,1);
t212 = t153 / 0.2e1;
t150 = sin(pkin(8));
t126 = pkin(1) * t150 + pkin(6);
t119 = t126 * qJD(1);
t176 = pkin(7) * qJD(1) + t119;
t187 = qJD(2) * t153;
t75 = t156 * t176 + t187;
t69 = t152 * t75;
t139 = t156 * qJD(2);
t74 = -t153 * t176 + t139;
t72 = qJD(3) * pkin(3) + t74;
t21 = t155 * t72 - t69;
t93 = t101 * qJ(5);
t13 = t21 - t93;
t202 = mrSges(6,3) * t100;
t78 = -mrSges(6,2) * t147 + t202;
t203 = mrSges(5,3) * t100;
t79 = -mrSges(5,2) * t147 + t203;
t233 = t78 + t79;
t80 = mrSges(6,1) * t147 - mrSges(6,3) * t101;
t81 = mrSges(5,1) * t147 - mrSges(5,3) * t101;
t232 = t80 + t81;
t151 = cos(pkin(8));
t127 = -pkin(1) * t151 - pkin(2);
t143 = t156 * pkin(3);
t116 = t127 - t143;
t117 = t126 * qJDD(1);
t186 = qJD(3) * t153;
t55 = qJD(3) * t139 + t153 * qJDD(2) + t156 * t117 - t119 * t186;
t87 = t119 * t156 + t187;
t56 = -qJD(3) * t87 + t156 * qJDD(2) - t117 * t153;
t231 = -t153 * t56 + t156 * t55;
t194 = qJ(5) * t100;
t71 = t155 * t75;
t22 = t152 * t72 + t71;
t14 = t22 + t194;
t174 = t22 * mrSges(5,3) + t14 * mrSges(6,3);
t190 = pkin(4) * t141 + t143;
t229 = -mrSges(3,1) - m(6) * (pkin(2) + t190) - m(5) * (t143 + pkin(2)) - m(4) * pkin(2) + t247;
t173 = mrSges(4,1) * t153 + mrSges(4,2) * t156;
t227 = t127 * qJD(1) * t173 + qJD(3) * (Ifges(4,5) * t156 - Ifges(4,6) * t153) / 0.2e1;
t158 = -pkin(7) - pkin(6);
t226 = m(4) * pkin(6) - m(5) * t158 - m(6) * (-qJ(5) + t158) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t218 = t101 / 0.2e1;
t148 = qJ(1) + pkin(8);
t135 = sin(t148);
t210 = g(2) * t135;
t205 = pkin(7) + t126;
t145 = qJDD(3) + qJDD(4);
t182 = qJD(1) * qJD(3);
t109 = qJDD(1) * t156 - t153 * t182;
t110 = qJDD(1) * t153 + t156 * t182;
t166 = t107 * qJD(4);
t48 = -qJD(1) * t166 + t109 * t155 - t110 * t152;
t29 = -mrSges(6,2) * t145 + mrSges(6,3) * t48;
t30 = -mrSges(5,2) * t145 + mrSges(5,3) * t48;
t204 = t29 + t30;
t34 = t155 * t74 - t69;
t103 = t205 * t153;
t104 = t205 * t156;
t59 = -t152 * t103 + t155 * t104;
t201 = Ifges(4,4) * t153;
t200 = Ifges(4,4) * t156;
t196 = t156 * Ifges(4,2);
t136 = cos(t148);
t192 = t136 * t140;
t189 = qJD(1) * t153;
t188 = qJD(1) * t156;
t185 = qJD(3) * t156;
t184 = qJD(4) * t152;
t183 = qJD(4) * t155;
t181 = pkin(3) * t186;
t165 = t106 * qJD(4);
t47 = qJD(1) * t165 + t109 * t152 + t110 * t155;
t180 = -t48 * mrSges(6,1) + t47 * mrSges(6,2);
t33 = -t152 * t74 - t71;
t178 = qJD(3) * t205;
t58 = -t155 * t103 - t104 * t152;
t175 = -t136 * t179 - t207 * t192;
t118 = t127 * qJDD(1);
t172 = t196 + t201;
t86 = -t119 * t153 + t139;
t170 = t153 * t87 + t156 * t86;
t32 = qJDD(3) * pkin(3) - pkin(7) * t110 + t56;
t35 = pkin(7) * t109 + t55;
t5 = t152 * t32 + t155 * t35 + t72 * t183 - t184 * t75;
t94 = t153 * t178;
t95 = t156 * t178;
t11 = -t103 * t183 - t104 * t184 - t152 * t95 - t155 * t94;
t167 = t153 * (Ifges(4,1) * t156 - t201);
t102 = t116 * qJD(1);
t162 = m(6) * (-pkin(3) * t153 - pkin(4) * t140) - t173;
t76 = -pkin(3) * t109 + t118;
t6 = -qJD(4) * t22 - t152 * t35 + t155 * t32;
t12 = -qJD(4) * t59 + t152 * t94 - t155 * t95;
t10 = pkin(4) * t147 + t13;
t2 = pkin(4) * t145 - qJ(5) * t47 - qJD(5) * t101 + t6;
t3 = qJ(5) * t48 + qJD(5) * t100 + t5;
t64 = -pkin(4) * t100 + qJD(5) + t102;
t160 = t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2) + t10 * t202 - t102 * (mrSges(5,1) * t101 + mrSges(5,2) * t100) - t64 * (mrSges(6,1) * t101 + mrSges(6,2) * t100) + t21 * t203 + t234 * t48 + t236 * t47 - (t238 * t100 - t243) * t101 / 0.2e1 + t242 * t218 - (t236 * t100 - t234 * t101) * t147 / 0.2e1 + (Ifges(6,3) + Ifges(5,3)) * t145 - (-t235 * t101 + t241 + t244) * t100 / 0.2e1;
t157 = cos(qJ(1));
t154 = sin(qJ(1));
t132 = Ifges(4,4) * t188;
t122 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t188;
t120 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t189;
t99 = Ifges(4,1) * t189 + Ifges(4,5) * qJD(3) + t132;
t98 = Ifges(4,6) * qJD(3) + qJD(1) * t172;
t92 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t110;
t91 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t109;
t77 = pkin(3) * t189 + pkin(4) * t101;
t73 = -pkin(4) * t106 + t116;
t68 = -qJD(3) * t107 - t166;
t67 = qJD(3) * t106 + t165;
t63 = -mrSges(5,1) * t100 + mrSges(5,2) * t101;
t62 = -mrSges(6,1) * t100 + mrSges(6,2) * t101;
t57 = -pkin(4) * t68 + t181;
t41 = qJ(5) * t106 + t59;
t40 = -qJ(5) * t107 + t58;
t28 = mrSges(5,1) * t145 - mrSges(5,3) * t47;
t27 = mrSges(6,1) * t145 - mrSges(6,3) * t47;
t18 = -t93 + t34;
t17 = t33 - t194;
t16 = -pkin(4) * t48 + qJDD(5) + t76;
t8 = -qJ(5) * t67 - qJD(5) * t107 + t12;
t7 = qJ(5) * t68 + qJD(5) * t106 + t11;
t1 = [(m(4) * t127 + t123) * t118 + t227 * qJD(3) + (Ifges(4,1) * t110 + Ifges(4,4) * t109) * t212 + (0.2e1 * Ifges(4,5) * t212 + Ifges(4,6) * t156) * qJDD(3) + (-t185 * t86 - t186 * t87 + t231) * mrSges(4,3) + (m(4) * (-qJD(3) * t170 + t231) + t156 * t91 - t153 * t92 - t120 * t185 - t122 * t186) * t126 + (-t76 * mrSges(5,1) - t16 * mrSges(6,1) + t5 * mrSges(5,3) + t3 * mrSges(6,3) + t234 * t145 + t235 * t48 + t237 * t47) * t106 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t151 - 0.2e1 * mrSges(3,2) * t150 + m(3) * (t150 ^ 2 + t151 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (t156 * (-Ifges(4,2) * t153 + t200) + t167) * t182 / 0.2e1 + (-mrSges(2,2) * t157 + t229 * t135 + t226 * t136 + t240 * t154) * g(3) + (mrSges(2,2) * t154 - t226 * t135 + t229 * t136 + t240 * t157) * g(2) + (t76 * mrSges(5,2) + t16 * mrSges(6,2) - t6 * mrSges(5,3) - t2 * mrSges(6,3) + t236 * t145 + t237 * t48 + t238 * t47) * t107 + t110 * (t153 * Ifges(4,1) + t200) / 0.2e1 + t99 * t185 / 0.2e1 - t98 * t186 / 0.2e1 + m(5) * (t102 * t181 + t11 * t22 + t116 * t76 + t12 * t21 + t5 * t59 + t58 * t6) + t73 * t180 + t109 * t172 / 0.2e1 + t127 * (-mrSges(4,1) * t109 + mrSges(4,2) * t110) + t116 * (-mrSges(5,1) * t48 + mrSges(5,2) * t47) + t7 * t78 + t11 * t79 + t8 * t80 + t12 * t81 + t63 * t181 + t58 * t28 + t59 * t30 + t57 * t62 + t40 * t27 + t41 * t29 + m(6) * (t10 * t8 + t14 * t7 + t16 * t73 + t2 * t40 + t3 * t41 + t57 * t64) + (t102 * mrSges(5,2) + t64 * mrSges(6,2) + t236 * t245 + t237 * t246 + t238 * t218 + t241 / 0.2e1 - t21 * mrSges(5,3) - t10 * mrSges(6,3)) * t67 + (t234 * t245 + t235 * t246 + t237 * t218 - t102 * mrSges(5,1) - t64 * mrSges(6,1) + t174 + t242 / 0.2e1) * t68 + t156 * (Ifges(4,4) * t110 + Ifges(4,2) * t109) / 0.2e1; m(3) * qJDD(2) + t153 * t91 + t156 * t92 + t232 * t68 + t233 * t67 + t204 * t107 + (t27 + t28) * t106 + (-t120 * t153 + t122 * t156) * qJD(3) + m(4) * (t153 * t55 + t156 * t56 + (-t153 * t86 + t156 * t87) * qJD(3)) + m(5) * (t106 * t6 + t107 * t5 + t21 * t68 + t22 * t67) + m(6) * (t10 * t68 + t106 * t2 + t107 * t3 + t14 * t67) + t230 * g(1); (t98 * t212 + (t196 * t212 - t167 / 0.2e1) * qJD(1) + t170 * mrSges(4,3) - (t132 + t99) * t156 / 0.2e1 - t227) * qJD(1) + (-m(6) * t190 + t247) * g(1) + (t136 * t162 + t175) * g(3) + (t140 * t207 - t162 + t179) * t210 + t160 - m(6) * (t10 * t17 + t14 * t18 + t64 * t77) + Ifges(4,5) * t110 + t87 * t120 - t86 * t122 + Ifges(4,6) * t109 - t77 * t62 - t18 * t78 - t34 * t79 - t17 * t80 - t33 * t81 + t56 * mrSges(4,1) + (-t63 * t189 + t155 * t28 + (m(6) * t3 + t204) * t152 + ((m(6) * t14 + t233) * t155 + (-m(6) * t10 - t232) * t152) * qJD(4) + (-g(1) * t156 - t21 * t184 + t22 * t183 + t152 * t5 + t155 * t6 + (-g(3) * t136 - qJD(1) * t102 + t210) * t153) * m(5)) * pkin(3) - t55 * mrSges(4,2) - m(5) * (t21 * t33 + t22 * t34) + Ifges(4,3) * qJDD(3) + t174 * t101 + (m(6) * t2 + t27) * (pkin(3) * t155 + pkin(4)); t169 * g(1) + (-pkin(4) * t62 + t174) * t101 + t175 * g(3) + (t179 + (m(6) * pkin(4) + t207) * t140) * t210 + t160 - t13 * t78 - t21 * t79 + t14 * t80 + t22 * t81 + (-(-t10 + t13) * t14 + (-g(1) * t141 - g(3) * t192 - t101 * t64 + t2) * pkin(4)) * m(6) + pkin(4) * t27; -t100 * t78 + t101 * t80 + (g(2) * t136 + g(3) * t135 + t10 * t101 - t14 * t100 + t16) * m(6) + t180;];
tau = t1;
