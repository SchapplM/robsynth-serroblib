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
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:03:10
% EndTime: 2019-12-05 18:03:27
% DurationCPUTime: 6.32s
% Computational Cost: add. (2644->349), mult. (5640->441), div. (0->0), fcn. (3499->12), ass. (0->157)
t237 = Ifges(5,4) + Ifges(6,4);
t151 = sin(qJ(3));
t154 = cos(qJ(3));
t123 = -mrSges(4,1) * t154 + mrSges(4,2) * t151;
t147 = qJ(3) + qJ(4);
t140 = sin(t147);
t141 = cos(t147);
t204 = mrSges(5,2) + mrSges(6,2);
t205 = mrSges(5,1) + mrSges(6,1);
t167 = t140 * t204 - t205 * t141;
t247 = -t123 - t167;
t238 = Ifges(5,1) + Ifges(6,1);
t236 = Ifges(5,5) + Ifges(6,5);
t235 = Ifges(5,2) + Ifges(6,2);
t234 = Ifges(5,6) + Ifges(6,6);
t150 = sin(qJ(4));
t153 = cos(qJ(4));
t106 = -t150 * t151 + t153 * t154;
t100 = t106 * qJD(1);
t246 = t100 / 0.2e1;
t145 = qJD(3) + qJD(4);
t245 = t145 / 0.2e1;
t244 = t237 * t100;
t177 = t204 * t141;
t107 = t150 * t154 + t151 * t153;
t101 = t107 * qJD(1);
t243 = t237 * t101;
t242 = t235 * t100 + t234 * t145 + t243;
t241 = t238 * t101 + t236 * t145 + t244;
t230 = m(3) + m(6) + m(5) + m(4);
t240 = pkin(1) * t230 + mrSges(2,1);
t212 = t151 / 0.2e1;
t148 = sin(pkin(8));
t126 = pkin(1) * t148 + pkin(6);
t119 = t126 * qJD(1);
t174 = pkin(7) * qJD(1) + t119;
t187 = qJD(2) * t151;
t75 = t154 * t174 + t187;
t69 = t150 * t75;
t139 = t154 * qJD(2);
t74 = -t151 * t174 + t139;
t72 = qJD(3) * pkin(3) + t74;
t21 = t153 * t72 - t69;
t93 = t101 * qJ(5);
t13 = t21 - t93;
t200 = mrSges(6,3) * t100;
t78 = -mrSges(6,2) * t145 + t200;
t201 = mrSges(5,3) * t100;
t79 = -mrSges(5,2) * t145 + t201;
t233 = t78 + t79;
t80 = mrSges(6,1) * t145 - t101 * mrSges(6,3);
t81 = mrSges(5,1) * t145 - t101 * mrSges(5,3);
t232 = t80 + t81;
t149 = cos(pkin(8));
t127 = -pkin(1) * t149 - pkin(2);
t142 = t154 * pkin(3);
t116 = t127 - t142;
t117 = t126 * qJDD(1);
t186 = qJD(3) * t151;
t55 = qJD(3) * t139 + t151 * qJDD(2) + t154 * t117 - t119 * t186;
t87 = t119 * t154 + t187;
t56 = -qJD(3) * t87 + t154 * qJDD(2) - t117 * t151;
t231 = -t151 * t56 + t154 * t55;
t192 = qJ(5) * t100;
t71 = t153 * t75;
t22 = t150 * t72 + t71;
t14 = t22 + t192;
t172 = t22 * mrSges(5,3) + t14 * mrSges(6,3);
t188 = pkin(4) * t141 + t142;
t229 = mrSges(3,1) + m(6) * (pkin(2) + t188) + m(5) * (t142 + pkin(2)) + m(4) * pkin(2) + t247;
t171 = mrSges(4,1) * t151 + mrSges(4,2) * t154;
t227 = t127 * qJD(1) * t171 + qJD(3) * (Ifges(4,5) * t154 - Ifges(4,6) * t151) / 0.2e1;
t156 = -pkin(7) - pkin(6);
t226 = m(4) * pkin(6) - m(5) * t156 - m(6) * (-qJ(5) + t156) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t218 = t101 / 0.2e1;
t146 = qJ(1) + pkin(8);
t136 = cos(t146);
t208 = g(3) * t136;
t203 = pkin(7) + t126;
t143 = qJDD(3) + qJDD(4);
t180 = qJD(1) * qJD(3);
t109 = qJDD(1) * t154 - t151 * t180;
t110 = qJDD(1) * t151 + t154 * t180;
t164 = t107 * qJD(4);
t48 = -qJD(1) * t164 + t109 * t153 - t110 * t150;
t29 = -mrSges(6,2) * t143 + mrSges(6,3) * t48;
t30 = -mrSges(5,2) * t143 + mrSges(5,3) * t48;
t202 = t29 + t30;
t34 = t153 * t74 - t69;
t103 = t203 * t151;
t104 = t203 * t154;
t59 = -t150 * t103 + t153 * t104;
t199 = Ifges(4,4) * t151;
t198 = Ifges(4,4) * t154;
t194 = t154 * Ifges(4,2);
t135 = sin(t146);
t190 = t135 * t140;
t185 = qJD(3) * t154;
t184 = qJD(4) * t150;
t183 = qJD(4) * t153;
t182 = t151 * qJD(1);
t181 = t154 * qJD(1);
t179 = pkin(3) * t186;
t163 = t106 * qJD(4);
t47 = qJD(1) * t163 + t109 * t150 + t110 * t153;
t178 = -t48 * mrSges(6,1) + t47 * mrSges(6,2);
t33 = -t150 * t74 - t71;
t176 = qJD(3) * t203;
t58 = -t153 * t103 - t104 * t150;
t173 = -t135 * t177 - t205 * t190;
t118 = t127 * qJDD(1);
t170 = t194 + t199;
t86 = -t119 * t151 + t139;
t168 = t87 * t151 + t86 * t154;
t32 = qJDD(3) * pkin(3) - pkin(7) * t110 + t56;
t35 = pkin(7) * t109 + t55;
t5 = t150 * t32 + t153 * t35 + t72 * t183 - t184 * t75;
t94 = t151 * t176;
t95 = t154 * t176;
t11 = -t103 * t183 - t104 * t184 - t150 * t95 - t153 * t94;
t165 = t151 * (Ifges(4,1) * t154 - t199);
t102 = t116 * qJD(1);
t160 = m(6) * (-pkin(3) * t151 - pkin(4) * t140) - t171;
t76 = -pkin(3) * t109 + t118;
t6 = -qJD(4) * t22 - t150 * t35 + t153 * t32;
t12 = -qJD(4) * t59 + t150 * t94 - t153 * t95;
t10 = pkin(4) * t145 + t13;
t2 = pkin(4) * t143 - qJ(5) * t47 - qJD(5) * t101 + t6;
t3 = qJ(5) * t48 + qJD(5) * t100 + t5;
t64 = -pkin(4) * t100 + qJD(5) + t102;
t158 = t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2) + t10 * t200 - t102 * (mrSges(5,1) * t101 + mrSges(5,2) * t100) - t64 * (mrSges(6,1) * t101 + mrSges(6,2) * t100) + t21 * t201 + t234 * t48 + t236 * t47 - (t238 * t100 - t243) * t101 / 0.2e1 + t242 * t218 - (t236 * t100 - t234 * t101) * t145 / 0.2e1 + (Ifges(6,3) + Ifges(5,3)) * t143 - (-t235 * t101 + t241 + t244) * t100 / 0.2e1;
t155 = cos(qJ(1));
t152 = sin(qJ(1));
t132 = Ifges(4,4) * t181;
t122 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t181;
t120 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t182;
t99 = Ifges(4,1) * t182 + Ifges(4,5) * qJD(3) + t132;
t98 = Ifges(4,6) * qJD(3) + qJD(1) * t170;
t92 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t110;
t91 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t109;
t77 = pkin(3) * t182 + pkin(4) * t101;
t73 = -pkin(4) * t106 + t116;
t68 = -qJD(3) * t107 - t164;
t67 = qJD(3) * t106 + t163;
t63 = -mrSges(5,1) * t100 + mrSges(5,2) * t101;
t62 = -mrSges(6,1) * t100 + mrSges(6,2) * t101;
t57 = -pkin(4) * t68 + t179;
t41 = qJ(5) * t106 + t59;
t40 = -qJ(5) * t107 + t58;
t28 = mrSges(5,1) * t143 - mrSges(5,3) * t47;
t27 = mrSges(6,1) * t143 - mrSges(6,3) * t47;
t18 = -t93 + t34;
t17 = t33 - t192;
t16 = -pkin(4) * t48 + qJDD(5) + t76;
t8 = -qJ(5) * t67 - qJD(5) * t107 + t12;
t7 = qJ(5) * t68 + qJD(5) * t106 + t11;
t1 = [(mrSges(5,2) * t76 + mrSges(6,2) * t16 - mrSges(5,3) * t6 - mrSges(6,3) * t2 + t143 * t236 + t237 * t48 + t238 * t47) * t107 + (-mrSges(5,1) * t76 - mrSges(6,1) * t16 + mrSges(5,3) * t5 + mrSges(6,3) * t3 + t143 * t234 + t235 * t48 + t237 * t47) * t106 + (t102 * mrSges(5,2) + t64 * mrSges(6,2) + t241 / 0.2e1 + t236 * t245 + t237 * t246 + t238 * t218 - t21 * mrSges(5,3) - t10 * mrSges(6,3)) * t67 + (t242 / 0.2e1 + t234 * t245 + t235 * t246 + t237 * t218 - t102 * mrSges(5,1) - t64 * mrSges(6,1) + t172) * t68 + (t154 * (-Ifges(4,2) * t151 + t198) + t165) * t180 / 0.2e1 + (m(4) * t127 + t123) * t118 + (mrSges(2,2) * t155 + t135 * t229 - t136 * t226 + t240 * t152) * g(3) + (-mrSges(2,2) * t152 + t135 * t226 + t136 * t229 + t240 * t155) * g(2) + (0.2e1 * Ifges(4,5) * t212 + Ifges(4,6) * t154) * qJDD(3) + (-t185 * t86 - t186 * t87 + t231) * mrSges(4,3) + (m(4) * (-qJD(3) * t168 + t231) - t120 * t185 - t122 * t186 + t154 * t91 - t151 * t92) * t126 + t227 * qJD(3) + (Ifges(4,1) * t110 + Ifges(4,4) * t109) * t212 + t110 * (t151 * Ifges(4,1) + t198) / 0.2e1 + t99 * t185 / 0.2e1 - t98 * t186 / 0.2e1 + m(5) * (t102 * t179 + t11 * t22 + t116 * t76 + t12 * t21 + t5 * t59 + t58 * t6) + t73 * t178 + t63 * t179 + t127 * (-mrSges(4,1) * t109 + mrSges(4,2) * t110) + t116 * (-mrSges(5,1) * t48 + mrSges(5,2) * t47) + t12 * t81 + m(6) * (t10 * t8 + t14 * t7 + t16 * t73 + t2 * t40 + t3 * t41 + t57 * t64) + t109 * t170 / 0.2e1 + t7 * t78 + t11 * t79 + t8 * t80 + t58 * t28 + t59 * t30 + t57 * t62 + t40 * t27 + t41 * t29 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t149 - 0.2e1 * mrSges(3,2) * t148 + m(3) * (t148 ^ 2 + t149 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t154 * (Ifges(4,4) * t110 + Ifges(4,2) * t109) / 0.2e1; m(3) * qJDD(2) + t151 * t91 + t154 * t92 + t232 * t68 + t233 * t67 + t202 * t107 + (t27 + t28) * t106 + (-t151 * t120 + t154 * t122) * qJD(3) + m(4) * (t151 * t55 + t154 * t56 + (-t151 * t86 + t154 * t87) * qJD(3)) + m(5) * (t106 * t6 + t107 * t5 + t21 * t68 + t22 * t67) + m(6) * (t10 * t68 + t106 * t2 + t107 * t3 + t14 * t67) - t230 * g(1); (t98 * t212 + (t194 * t212 - t165 / 0.2e1) * qJD(1) + t168 * mrSges(4,3) - (t132 + t99) * t154 / 0.2e1 - t227) * qJD(1) + (-m(6) * t188 - t247) * g(1) + (t135 * t160 + t173) * g(2) + (t140 * t205 - t160 + t177) * t208 - m(6) * (t10 * t17 + t14 * t18 + t64 * t77) + Ifges(4,6) * t109 + Ifges(4,5) * t110 + t87 * t120 - t86 * t122 - t33 * t81 - t77 * t62 - t18 * t78 - t34 * t79 - t17 * t80 - t55 * mrSges(4,2) + t56 * mrSges(4,1) - m(5) * (t21 * t33 + t22 * t34) + t158 + (-t63 * t182 + t153 * t28 + (m(6) * t3 + t202) * t150 + ((m(6) * t14 + t233) * t153 + (-m(6) * t10 - t232) * t150) * qJD(4) + (-g(1) * t154 - t21 * t184 + t22 * t183 + t150 * t5 + t153 * t6 + (-g(2) * t135 - qJD(1) * t102 + t208) * t151) * m(5)) * pkin(3) + Ifges(4,3) * qJDD(3) + t172 * t101 + (m(6) * t2 + t27) * (pkin(3) * t153 + pkin(4)); t167 * g(1) + (-pkin(4) * t62 + t172) * t101 + t173 * g(2) + (t177 + (m(6) * pkin(4) + t205) * t140) * t208 + t22 * t81 - t13 * t78 - t21 * t79 + t14 * t80 + pkin(4) * t27 + t158 + (-(-t10 + t13) * t14 + (-g(1) * t141 - g(2) * t190 - t101 * t64 + t2) * pkin(4)) * m(6); -t100 * t78 + t101 * t80 + (-g(2) * t136 - g(3) * t135 + t10 * t101 - t14 * t100 + t16) * m(6) + t178;];
tau = t1;
