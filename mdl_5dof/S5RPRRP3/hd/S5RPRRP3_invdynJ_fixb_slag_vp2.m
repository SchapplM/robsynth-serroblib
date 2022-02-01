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
% m [6x1]
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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:29:35
% EndTime: 2022-01-23 09:29:46
% DurationCPUTime: 6.24s
% Computational Cost: add. (2644->344), mult. (5640->434), div. (0->0), fcn. (3499->12), ass. (0->156)
t233 = Ifges(5,4) + Ifges(6,4);
t148 = sin(qJ(3));
t151 = cos(qJ(3));
t119 = -mrSges(4,1) * t151 + mrSges(4,2) * t148;
t144 = qJ(3) + qJ(4);
t136 = sin(t144);
t137 = cos(t144);
t199 = mrSges(5,2) + mrSges(6,2);
t200 = mrSges(5,1) + mrSges(6,1);
t165 = t136 * t199 - t200 * t137;
t243 = t119 + t165;
t234 = Ifges(5,1) + Ifges(6,1);
t232 = Ifges(5,5) + Ifges(6,5);
t231 = Ifges(5,2) + Ifges(6,2);
t230 = Ifges(5,6) + Ifges(6,6);
t147 = sin(qJ(4));
t150 = cos(qJ(4));
t106 = -t147 * t148 + t150 * t151;
t100 = t106 * qJD(1);
t242 = t100 / 0.2e1;
t142 = qJD(3) + qJD(4);
t241 = t142 / 0.2e1;
t240 = t233 * t100;
t107 = t147 * t151 + t148 * t150;
t101 = t107 * qJD(1);
t239 = t233 * t101;
t238 = t231 * t100 + t230 * t142 + t239;
t237 = t234 * t101 + t232 * t142 + t240;
t143 = qJ(1) + pkin(8);
t131 = sin(t143);
t132 = cos(t143);
t226 = g(1) * t132 + g(2) * t131;
t225 = m(3) + m(6) + m(5) + m(4);
t235 = pkin(1) * t225 + mrSges(2,1);
t208 = t148 / 0.2e1;
t145 = sin(pkin(8));
t122 = pkin(1) * t145 + pkin(6);
t115 = t122 * qJD(1);
t171 = pkin(7) * qJD(1) + t115;
t184 = qJD(2) * t148;
t75 = t151 * t171 + t184;
t69 = t147 * t75;
t135 = t151 * qJD(2);
t74 = -t148 * t171 + t135;
t72 = qJD(3) * pkin(3) + t74;
t21 = t150 * t72 - t69;
t93 = t101 * qJ(5);
t13 = t21 - t93;
t195 = mrSges(6,3) * t100;
t78 = -mrSges(6,2) * t142 + t195;
t196 = mrSges(5,3) * t100;
t79 = -mrSges(5,2) * t142 + t196;
t229 = t78 + t79;
t80 = mrSges(6,1) * t142 - mrSges(6,3) * t101;
t81 = mrSges(5,1) * t142 - t101 * mrSges(5,3);
t228 = t80 + t81;
t146 = cos(pkin(8));
t123 = -pkin(1) * t146 - pkin(2);
t138 = t151 * pkin(3);
t112 = t123 - t138;
t113 = t122 * qJDD(1);
t183 = qJD(3) * t148;
t55 = qJD(3) * t135 + t148 * qJDD(2) + t151 * t113 - t115 * t183;
t87 = t115 * t151 + t184;
t56 = -qJD(3) * t87 + t151 * qJDD(2) - t113 * t148;
t227 = -t148 * t56 + t151 * t55;
t186 = qJ(5) * t100;
t71 = t150 * t75;
t22 = t147 * t72 + t71;
t14 = t22 + t186;
t170 = t22 * mrSges(5,3) + t14 * mrSges(6,3);
t125 = pkin(4) * t137;
t185 = t125 + t138;
t224 = -m(6) * (pkin(2) + t185) - m(5) * (t138 + pkin(2)) - m(4) * pkin(2) - mrSges(3,1) + t243;
t169 = mrSges(4,1) * t148 + mrSges(4,2) * t151;
t223 = t123 * qJD(1) * t169 + qJD(3) * (Ifges(4,5) * t151 - Ifges(4,6) * t148) / 0.2e1;
t153 = -pkin(7) - pkin(6);
t222 = -m(4) * pkin(6) + m(5) * t153 + m(6) * (-qJ(5) + t153) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t214 = t101 / 0.2e1;
t205 = pkin(4) * t101;
t198 = pkin(7) + t122;
t140 = qJDD(3) + qJDD(4);
t177 = qJD(1) * qJD(3);
t109 = qJDD(1) * t151 - t148 * t177;
t110 = qJDD(1) * t148 + t151 * t177;
t162 = t107 * qJD(4);
t48 = -qJD(1) * t162 + t109 * t150 - t110 * t147;
t29 = -mrSges(6,2) * t140 + mrSges(6,3) * t48;
t30 = -mrSges(5,2) * t140 + mrSges(5,3) * t48;
t197 = t29 + t30;
t34 = t150 * t74 - t69;
t103 = t198 * t148;
t104 = t198 * t151;
t59 = -t147 * t103 + t150 * t104;
t194 = Ifges(4,4) * t148;
t193 = Ifges(4,4) * t151;
t189 = t151 * Ifges(4,2);
t182 = qJD(3) * t151;
t181 = qJD(4) * t147;
t180 = qJD(4) * t150;
t179 = t148 * qJD(1);
t178 = t151 * qJD(1);
t176 = pkin(3) * t183;
t161 = t106 * qJD(4);
t47 = qJD(1) * t161 + t109 * t147 + t110 * t150;
t175 = -t48 * mrSges(6,1) + t47 * mrSges(6,2);
t174 = t199 * t137;
t33 = -t147 * t74 - t71;
t173 = qJD(3) * t198;
t58 = -t150 * t103 - t104 * t147;
t114 = t123 * qJDD(1);
t168 = t189 + t194;
t86 = -t115 * t148 + t135;
t166 = t87 * t148 + t86 * t151;
t32 = qJDD(3) * pkin(3) - pkin(7) * t110 + t56;
t35 = pkin(7) * t109 + t55;
t5 = t147 * t32 + t150 * t35 + t72 * t180 - t181 * t75;
t94 = t148 * t173;
t95 = t151 * t173;
t11 = -t103 * t180 - t104 * t181 - t147 * t95 - t150 * t94;
t163 = t148 * (Ifges(4,1) * t151 - t194);
t102 = t112 * qJD(1);
t76 = -pkin(3) * t109 + t114;
t6 = -qJD(4) * t22 - t147 * t35 + t150 * t32;
t12 = -qJD(4) * t59 + t147 * t94 - t150 * t95;
t10 = pkin(4) * t142 + t13;
t2 = pkin(4) * t140 - qJ(5) * t47 - qJD(5) * t101 + t6;
t3 = qJ(5) * t48 + qJD(5) * t100 + t5;
t64 = -pkin(4) * t100 + qJD(5) + t102;
t155 = t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2) + t10 * t195 - t102 * (mrSges(5,1) * t101 + mrSges(5,2) * t100) - t64 * (mrSges(6,1) * t101 + mrSges(6,2) * t100) + t21 * t196 + t230 * t48 + t232 * t47 - (t234 * t100 - t239) * t101 / 0.2e1 + t238 * t214 - (t232 * t100 - t230 * t101) * t142 / 0.2e1 + (Ifges(6,3) + Ifges(5,3)) * t140 - (-t231 * t101 + t237 + t240) * t100 / 0.2e1;
t152 = cos(qJ(1));
t149 = sin(qJ(1));
t128 = Ifges(4,4) * t178;
t118 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t178;
t116 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t179;
t99 = Ifges(4,1) * t179 + Ifges(4,5) * qJD(3) + t128;
t98 = Ifges(4,6) * qJD(3) + qJD(1) * t168;
t92 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t110;
t91 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t109;
t77 = pkin(3) * t179 + t205;
t73 = -pkin(4) * t106 + t112;
t68 = -qJD(3) * t107 - t162;
t67 = qJD(3) * t106 + t161;
t63 = -mrSges(5,1) * t100 + mrSges(5,2) * t101;
t62 = -mrSges(6,1) * t100 + mrSges(6,2) * t101;
t57 = -pkin(4) * t68 + t176;
t41 = qJ(5) * t106 + t59;
t40 = -qJ(5) * t107 + t58;
t28 = mrSges(5,1) * t140 - mrSges(5,3) * t47;
t27 = mrSges(6,1) * t140 - mrSges(6,3) * t47;
t18 = -t93 + t34;
t17 = t33 - t186;
t16 = -pkin(4) * t48 + qJDD(5) + t76;
t8 = -qJ(5) * t67 - qJD(5) * t107 + t12;
t7 = qJ(5) * t68 + qJD(5) * t106 + t11;
t1 = [(0.2e1 * Ifges(4,5) * t208 + Ifges(4,6) * t151) * qJDD(3) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t146 - 0.2e1 * mrSges(3,2) * t145 + m(3) * (t145 ^ 2 + t146 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (-t182 * t86 - t183 * t87 + t227) * mrSges(4,3) + (t151 * t91 - t148 * t92 - t118 * t183 - t116 * t182 + m(4) * (-qJD(3) * t166 + t227)) * t122 + t223 * qJD(3) + (mrSges(2,2) * t152 - t131 * t224 + t132 * t222 + t235 * t149) * g(1) + (mrSges(2,2) * t149 + t131 * t222 + t132 * t224 - t235 * t152) * g(2) + (t151 * (-Ifges(4,2) * t148 + t193) + t163) * t177 / 0.2e1 + (m(4) * t123 + t119) * t114 + (Ifges(4,1) * t110 + Ifges(4,4) * t109) * t208 + (mrSges(5,2) * t76 + mrSges(6,2) * t16 - mrSges(5,3) * t6 - mrSges(6,3) * t2 + t140 * t232 + t233 * t48 + t234 * t47) * t107 + (-mrSges(5,1) * t76 - mrSges(6,1) * t16 + mrSges(5,3) * t5 + mrSges(6,3) * t3 + t140 * t230 + t231 * t48 + t233 * t47) * t106 + t110 * (t148 * Ifges(4,1) + t193) / 0.2e1 + t151 * (Ifges(4,4) * t110 + Ifges(4,2) * t109) / 0.2e1 + m(5) * (t102 * t176 + t11 * t22 + t112 * t76 + t12 * t21 + t5 * t59 + t58 * t6) + t73 * t175 + m(6) * (t10 * t8 + t14 * t7 + t16 * t73 + t2 * t40 + t3 * t41 + t57 * t64) + (t232 * t241 + t233 * t242 + t234 * t214 + t237 / 0.2e1 + t102 * mrSges(5,2) + t64 * mrSges(6,2) - t21 * mrSges(5,3) - t10 * mrSges(6,3)) * t67 + (t230 * t241 + t231 * t242 + t233 * t214 - t102 * mrSges(5,1) - t64 * mrSges(6,1) + t170 + t238 / 0.2e1) * t68 + t123 * (-mrSges(4,1) * t109 + mrSges(4,2) * t110) + t112 * (-mrSges(5,1) * t48 + mrSges(5,2) * t47) + t7 * t78 + t11 * t79 + t8 * t80 + t12 * t81 + t57 * t62 + t58 * t28 + t59 * t30 + t40 * t27 + t41 * t29 + t99 * t182 / 0.2e1 - t98 * t183 / 0.2e1 + t109 * t168 / 0.2e1 + t63 * t176; m(3) * qJDD(2) + t148 * t91 + t151 * t92 + t228 * t68 + t229 * t67 + t197 * t107 + (t27 + t28) * t106 + (-t116 * t148 + t118 * t151) * qJD(3) + m(4) * (t148 * t55 + t151 * t56 + (-t148 * t86 + t151 * t87) * qJD(3)) + m(5) * (t106 * t6 + t107 * t5 + t21 * t68 + t22 * t67) + m(6) * (t10 * t68 + t106 * t2 + t107 * t3 + t14 * t67) - t225 * g(3); t155 + (-t63 * t179 + t150 * t28 + (m(6) * t3 + t197) * t147 + ((m(6) * t14 + t229) * t150 + (-m(6) * t10 - t228) * t147) * qJD(4) + (-t21 * t181 + t22 * t180 + t147 * t5 + t150 * t6 - g(3) * t151 + (-qJD(1) * t102 + t226) * t148) * m(5)) * pkin(3) + t87 * t116 - t86 * t118 + Ifges(4,6) * t109 + Ifges(4,5) * t110 - m(5) * (t21 * t33 + t22 * t34) - t77 * t62 - t18 * t78 - t34 * t79 - t17 * t80 - t33 * t81 - t55 * mrSges(4,2) + t170 * t101 + t56 * mrSges(4,1) + (t98 * t208 + (t189 * t208 - t163 / 0.2e1) * qJD(1) + t166 * mrSges(4,3) - (t128 + t99) * t151 / 0.2e1 - t223) * qJD(1) + (-m(6) * t185 + t243) * g(3) - m(6) * (t10 * t17 + t14 * t18 + t64 * t77) + Ifges(4,3) * qJDD(3) + t226 * (-m(6) * (-pkin(3) * t148 - pkin(4) * t136) + t136 * t200 + t169 + t174) + (m(6) * t2 + t27) * (pkin(3) * t150 + pkin(4)); t165 * g(3) + (-pkin(4) * t62 + t170) * t101 + t155 + (-g(3) * t125 - t64 * t205 - (-t10 + t13) * t14 + t2 * pkin(4)) * m(6) - t13 * t78 - t21 * t79 + t14 * t80 + t22 * t81 + pkin(4) * t27 + t226 * (t174 + (m(6) * pkin(4) + t200) * t136); -t100 * t78 + t101 * t80 + (-g(1) * t131 + g(2) * t132 + t10 * t101 - t14 * t100 + t16) * m(6) + t175;];
tau = t1;
