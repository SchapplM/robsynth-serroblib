% Calculate vector of inverse dynamics joint torques for
% S5PRPRP1
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:19
% EndTime: 2019-12-05 15:28:27
% DurationCPUTime: 3.61s
% Computational Cost: add. (1310->262), mult. (2886->321), div. (0->0), fcn. (1894->8), ass. (0->107)
t161 = mrSges(5,1) + mrSges(6,1);
t157 = Ifges(5,1) + Ifges(6,1);
t155 = Ifges(6,4) + Ifges(5,5);
t63 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t156 = -Ifges(5,4) + Ifges(6,5);
t82 = sin(pkin(8));
t83 = cos(pkin(8));
t144 = mrSges(4,3) * (t82 ^ 2 + t83 ^ 2);
t81 = pkin(7) + qJ(2);
t75 = sin(t81);
t132 = g(2) * t75;
t77 = cos(t81);
t160 = g(1) * t77 + t132;
t131 = cos(qJ(4));
t85 = sin(qJ(4));
t53 = t131 * t82 + t85 * t83;
t135 = t53 / 0.2e1;
t159 = m(5) + m(6);
t68 = pkin(3) * t83 + pkin(2);
t158 = m(5) * t68;
t154 = -Ifges(5,6) + Ifges(6,6);
t106 = t131 * t83;
t90 = -t85 * t82 + t106;
t49 = t90 * qJD(4);
t28 = qJD(2) * t49 + qJDD(2) * t53;
t15 = -qJDD(4) * mrSges(6,1) + t28 * mrSges(6,2);
t153 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t28 + t15;
t113 = qJDD(2) * t82;
t50 = t53 * qJD(4);
t29 = qJD(2) * t50 - qJDD(2) * t106 + t113 * t85;
t17 = -mrSges(6,2) * t29 + qJDD(4) * mrSges(6,3);
t152 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t29 + t17;
t47 = t90 * qJD(2);
t129 = Ifges(6,5) * t47;
t44 = Ifges(5,4) * t47;
t48 = t53 * qJD(2);
t151 = t155 * qJD(4) + t157 * t48 - t129 + t44;
t125 = t47 * mrSges(5,3);
t126 = t47 * mrSges(6,2);
t40 = qJD(4) * mrSges(6,3) + t126;
t117 = -qJD(4) * mrSges(5,2) + t125 + t40;
t123 = t48 * mrSges(5,3);
t124 = t48 * mrSges(6,2);
t116 = t161 * qJD(4) - t123 - t124;
t148 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3);
t71 = t83 * qJDD(1);
t45 = -t63 * t82 + t71;
t46 = t82 * qJDD(1) + t83 * t63;
t147 = -t45 * t82 + t46 * t83;
t80 = pkin(8) + qJ(4);
t74 = sin(t80);
t76 = cos(t80);
t97 = mrSges(5,1) * t76 - mrSges(5,2) * t74;
t98 = -mrSges(4,1) * t83 + mrSges(4,2) * t82;
t145 = -m(4) * pkin(2) - mrSges(3,1) - t97 + t98;
t118 = pkin(6) + qJ(3);
t60 = t118 * t82;
t73 = t83 * qJD(1);
t41 = -qJD(2) * t60 + t73;
t114 = qJ(3) * qJD(2);
t55 = t82 * qJD(1) + t83 * t114;
t42 = pkin(6) * qJD(2) * t83 + t55;
t13 = t131 * t42 + t85 * t41;
t34 = t71 + (-pkin(6) * qJDD(2) - t63) * t82;
t112 = qJDD(2) * t83;
t35 = pkin(6) * t112 + t46;
t5 = -qJD(4) * t13 + t131 * t34 - t85 * t35;
t9 = qJD(4) * qJ(5) + t13;
t143 = -m(6) * t9 - t117;
t107 = t131 * t41;
t122 = t85 * t42;
t12 = t107 - t122;
t7 = -qJD(4) * pkin(4) + qJD(5) - t12;
t142 = -m(6) * t7 + t116;
t140 = t47 / 0.2e1;
t139 = -t47 / 0.2e1;
t137 = t48 / 0.2e1;
t134 = m(2) + m(3);
t130 = Ifges(5,4) * t48;
t111 = m(4) + t159;
t108 = qJD(4) * t107 + t131 * t35 + t85 * t34;
t102 = t29 * mrSges(5,1) + t28 * mrSges(5,2);
t101 = t29 * mrSges(6,1) - t28 * mrSges(6,3);
t99 = -mrSges(4,1) * t112 + mrSges(4,2) * t113;
t95 = t76 * mrSges(6,1) + t74 * mrSges(6,3);
t94 = pkin(4) * t76 + qJ(5) * t74;
t93 = -(-t114 * t82 + t73) * t82 + t55 * t83;
t61 = t118 * t83;
t91 = -t131 * t60 - t85 * t61;
t31 = t131 * t61 - t85 * t60;
t59 = -qJD(2) * t68 + qJD(3);
t57 = -qJDD(2) * t68 + qJDD(3);
t88 = m(6) * t94 + t95;
t69 = -qJDD(2) * pkin(2) + qJDD(3);
t43 = Ifges(6,5) * t48;
t27 = -pkin(4) * t90 - qJ(5) * t53 - t68;
t23 = -mrSges(6,1) * t47 - mrSges(6,3) * t48;
t22 = pkin(4) * t48 - qJ(5) * t47;
t19 = t47 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t130;
t18 = Ifges(6,6) * qJD(4) - t47 * Ifges(6,3) + t43;
t8 = -pkin(4) * t47 - qJ(5) * t48 + t59;
t6 = pkin(4) * t50 - qJ(5) * t49 - qJD(5) * t53;
t4 = -qJD(4) * t122 + t108;
t3 = pkin(4) * t29 - qJ(5) * t28 - qJD(5) * t48 + t57;
t2 = -qJDD(4) * pkin(4) + qJDD(5) - t5;
t1 = qJDD(4) * qJ(5) + (qJD(5) - t122) * qJD(4) + t108;
t10 = [t152 * t53 - t153 * t90 - t116 * t50 + t117 * t49 + t134 * qJDD(1) + m(4) * (t45 * t83 + t46 * t82) + m(5) * (-t12 * t50 + t13 * t49 + t4 * t53 + t5 * t90) + m(6) * (t1 * t53 - t2 * t90 + t49 * t9 + t50 * t7) + (-t111 - t134) * g(3); t27 * t101 - t68 * t102 - pkin(2) * t99 + t69 * t98 + (t154 * t50 + t155 * t49) * qJD(4) / 0.2e1 + (t155 * qJDD(4) + t157 * t28) * t135 + (t156 * t50 + t157 * t49) * t137 + t151 * t49 / 0.2e1 + (m(5) * t4 + m(6) * t1 + t152) * t31 - (-m(5) * t5 + m(6) * t2 + t153) * t91 + (t59 * mrSges(5,2) - t8 * mrSges(6,3) + Ifges(5,4) * t140 + Ifges(6,5) * t139) * t49 + (t59 * mrSges(5,1) + t8 * mrSges(6,1) + t18 / 0.2e1 - t19 / 0.2e1 + Ifges(6,3) * t139 - Ifges(5,2) * t140) * t50 + m(4) * (-pkin(2) * t69 + qJ(3) * t147 + t93 * qJD(3)) + t147 * mrSges(4,3) + (m(5) * t13 - t143) * (qJD(3) * t90 + qJD(4) * t91) + (-m(5) * t12 - t142) * (qJD(3) * t53 + qJD(4) * t31) + (-(Ifges(6,3) + Ifges(5,2)) * t90 + 0.2e1 * t156 * t135) * t29 + m(6) * (t27 * t3 + t6 * t8) + t63 * t144 + t6 * t23 + (Ifges(4,4) * t82 + Ifges(4,2) * t83) * t112 + (Ifges(4,1) * t82 + Ifges(4,4) * t83) * t113 + Ifges(3,3) * qJDD(2) + (-mrSges(5,1) * t90 + mrSges(5,2) * t53 - t158) * t57 + (-t12 * t49 - t13 * t50 + t4 * t90 - t5 * t53 - t132) * mrSges(5,3) + (t1 * t90 + t2 * t53 + t49 * t7 - t50 * t9 - t132) * mrSges(6,2) + (-t154 * t90 + t155 * t53) * qJDD(4) / 0.2e1 + (-t156 * t90 + t157 * t53) * t28 / 0.2e1 + t3 * (-mrSges(6,1) * t90 - mrSges(6,3) * t53) - t90 * (Ifges(6,5) * t28 + Ifges(6,6) * qJDD(4)) / 0.2e1 + t90 * (Ifges(5,4) * t28 + Ifges(5,6) * qJDD(4)) / 0.2e1 + ((-t118 * t159 - mrSges(6,2) - mrSges(5,3) + t148) * t77 + (-m(6) * (-t68 - t94) + t95 + t158 - t145) * t75) * g(1) + (t148 * t75 - t159 * (t118 * t75 + t77 * t68) + (-t88 + t145) * t77) * g(2); t116 * t48 - t117 * t47 - qJD(2) ^ 2 * t144 + t99 + t101 + t102 + (-g(1) * t75 + g(2) * t77) * t111 + (-t47 * t9 - t48 * t7 + t3) * m(6) + (t12 * t48 - t13 * t47 + t57) * m(5) + (-qJD(2) * t93 + t69) * m(4); (-pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t9 - t22 * t8) * m(6) + (-Ifges(5,2) * t48 + t151 + t44) * t139 + t154 * t29 + t155 * t28 + t160 * ((-m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3)) * t76 + (m(6) * pkin(4) + t161) * t74) + (t123 + t142) * t13 + (Ifges(6,2) + Ifges(5,3)) * qJDD(4) + (-t97 - t88) * g(3) + qJD(5) * t40 - t22 * t23 - pkin(4) * t15 + qJ(5) * t17 + t1 * mrSges(6,3) - t2 * mrSges(6,1) - t4 * mrSges(5,2) + t5 * mrSges(5,1) + t9 * t124 - t7 * t126 + t19 * t137 - (t154 * t48 + t155 * t47) * qJD(4) / 0.2e1 - (t157 * t47 - t130 + t18 + t43) * t48 / 0.2e1 - t59 * (mrSges(5,1) * t48 + mrSges(5,2) * t47) - t8 * (mrSges(6,1) * t48 - mrSges(6,3) * t47) + (t125 + t143) * t12 + (Ifges(6,3) * t48 + t129) * t140; -qJD(4) * t40 + t48 * t23 + (g(3) * t76 - t9 * qJD(4) - t160 * t74 + t8 * t48 + t2) * m(6) + t15;];
tau = t10;
