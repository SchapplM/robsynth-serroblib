% Calculate vector of inverse dynamics joint torques for
% S5RPPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:04
% EndTime: 2019-12-31 18:01:08
% DurationCPUTime: 2.00s
% Computational Cost: add. (2097->234), mult. (3173->305), div. (0->0), fcn. (1614->10), ass. (0->112)
t86 = sin(qJ(5));
t89 = cos(qJ(5));
t61 = -mrSges(6,1) * t89 + mrSges(6,2) * t86;
t163 = m(6) * pkin(4) + mrSges(5,1) - t61;
t162 = qJD(1) - qJD(4);
t148 = t86 / 0.2e1;
t84 = sin(pkin(8));
t122 = qJ(2) * t84;
t92 = -pkin(1) - pkin(2);
t70 = qJD(1) * t92 + qJD(2);
t85 = cos(pkin(8));
t60 = t85 * t70;
t30 = t60 + (-pkin(3) - t122) * qJD(1);
t118 = qJ(2) * qJD(1);
t34 = t118 * t85 + t70 * t84;
t87 = sin(qJ(4));
t90 = cos(qJ(4));
t13 = t30 * t90 - t34 * t87;
t101 = mrSges(6,1) * t86 + mrSges(6,2) * t89;
t9 = pkin(4) * t162 - t13;
t161 = t9 * t101;
t159 = mrSges(2,1) + mrSges(3,1);
t158 = -mrSges(3,3) + mrSges(2,2);
t46 = t84 * t87 - t90 * t85;
t157 = t162 * t46;
t48 = t84 * t90 + t85 * t87;
t124 = t162 * t48;
t58 = t85 * t92 - t122;
t50 = -pkin(3) + t58;
t59 = qJ(2) * t85 + t84 * t92;
t21 = t87 * t50 + t90 * t59;
t116 = qJD(1) * qJD(2);
t72 = qJDD(1) * qJ(2) + t116;
t130 = t86 * mrSges(6,3);
t56 = qJD(5) * mrSges(6,1) + t130 * t162;
t127 = t89 * mrSges(6,3);
t57 = -qJD(5) * mrSges(6,2) - t127 * t162;
t156 = t86 * t56 - t89 * t57;
t121 = qJD(5) * t86;
t82 = -qJDD(1) + qJDD(4);
t40 = t121 * t162 + t82 * t89;
t28 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t40;
t120 = qJD(5) * t89;
t41 = -t120 * t162 + t82 * t86;
t29 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t41;
t155 = t89 * t28 - t86 * t29;
t17 = -mrSges(6,1) * t40 + mrSges(6,2) * t41;
t14 = t30 * t87 + t34 * t90;
t66 = qJDD(1) * t92 + qJDD(2);
t26 = t85 * t66 - t72 * t84;
t24 = -qJDD(1) * pkin(3) + t26;
t27 = t66 * t84 + t72 * t85;
t6 = -qJD(4) * t14 + t24 * t90 - t27 * t87;
t4 = -pkin(4) * t82 - t6;
t154 = m(6) * t4 + t17;
t114 = pkin(8) + qJ(4);
t106 = sin(t114);
t107 = cos(t114);
t88 = sin(qJ(1));
t91 = cos(qJ(1));
t35 = -t106 * t88 - t107 * t91;
t36 = t106 * t91 - t107 * t88;
t153 = t36 * mrSges(5,2) + t163 * t35;
t152 = t35 * mrSges(5,2) - t163 * t36;
t37 = t61 * t162;
t151 = -m(6) * t9 - mrSges(5,1) * t162 + t37;
t10 = -pkin(7) * t162 + t14;
t7 = qJD(3) * t89 - t10 * t86;
t8 = qJD(3) * t86 + t10 * t89;
t104 = -t7 * t86 + t8 * t89;
t150 = -m(6) * t104 - mrSges(5,2) * t162 + t156;
t5 = t13 * qJD(4) + t87 * t24 + t90 * t27;
t3 = pkin(7) * t82 + t5;
t1 = qJD(5) * t7 + qJDD(3) * t86 + t3 * t89;
t105 = t7 * t89 + t8 * t86;
t2 = -qJD(5) * t8 + qJDD(3) * t89 - t3 * t86;
t95 = -qJD(5) * t105 + t1 * t89 - t2 * t86;
t149 = m(6) * t95 - t56 * t120 - t57 * t121 + t155;
t146 = m(4) + m(5);
t145 = mrSges(4,1) * t84;
t144 = mrSges(4,2) * t85;
t143 = Ifges(6,4) * t86;
t142 = Ifges(6,4) * t89;
t141 = Ifges(6,2) * t89;
t138 = t82 * mrSges(5,1);
t137 = t82 * mrSges(5,2);
t134 = t162 * t86;
t133 = t162 * t89;
t132 = t84 * t88;
t131 = t84 * t91;
t123 = t91 * pkin(1) + t88 * qJ(2);
t119 = qJDD(1) * pkin(1);
t117 = m(6) + t146;
t110 = -m(6) * pkin(7) - mrSges(6,3);
t80 = t91 * qJ(2);
t109 = -pkin(1) * t88 + t80;
t76 = pkin(3) * t85 + pkin(2);
t108 = pkin(3) * t132 + t91 * t76 + t123;
t100 = Ifges(6,1) * t89 - t143;
t63 = t141 + t143;
t98 = -(-t118 * t84 + t60) * t84 + t34 * t85;
t20 = t50 * t90 - t59 * t87;
t52 = qJD(5) * (Ifges(6,5) * t89 - Ifges(6,6) * t86);
t31 = Ifges(6,6) * qJD(5) - t162 * t63;
t71 = Ifges(6,4) * t133;
t32 = -Ifges(6,1) * t134 + Ifges(6,5) * qJD(5) - t71;
t94 = Ifges(5,3) * t82 + t4 * t61 + t6 * mrSges(5,1) + qJD(5) * t161 + t1 * t127 + t40 * t63 / 0.2e1 + t41 * (Ifges(6,1) * t86 + t142) / 0.2e1 - t5 * mrSges(5,2) + (Ifges(6,1) * t41 + Ifges(6,4) * t40) * t148 + t89 * (Ifges(6,4) * t41 + Ifges(6,2) * t40) / 0.2e1 - t31 * t121 / 0.2e1 + t32 * t120 / 0.2e1 - t2 * t130 + (-t120 * t7 - t121 * t8) * mrSges(6,3) + (0.2e1 * Ifges(6,5) * t148 + Ifges(6,6) * t89) * qJDD(5) + (t52 - t100 * t134 - (-Ifges(6,2) * t86 + t142) * t133) * qJD(5) / 0.2e1;
t78 = qJDD(2) - t119;
t74 = pkin(3) * t131;
t49 = t85 * t91 + t132;
t47 = -t85 * t88 + t131;
t11 = [t116 * t144 + t116 * t145 + t20 * t138 - t21 * t137 + m(3) * (-pkin(1) * t78 + (t72 + t116) * qJ(2)) + m(4) * (t98 * qJD(2) + t26 * t58 + t27 * t59) - t26 * mrSges(4,1) + t27 * mrSges(4,2) - t94 + m(5) * (t20 * t6 + t21 * t5) + t154 * (pkin(4) - t20) + (-m(5) * t13 - t151) * (t48 * qJD(2) + qJD(4) * t21) + (t119 - t78) * mrSges(3,1) + (m(5) * t14 - t150) * (-qJD(2) * t46 + qJD(4) * t20) + 0.2e1 * t72 * mrSges(3,3) + t149 * (-pkin(7) + t21) + (-t58 * mrSges(4,1) + t59 * mrSges(4,2) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1) + (-t49 * mrSges(4,1) + t47 * mrSges(4,2) - m(6) * (pkin(7) * t36 + t108) - t36 * mrSges(6,3) - m(5) * t108 + (-m(4) - m(3)) * t123 + (-m(4) * pkin(2) - t159) * t91 + t158 * t88 + t153) * g(2) + (-m(4) * t80 - t47 * mrSges(4,1) - t49 * mrSges(4,2) - m(6) * (pkin(7) * t35 + t109 + t74) - t35 * mrSges(6,3) - m(3) * t109 - m(5) * (t74 + t80) + t158 * t91 + (-m(4) * t92 + m(6) * t76 - m(5) * (-pkin(1) - t76) + t159) * t88 + t152) * g(1); (t17 - t138) * t46 + t124 * t37 + (-mrSges(4,1) * t85 + mrSges(4,2) * t84 - mrSges(3,1)) * qJDD(1) + (-m(3) * qJ(2) - mrSges(3,3) - t144 - t145) * qJD(1) ^ 2 - (t124 * mrSges(5,1) - mrSges(5,2) * t157) * t162 + (-t137 + (-t89 * t56 - t86 * t57) * qJD(5) + t155) * t48 + m(3) * t78 - t157 * t156 + (-g(1) * t88 + g(2) * t91) * (m(3) + t117) + (t104 * t157 - t124 * t9 + t4 * t46 + t95 * t48) * m(6) + (t124 * t13 + t14 * t157 - t46 * t6 + t48 * t5) * m(5) + (-qJD(1) * t98 + t26 * t85 + t27 * t84) * m(4); m(6) * (qJD(5) * t104 + t1 * t86 + t2 * t89) + t57 * t120 + t86 * t28 - t56 * t121 + t89 * t29 + t146 * qJDD(3) + t117 * g(3); t94 + (-t110 * t36 - t153) * g(2) + (-t110 * t35 - t152) * g(1) - t154 * pkin(4) + t151 * t14 + t150 * t13 + t149 * pkin(7); t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t41 + Ifges(6,6) * t40 + Ifges(6,3) * qJDD(5) - g(3) * t61 + t8 * t56 - t7 * t57 - (-t161 + t31 * t148 - t52 / 0.2e1 - (-t86 * t100 / 0.2e1 + t141 * t148) * t162 + t105 * mrSges(6,3) - (t32 - t71) * t89 / 0.2e1) * t162 + (-t35 * g(1) - t36 * g(2)) * t101;];
tau = t11;
