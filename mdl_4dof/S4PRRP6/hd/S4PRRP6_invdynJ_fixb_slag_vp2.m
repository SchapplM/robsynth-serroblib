% Calculate vector of inverse dynamics joint torques for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:13
% EndTime: 2019-12-31 16:30:18
% DurationCPUTime: 2.58s
% Computational Cost: add. (565->221), mult. (1330->293), div. (0->0), fcn. (683->6), ass. (0->96)
t145 = mrSges(4,1) + mrSges(5,1);
t150 = -mrSges(4,3) + mrSges(3,2);
t143 = Ifges(5,4) + Ifges(4,5);
t142 = Ifges(4,6) - Ifges(5,6);
t58 = sin(qJ(3));
t60 = cos(qJ(3));
t77 = t60 * mrSges(5,1) + t58 * mrSges(5,3);
t79 = mrSges(4,1) * t60 - mrSges(4,2) * t58;
t149 = -mrSges(3,1) - t77 - t79;
t100 = t58 * qJD(2);
t88 = mrSges(5,2) * t100;
t139 = -mrSges(4,3) * t100 + qJD(3) * t145 - t88;
t99 = t60 * qJD(2);
t93 = mrSges(5,2) * t99;
t39 = qJD(3) * mrSges(5,3) + t93;
t109 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t99 + t39;
t59 = sin(qJ(2));
t106 = qJD(1) * t59;
t40 = qJD(2) * pkin(5) + t106;
t113 = t40 * t60;
t20 = qJD(3) * qJ(4) + t113;
t101 = t20 * qJD(3);
t102 = qJD(3) * t60;
t114 = t40 * t58;
t15 = -qJD(3) * pkin(3) + qJD(4) + t114;
t61 = cos(qJ(2));
t98 = qJD(1) * qJD(2);
t47 = t61 * t98;
t34 = t59 * qJDD(1) + t47;
t26 = qJDD(2) * pkin(5) + t34;
t9 = t60 * t26;
t2 = qJDD(3) * qJ(4) + t9 + (qJD(4) - t114) * qJD(3);
t5 = -t102 * t40 - t26 * t58;
t3 = -qJDD(3) * pkin(3) + qJDD(4) - t5;
t82 = t2 * t60 + t3 * t58;
t148 = -t58 * t101 + t15 * t102 + t82;
t146 = -m(5) - m(4);
t123 = g(3) * t59;
t144 = mrSges(4,2) - mrSges(5,3);
t49 = Ifges(4,4) * t99;
t115 = Ifges(5,5) * t60;
t75 = Ifges(5,1) * t58 - t115;
t141 = Ifges(4,1) * t100 + qJD(2) * t75 + t143 * qJD(3) + t49;
t28 = t77 * qJD(2);
t140 = t79 * qJD(2) + t28;
t105 = qJD(1) * t61;
t41 = -qJD(2) * pkin(2) - t105;
t76 = t58 * mrSges(5,1) - t60 * mrSges(5,3);
t78 = mrSges(4,1) * t58 + mrSges(4,2) * t60;
t71 = pkin(3) * t60 + qJ(4) * t58;
t35 = -pkin(2) - t71;
t8 = qJD(2) * t35 - t105;
t138 = t41 * t78 + t8 * t76;
t137 = -t142 * t58 + t143 * t60;
t103 = qJD(3) * t58;
t4 = -t103 * t40 + t9;
t81 = t4 * t60 - t5 * t58;
t131 = t40 * (t58 ^ 2 + t60 ^ 2);
t97 = qJD(2) * qJD(3);
t31 = -t60 * qJDD(2) + t58 * t97;
t13 = -mrSges(5,2) * t31 + qJDD(3) * mrSges(5,3);
t130 = (-qJDD(3) * mrSges(4,2) - mrSges(4,3) * t31 + t13) * t60;
t32 = qJDD(2) * t58 + t60 * t97;
t12 = -qJDD(3) * mrSges(5,1) + t32 * mrSges(5,2);
t129 = (-qJDD(3) * mrSges(4,1) + mrSges(4,3) * t32 + t12) * t58;
t62 = qJD(2) ^ 2;
t127 = t58 / 0.2e1;
t118 = Ifges(4,4) * t58;
t117 = Ifges(4,4) * t60;
t116 = Ifges(5,5) * t58;
t112 = t58 * t61;
t111 = t60 * t61;
t46 = t59 * t98;
t84 = -t97 / 0.2e1;
t33 = qJDD(1) * t61 - t46;
t74 = t60 * Ifges(4,2) + t118;
t70 = pkin(3) * t58 - qJ(4) * t60;
t69 = t15 * t60 - t20 * t58;
t66 = t58 * (Ifges(4,1) * t60 - t118);
t65 = t60 * (Ifges(5,3) * t58 + t115);
t25 = -qJDD(2) * pkin(2) - t33;
t57 = cos(pkin(6));
t56 = sin(pkin(6));
t48 = Ifges(5,5) * t100;
t30 = t70 * qJD(2);
t22 = Ifges(4,6) * qJD(3) + qJD(2) * t74;
t21 = Ifges(5,6) * qJD(3) - Ifges(5,3) * t99 + t48;
t19 = t111 * t57 + t56 * t58;
t18 = t112 * t57 - t56 * t60;
t17 = t111 * t56 - t57 * t58;
t16 = t112 * t56 + t57 * t60;
t14 = qJD(3) * t70 - qJD(4) * t58;
t7 = mrSges(4,1) * t31 + mrSges(4,2) * t32;
t6 = mrSges(5,1) * t31 - mrSges(5,3) * t32;
t1 = pkin(3) * t31 - qJ(4) * t32 - qJD(4) * t100 + t25;
t10 = [m(2) * qJDD(1) + (-m(2) - m(3) + t146) * g(3) + (qJDD(2) * mrSges(3,1) - t62 * mrSges(3,2) - t6 - t7 + (t109 * t60 - t139 * t58) * qJD(2) + m(3) * t33 + m(5) * (t100 * t15 + t20 * t99 - t1) + m(4) * (qJD(2) * t131 - t25)) * t61 + (-t62 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t130 + t129 - t140 * qJD(2) + (-t109 * t58 - t139 * t60) * qJD(3) + m(3) * t34 + m(5) * (qJD(2) * t8 + t148) + m(4) * (qJD(2) * t41 + t81)) * t59; (t150 * t59 + t146 * (t61 * pkin(2) + t59 * pkin(5)) + (-m(5) * t71 + t149) * t61) * g(3) + (-t22 / 0.2e1 + t21 / 0.2e1) * t103 + (t58 * (Ifges(5,1) * t60 + t116) + t60 * (-Ifges(4,2) * t58 + t117) + t66) * t97 / 0.2e1 + (Ifges(4,1) * t58 + t117 + t75) * t32 / 0.2e1 + (g(1) * t57 + g(2) * t56) * ((m(4) * pkin(2) - m(5) * t35 - t149) * t59 + (t146 * pkin(5) - mrSges(5,2) + t150) * t61) - t25 * t79 - t1 * t77 + (t116 / 0.2e1 - t74 / 0.2e1 + (-Ifges(5,3) - Ifges(4,2) / 0.2e1) * t60 + (Ifges(5,5) - Ifges(4,4)) * t127) * t31 + (t46 + t33) * mrSges(3,1) - t60 * (Ifges(5,5) * t32 + Ifges(5,6) * qJDD(3)) / 0.2e1 + t60 * (Ifges(4,4) * t32 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t140 * t106 + t141 * t102 / 0.2e1 + ((Ifges(4,1) + Ifges(5,1)) * t32 + t143 * qJDD(3)) * t127 + (t142 * t60 + t143 * t58) * qJDD(3) / 0.2e1 + t137 * qJD(3) ^ 2 / 0.2e1 + t138 * qJD(3) + t81 * mrSges(4,3) + (-pkin(2) * t25 + pkin(5) * t81 - (t131 * t61 + t41 * t59) * qJD(1)) * m(4) + t35 * t6 + (t1 * t35 + t14 * t8 + (qJD(3) * t69 + t82) * pkin(5) - (t59 * t8 + (t15 * t58 + t20 * t60) * t61) * qJD(1)) * m(5) + (-t123 + t148) * mrSges(5,2) + t109 * (-pkin(5) * t103 - t60 * t105) - t139 * (pkin(5) * t102 - t58 * t105) + t65 * t84 + pkin(5) * t129 + pkin(5) * t130 + Ifges(3,3) * qJDD(2) + (t47 - t34) * mrSges(3,2) - pkin(2) * t7 - t14 * t28; t20 * t88 + (t78 + t76) * t123 + t22 * t100 / 0.2e1 - t15 * t93 - (Ifges(5,1) * t99 + t21 + t48) * t100 / 0.2e1 + (Ifges(5,2) + Ifges(4,3)) * qJDD(3) + (-t66 / 0.2e1 + t65 / 0.2e1) * t62 + t139 * t113 - (-Ifges(4,2) * t100 + t141 + t49) * t99 / 0.2e1 - t142 * t31 + t143 * t32 + (t144 * t19 + t145 * t18) * g(1) + (t144 * t17 + t145 * t16) * g(2) + t137 * t84 - t138 * qJD(2) + t109 * t114 + qJD(4) * t39 + (-t3 * pkin(3) + t2 * qJ(4) + t20 * qJD(4) - t8 * t30 - t40 * t69 + t70 * t123 - g(2) * (-pkin(3) * t16 + qJ(4) * t17) - g(1) * (-pkin(3) * t18 + qJ(4) * t19)) * m(5) + t2 * mrSges(5,3) - t3 * mrSges(5,1) - t4 * mrSges(4,2) + t5 * mrSges(4,1) - pkin(3) * t12 + qJ(4) * t13 + t30 * t28; -t28 * t100 - qJD(3) * t39 + (-g(1) * t18 - g(2) * t16 + t100 * t8 - t123 * t58 - t101 + t3) * m(5) + t12;];
tau = t10;
