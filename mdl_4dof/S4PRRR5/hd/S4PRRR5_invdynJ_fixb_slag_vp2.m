% Calculate vector of inverse dynamics joint torques for
% S4PRRR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:35
% EndTime: 2019-12-31 16:33:38
% DurationCPUTime: 1.33s
% Computational Cost: add. (893->176), mult. (1623->253), div. (0->0), fcn. (978->10), ass. (0->93)
t77 = cos(qJ(4));
t137 = mrSges(5,1) * t77;
t147 = m(5) * pkin(3);
t71 = qJ(2) + qJ(3);
t65 = sin(t71);
t66 = cos(t71);
t83 = mrSges(4,2) * t66 + (mrSges(4,1) + t137 + t147) * t65;
t160 = m(4) + m(5);
t74 = sin(qJ(4));
t121 = t74 * mrSges(5,2);
t158 = -t65 * t121 + t66 * (-m(5) * pkin(6) - mrSges(5,3));
t146 = t74 / 0.2e1;
t76 = sin(qJ(2));
t79 = cos(qJ(2));
t157 = mrSges(3,2) * t79 + (t160 * pkin(2) + mrSges(3,1)) * t76 + t83;
t95 = t121 - t137;
t113 = qJD(1) * t76;
t75 = sin(qJ(3));
t105 = t75 * t113;
t54 = qJD(2) * pkin(2) + qJD(1) * t79;
t78 = cos(qJ(3));
t26 = t54 * t78 - t105;
t68 = qJD(2) + qJD(3);
t20 = -pkin(3) * t68 - t26;
t94 = mrSges(5,1) * t74 + mrSges(5,2) * t77;
t155 = t20 * t94 + qJD(4) * (Ifges(5,5) * t77 - Ifges(5,6) * t74) / 0.2e1;
t111 = qJD(4) * t74;
t27 = t113 * t78 + t54 * t75;
t21 = pkin(6) * t68 + t27;
t67 = qJDD(2) + qJDD(3);
t112 = qJD(3) * t78;
t109 = qJD(1) * qJD(2);
t103 = t76 * t109;
t47 = t79 * qJDD(1) - t103;
t38 = qJDD(2) * pkin(2) + t47;
t102 = t79 * t109;
t48 = qJDD(1) * t76 + t102;
t8 = -qJD(3) * t105 + t54 * t112 + t75 * t38 + t78 * t48;
t5 = pkin(6) * t67 + t8;
t2 = -t111 * t21 + t5 * t77;
t138 = t2 * t77;
t110 = qJD(4) * t77;
t3 = -t110 * t21 - t5 * t74;
t97 = -t3 * t74 + t138;
t153 = m(5) * t97;
t151 = (-mrSges(4,1) + t95) * t66 + (mrSges(4,2) - mrSges(5,3)) * t65;
t31 = t95 * t68;
t150 = -t68 * mrSges(4,1) + t31;
t32 = -t111 * t68 + t67 * t77;
t33 = t110 * t68 + t67 * t74;
t149 = t77 * (-qJDD(4) * mrSges(5,2) + mrSges(5,3) * t32) - t74 * (qJDD(4) * mrSges(5,1) - mrSges(5,3) * t33);
t120 = t74 * mrSges(5,3);
t45 = qJD(4) * mrSges(5,1) - t120 * t68;
t126 = t68 * t77;
t46 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t126;
t148 = -t68 * mrSges(4,2) - t74 * t45 + t77 * t46;
t9 = -qJD(3) * t27 + t38 * t78 - t48 * t75;
t144 = pkin(2) * t75;
t142 = pkin(2) * t78;
t141 = pkin(2) * t79;
t135 = Ifges(5,1) * t74;
t134 = Ifges(5,4) * t74;
t133 = Ifges(5,4) * t77;
t132 = Ifges(5,2) * t77;
t129 = t67 * mrSges(4,2);
t72 = sin(pkin(7));
t125 = t72 * t74;
t124 = t72 * t77;
t73 = cos(pkin(7));
t123 = t73 * t74;
t122 = t73 * t77;
t115 = t66 * pkin(3) + t65 * pkin(6);
t114 = pkin(2) * qJD(3);
t104 = (t74 ^ 2 + t77 ^ 2) * t21;
t99 = t158 * t72;
t98 = t158 * t73;
t93 = t132 + t134;
t90 = t45 * t77 + t46 * t74;
t42 = t75 * t79 + t76 * t78;
t41 = t75 * t76 - t78 * t79;
t87 = t74 * (Ifges(5,1) * t77 - t134);
t82 = -qJD(4) * t90 + t149;
t28 = Ifges(5,6) * qJD(4) + t68 * t93;
t55 = Ifges(5,4) * t126;
t29 = Ifges(5,5) * qJD(4) + t135 * t68 + t55;
t6 = -pkin(3) * t67 - t9;
t81 = t9 * mrSges(4,1) - t8 * mrSges(4,2) + mrSges(5,3) * t138 - t3 * t120 + t6 * t95 + (Ifges(5,1) * t33 + Ifges(5,4) * t32) * t146 + t77 * (Ifges(5,4) * t33 + Ifges(5,2) * t32) / 0.2e1 + t32 * t93 / 0.2e1 + t33 * (t133 + t135) / 0.2e1 - t28 * t111 / 0.2e1 + Ifges(4,3) * t67 + (t29 + t68 * (-Ifges(5,2) * t74 + t133)) * t110 / 0.2e1 + (0.2e1 * Ifges(5,5) * t146 + Ifges(5,6) * t77) * qJDD(4) + (t87 * t68 / 0.2e1 + t155) * qJD(4);
t80 = qJD(2) ^ 2;
t62 = -pkin(3) - t142;
t15 = t68 * t42;
t14 = t68 * t41;
t13 = -mrSges(5,1) * t32 + mrSges(5,2) * t33;
t1 = [m(2) * qJDD(1) + t41 * t13 + t15 * t31 + (-qJDD(2) * t76 - t79 * t80) * mrSges(3,2) + (-t15 * t68 - t41 * t67) * mrSges(4,1) + (qJDD(2) * t79 - t76 * t80) * mrSges(3,1) - t148 * t14 + (t82 - t129) * t42 + m(3) * (t47 * t79 + t48 * t76) + m(4) * (-t14 * t27 - t15 * t26 - t41 * t9 + t42 * t8) + m(5) * (-t104 * t14 + t15 * t20 + t41 * t6 + t42 * t97) + (-m(2) - m(3) - t160) * g(3); t81 - t129 * t144 + m(5) * (t6 * t62 + (t104 * t78 + t20 * t75) * t114) + t62 * t13 + t67 * mrSges(4,1) * t142 + Ifges(3,3) * qJDD(2) + t150 * t75 * t114 + (t102 - t48) * mrSges(3,2) + (t103 + t47) * mrSges(3,1) + (-m(4) * t141 - m(5) * (t115 + t141) - t79 * mrSges(3,1) + t76 * mrSges(3,2) + t151) * g(3) + (-t110 * t45 - t111 * t46 + t149 + t153) * (pkin(6) + t144) + ((m(4) * t26 - m(5) * t20 - t150) * t42 - (-m(4) * t27 - m(5) * t104 - t148) * t41) * qJD(1) + (t157 * t72 + t99) * g(2) + (t157 * t73 + t98) * g(1) + (m(4) * (t75 * t8 + t78 * t9 + (-t26 * t75 + t27 * t78) * qJD(3)) + t148 * t112) * pkin(2); t81 + (t82 + t153) * pkin(6) - m(5) * (t104 * t26 + t20 * t27) + (t73 * t83 + t98) * g(1) + (t72 * t83 + t99) * g(2) - t150 * t27 - t148 * t26 - t6 * t147 + (-m(5) * t115 + t151) * g(3) - pkin(3) * t13; Ifges(5,5) * t33 + Ifges(5,6) * t32 + Ifges(5,3) * qJDD(4) - t2 * mrSges(5,2) + t3 * mrSges(5,1) - g(1) * ((-t123 * t66 + t124) * mrSges(5,1) + (-t122 * t66 - t125) * mrSges(5,2)) - g(2) * ((-t125 * t66 - t122) * mrSges(5,1) + (-t124 * t66 + t123) * mrSges(5,2)) + g(3) * t94 * t65 + t90 * t21 + (t28 * t146 + (-t87 / 0.2e1 + t132 * t146) * t68 - (t29 + t55) * t77 / 0.2e1 - t155) * t68;];
tau = t1;
