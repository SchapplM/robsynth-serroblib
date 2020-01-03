% Calculate vector of inverse dynamics joint torques for
% S4RPRP5
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:46
% EndTime: 2019-12-31 16:44:53
% DurationCPUTime: 3.98s
% Computational Cost: add. (977->231), mult. (2424->287), div. (0->0), fcn. (1544->8), ass. (0->97)
t143 = mrSges(4,1) + mrSges(5,1);
t138 = Ifges(4,1) + Ifges(5,1);
t130 = Ifges(5,4) + Ifges(4,5);
t67 = sin(pkin(6));
t68 = cos(pkin(6));
t102 = t67 ^ 2 + t68 ^ 2;
t96 = qJD(1) * qJD(2);
t56 = qJ(2) * qJDD(1) + t96;
t137 = Ifges(4,4) - Ifges(5,5);
t70 = sin(qJ(3));
t116 = cos(qJ(3));
t93 = t116 * t68;
t77 = -t70 * t67 + t93;
t41 = t77 * qJD(1);
t110 = t41 * mrSges(4,3);
t111 = t41 * mrSges(5,2);
t33 = qJD(3) * mrSges(5,3) + t111;
t128 = -qJD(3) * mrSges(4,2) + t110 + t33;
t103 = pkin(5) + qJ(2);
t52 = t103 * t67;
t48 = qJD(1) * t52;
t53 = t103 * t68;
t49 = qJD(1) * t53;
t27 = t116 * t49 - t70 * t48;
t17 = qJD(3) * qJ(4) + t27;
t142 = -m(5) * t17 - t128;
t47 = t116 * t67 + t70 * t68;
t119 = t47 / 0.2e1;
t141 = -m(4) - m(5);
t60 = pkin(2) * t68 + pkin(1);
t140 = m(4) * t60;
t71 = sin(qJ(1));
t117 = g(2) * t71;
t129 = Ifges(5,6) - Ifges(4,6);
t113 = Ifges(5,5) * t41;
t39 = Ifges(4,4) * t41;
t42 = t47 * qJD(1);
t136 = t130 * qJD(3) + t138 * t42 - t113 + t39;
t109 = t42 * mrSges(5,2);
t127 = -mrSges(4,3) * t42 + t143 * qJD(3) - t109;
t72 = cos(qJ(1));
t135 = g(1) * t72 + t117;
t66 = pkin(6) + qJ(3);
t62 = sin(t66);
t63 = cos(t66);
t83 = mrSges(4,1) * t63 - mrSges(4,2) * t62;
t84 = -mrSges(3,1) * t68 + mrSges(3,2) * t67;
t133 = -m(3) * pkin(1) - mrSges(2,1) - t83 + t84;
t90 = m(3) * qJ(2) + mrSges(3,3);
t132 = mrSges(2,2) - mrSges(4,3) - t90;
t86 = pkin(5) * qJDD(1) + t56;
t36 = t86 * t67;
t37 = t86 * t68;
t5 = -qJD(3) * t27 - t116 * t36 - t70 * t37;
t124 = t41 / 0.2e1;
t123 = -t41 / 0.2e1;
t121 = t42 / 0.2e1;
t114 = Ifges(4,4) * t42;
t112 = t27 * mrSges(4,3);
t107 = t70 * t49;
t104 = qJD(3) / 0.2e1;
t100 = qJDD(1) * t67;
t99 = qJDD(1) * t68;
t94 = t116 * t48;
t95 = -qJD(3) * t94 + t116 * t37 - t70 * t36;
t43 = t77 * qJD(3);
t24 = qJD(1) * t43 + qJDD(1) * t47;
t44 = t47 * qJD(3);
t25 = qJD(1) * t44 - qJDD(1) * t93 + t100 * t70;
t89 = t25 * mrSges(4,1) + t24 * mrSges(4,2);
t88 = t25 * mrSges(5,1) - t24 * mrSges(5,3);
t10 = -qJDD(3) * mrSges(5,1) + t24 * mrSges(5,2);
t85 = -mrSges(3,1) * t99 + mrSges(3,2) * t100;
t81 = t63 * mrSges(5,1) + t62 * mrSges(5,3);
t80 = pkin(3) * t63 + qJ(4) * t62;
t26 = -t94 - t107;
t78 = -t116 * t52 - t70 * t53;
t29 = t116 * t53 - t70 * t52;
t51 = -qJD(1) * t60 + qJD(2);
t50 = -qJDD(1) * t60 + qJDD(2);
t75 = m(5) * t80 + t81;
t61 = -qJDD(1) * pkin(1) + qJDD(2);
t38 = Ifges(5,5) * t42;
t23 = -pkin(3) * t77 - qJ(4) * t47 - t60;
t19 = -mrSges(5,1) * t41 - mrSges(5,3) * t42;
t18 = pkin(3) * t42 - qJ(4) * t41;
t16 = -qJD(3) * pkin(3) + qJD(4) - t26;
t13 = t41 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t114;
t12 = Ifges(5,6) * qJD(3) - t41 * Ifges(5,3) + t38;
t11 = -mrSges(5,2) * t25 + qJDD(3) * mrSges(5,3);
t7 = -pkin(3) * t41 - qJ(4) * t42 + t51;
t6 = pkin(3) * t44 - qJ(4) * t43 - qJD(4) * t47;
t4 = -qJD(3) * t107 + t95;
t3 = -qJDD(3) * pkin(3) + qJDD(4) - t5;
t2 = qJDD(3) * qJ(4) + (qJD(4) - t107) * qJD(3) + t95;
t1 = pkin(3) * t25 - qJ(4) * t24 - qJD(4) * t42 + t50;
t8 = [(-t17 * t44 + t2 * t77 + t3 * t47 - t117) * mrSges(5,2) + (m(4) * t4 + m(5) * t2 - qJDD(3) * mrSges(4,2) + t11) * t29 + (-t24 * t78 + t4 * t77 - t5 * t47) * mrSges(4,3) + 0.2e1 * t102 * t56 * mrSges(3,3) + m(5) * (t1 * t23 + t6 * t7) + (m(4) * t5 - m(5) * t3 + qJDD(3) * mrSges(4,1) - t10) * t78 + t1 * (-mrSges(5,1) * t77 - mrSges(5,3) * t47) - t77 * (Ifges(5,5) * t24 + Ifges(5,6) * qJDD(3)) / 0.2e1 + t77 * (Ifges(4,4) * t24 + Ifges(4,6) * qJDD(3)) / 0.2e1 + (-t29 * mrSges(4,3) + (-Ifges(5,3) - Ifges(4,2)) * t77 - 0.2e1 * t137 * t119) * t25 + (m(4) * t27 - t142) * (qJD(2) * t77 + qJD(3) * t78) + (-t26 * mrSges(4,3) + t136 / 0.2e1 + mrSges(4,2) * t51 - mrSges(5,3) * t7 + Ifges(4,4) * t124 + Ifges(5,5) * t123 + t130 * t104 + t138 * t121 + t16 * mrSges(5,2)) * t43 + t6 * t19 + t23 * t88 - t60 * t89 + t61 * t84 - pkin(1) * t85 + (Ifges(3,4) * t67 + Ifges(3,2) * t68) * t99 + (Ifges(3,1) * t67 + Ifges(3,4) * t68) * t100 + m(3) * (-pkin(1) * t61 + (t56 + t96) * qJ(2) * t102) + (-m(4) * t26 + m(5) * t16 - t127) * (qJD(2) * t47 + qJD(3) * t29) + (-t129 * t77 + t130 * t47) * qJDD(3) / 0.2e1 + (t51 * mrSges(4,1) + t7 * mrSges(5,1) + t12 / 0.2e1 - t13 / 0.2e1 - t112 + Ifges(5,3) * t123 - Ifges(4,2) * t124 - t137 * t121 + t129 * t104) * t44 + (t130 * qJDD(3) + t138 * t24) * t119 + (t137 * t77 + t138 * t47) * t24 / 0.2e1 + (-mrSges(4,1) * t77 + mrSges(4,2) * t47 - t140) * t50 + (t141 * (t103 * t71 + t72 * t60) + t132 * t71 + (-t75 + t133) * t72) * g(2) + ((t141 * t103 - mrSges(5,2) + t132) * t72 + (-m(5) * (-t60 - t80) + t81 + t140 - t133) * t71) * g(1) + Ifges(2,3) * qJDD(1); t127 * t42 - t128 * t41 + m(3) * t61 + t85 + t88 + t89 + (-g(1) * t71 + g(2) * t72) * (m(3) - t141) - t90 * qJD(1) ^ 2 * t102 + (-t16 * t42 - t17 * t41 + t1) * m(5) + (t26 * t42 - t27 * t41 + t50) * m(4); t135 * ((-m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3)) * t63 + (m(5) * pkin(3) + t143) * t62) + (-Ifges(4,2) * t42 + t136 + t39) * t123 + t127 * t27 + t129 * t25 + t130 * t24 - (t129 * t42 + t130 * t41) * qJD(3) / 0.2e1 - t51 * (mrSges(4,1) * t42 + mrSges(4,2) * t41) - t7 * (mrSges(5,1) * t42 - mrSges(5,3) * t41) + (Ifges(5,3) * t42 + t113) * t124 + (t110 + t142) * t26 + (-pkin(3) * t3 + qJ(4) * t2 + qJD(4) * t17 - t16 * t27 - t18 * t7) * m(5) + (Ifges(5,2) + Ifges(4,3)) * qJDD(3) + qJD(4) * t33 + qJ(4) * t11 - t18 * t19 - pkin(3) * t10 + (-t83 - t75) * g(3) + t2 * mrSges(5,3) - t3 * mrSges(5,1) - t4 * mrSges(4,2) + t5 * mrSges(4,1) + t13 * t121 + t17 * t109 - t16 * t111 + t42 * t112 - (t138 * t41 - t114 + t12 + t38) * t42 / 0.2e1; -qJD(3) * t33 + t42 * t19 + (g(3) * t63 - t17 * qJD(3) - t135 * t62 + t7 * t42 + t3) * m(5) + t10;];
tau = t8;
