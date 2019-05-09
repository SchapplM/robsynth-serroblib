% Calculate vector of cutting forces with Newton-Euler
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 02:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:44:07
% EndTime: 2019-05-05 02:44:12
% DurationCPUTime: 2.01s
% Computational Cost: add. (23251->178), mult. (50255->230), div. (0->0), fcn. (34284->12), ass. (0->96)
t97 = cos(qJ(3));
t120 = qJD(2) * t97;
t94 = sin(qJ(3));
t121 = qJD(2) * t94;
t122 = cos(pkin(11));
t88 = sin(pkin(11));
t68 = -t122 * t120 + t88 * t121;
t69 = (t122 * t94 + t88 * t97) * qJD(2);
t44 = pkin(4) * t68 - qJ(5) * t69;
t141 = (2 * qJD(4)) + t44;
t100 = qJD(2) ^ 2;
t89 = sin(pkin(10));
t91 = cos(pkin(10));
t78 = g(1) * t89 - g(2) * t91;
t92 = cos(pkin(6));
t130 = t78 * t92;
t79 = -g(1) * t91 - g(2) * t89;
t87 = -g(3) + qJDD(1);
t90 = sin(pkin(6));
t95 = sin(qJ(2));
t98 = cos(qJ(2));
t139 = (t87 * t90 + t130) * t98 - t95 * t79;
t104 = -qJDD(2) * pkin(2) - t139;
t118 = qJD(2) * qJD(3);
t77 = qJDD(2) * t97 - t94 * t118;
t80 = qJD(3) * pkin(3) - qJ(4) * t121;
t86 = t97 ^ 2;
t102 = -pkin(3) * t77 + qJDD(4) + t80 * t121 + (-qJ(4) * t86 - pkin(8)) * t100 + t104;
t119 = qJD(3) * t68;
t135 = -2 * qJD(5);
t116 = t97 * t118;
t76 = qJDD(2) * t94 + t116;
t51 = t122 * t76 + t88 * t77;
t101 = (-t51 + t119) * qJ(5) + t102 + (pkin(4) * qJD(3) + t135) * t69;
t129 = t90 * t95;
t117 = t87 * t129 + t95 * t130 + t98 * t79;
t41 = -pkin(2) * t100 + qJDD(2) * pkin(8) + t117;
t61 = -t78 * t90 + t87 * t92;
t114 = -t94 * t41 + t97 * t61;
t27 = (-t76 + t116) * qJ(4) + (t100 * t94 * t97 + qJDD(3)) * pkin(3) + t114;
t125 = t97 * t41 + t94 * t61;
t28 = -pkin(3) * t100 * t86 + qJ(4) * t77 - qJD(3) * t80 + t125;
t111 = t122 * t27 - t88 * t28;
t99 = qJD(3) ^ 2;
t20 = -qJDD(3) * pkin(4) - t99 * qJ(5) + t141 * t69 + qJDD(5) - t111;
t16 = (t68 * t69 - qJDD(3)) * pkin(9) + (t51 + t119) * pkin(5) + t20;
t50 = -t122 * t77 + t76 * t88;
t60 = pkin(5) * t69 - qJD(3) * pkin(9);
t67 = t68 ^ 2;
t21 = -pkin(5) * t67 - t60 * t69 + (pkin(4) + pkin(9)) * t50 + t101;
t93 = sin(qJ(6));
t96 = cos(qJ(6));
t52 = -qJD(3) * t93 + t68 * t96;
t34 = t52 * qJD(6) + qJDD(3) * t96 + t50 * t93;
t53 = qJD(3) * t96 + t68 * t93;
t35 = -mrSges(7,1) * t52 + mrSges(7,2) * t53;
t66 = qJD(6) + t69;
t38 = -mrSges(7,2) * t66 + mrSges(7,3) * t52;
t49 = qJDD(6) + t51;
t14 = m(7) * (t16 * t96 - t21 * t93) - t34 * mrSges(7,3) + t49 * mrSges(7,1) - t53 * t35 + t66 * t38;
t33 = -t53 * qJD(6) - qJDD(3) * t93 + t50 * t96;
t39 = mrSges(7,1) * t66 - mrSges(7,3) * t53;
t15 = m(7) * (t16 * t93 + t21 * t96) + t33 * mrSges(7,3) - t49 * mrSges(7,2) + t52 * t35 - t66 * t39;
t59 = mrSges(6,1) * t69 + qJD(3) * mrSges(6,2);
t110 = t93 * t14 - t96 * t15 - m(6) * (t50 * pkin(4) + t101) + t69 * t59 + t51 * mrSges(6,3);
t58 = mrSges(6,1) * t68 - qJD(3) * mrSges(6,3);
t123 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t68 - t58;
t128 = mrSges(5,1) - mrSges(6,2);
t57 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t69;
t103 = m(5) * t102 + t51 * mrSges(5,2) + t123 * t68 + t128 * t50 + t69 * t57 - t110;
t81 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t121;
t82 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t120;
t140 = t103 + (t94 * t81 - t97 * t82) * qJD(2) - t77 * mrSges(4,1) + t76 * mrSges(4,2) + m(4) * (-pkin(8) * t100 + t104);
t137 = -2 * qJD(4);
t11 = m(3) * t139 + qJDD(2) * mrSges(3,1) - t100 * mrSges(3,2) - t140;
t133 = t11 * t98;
t127 = -mrSges(5,3) - mrSges(6,1);
t126 = t122 * t28 + t88 * t27;
t46 = -mrSges(6,2) * t68 - mrSges(6,3) * t69;
t124 = -mrSges(5,1) * t68 - mrSges(5,2) * t69 - t46;
t108 = pkin(4) * t99 - qJDD(3) * qJ(5) - t126;
t64 = t68 * t137;
t107 = -t33 * mrSges(7,1) - t52 * t38 + m(7) * (-t50 * pkin(5) - pkin(9) * t67 - t44 * t68 + t64 + ((2 * qJD(5)) + t60) * qJD(3) - t108) + t53 * t39 + t34 * mrSges(7,2);
t105 = -m(6) * (qJD(3) * t135 + t141 * t68 + t108) + t107;
t12 = m(5) * (t64 + t126) + t124 * t68 + t127 * t50 + (-mrSges(5,2) + mrSges(6,3)) * qJDD(3) + (-t57 + t59) * qJD(3) + t105;
t75 = (-mrSges(4,1) * t97 + mrSges(4,2) * t94) * qJD(2);
t109 = -m(6) * t20 - t96 * t14 - t93 * t15;
t9 = m(5) * t111 + (m(5) * t137 + t124) * t69 + t127 * t51 + t128 * qJDD(3) + t123 * qJD(3) + t109;
t7 = m(4) * t114 + qJDD(3) * mrSges(4,1) - t76 * mrSges(4,3) + qJD(3) * t82 + t88 * t12 - t75 * t121 + t122 * t9;
t8 = m(4) * t125 - qJDD(3) * mrSges(4,2) + t77 * mrSges(4,3) - qJD(3) * t81 + t122 * t12 + t75 * t120 - t88 * t9;
t4 = m(3) * t117 - t100 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t94 * t7 + t97 * t8;
t6 = m(3) * t61 + t7 * t97 + t8 * t94;
t115 = m(2) * t87 + t4 * t129 + t90 * t133 + t92 * t6;
t2 = m(2) * t79 - t11 * t95 + t4 * t98;
t1 = m(2) * t78 - t6 * t90 + (t4 * t95 + t133) * t92;
t3 = [-m(1) * g(1) - t1 * t89 + t2 * t91, t2, t4, t8, t12, -t50 * mrSges(6,2) - t68 * t58 - t110, t15; -m(1) * g(2) + t1 * t91 + t2 * t89, t1, t11, t7, t9, t50 * mrSges(6,1) - qJDD(3) * mrSges(6,3) - qJD(3) * t59 + t68 * t46 - t105, t14; -m(1) * g(3) + t115, t115, t6, t140, t103, t51 * mrSges(6,1) + qJDD(3) * mrSges(6,2) + qJD(3) * t58 + t69 * t46 - t109, t107;];
f_new  = t3;
