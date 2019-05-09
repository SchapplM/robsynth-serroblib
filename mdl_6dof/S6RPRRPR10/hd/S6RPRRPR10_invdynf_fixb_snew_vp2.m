% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPR10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-05-05 23:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR10_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR10_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR10_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:53:05
% EndTime: 2019-05-05 23:53:09
% DurationCPUTime: 1.18s
% Computational Cost: add. (13078->178), mult. (24971->212), div. (0->0), fcn. (15003->8), ass. (0->89)
t90 = sin(qJ(1));
t93 = cos(qJ(1));
t109 = -t93 * g(1) - t90 * g(2);
t136 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t109;
t92 = cos(qJ(3));
t119 = qJD(1) * t92;
t130 = cos(qJ(4));
t88 = sin(qJ(4));
t69 = -t130 * qJD(3) + t88 * t119;
t89 = sin(qJ(3));
t120 = qJD(1) * t89;
t80 = qJD(4) + t120;
t129 = t69 * t80;
t118 = qJD(1) * qJD(3);
t113 = t89 * t118;
t74 = qJDD(1) * t92 - t113;
t41 = -t69 * qJD(4) + t88 * qJDD(3) + t130 * t74;
t115 = t90 * g(1) - t93 * g(2);
t95 = qJD(1) ^ 2;
t102 = -t95 * qJ(2) + qJDD(2) - t115;
t131 = -pkin(1) - pkin(7);
t57 = t131 * qJDD(1) + t102;
t121 = t89 * g(3) + t92 * t57;
t72 = (pkin(3) * t89 - pkin(8) * t92) * qJD(1);
t94 = qJD(3) ^ 2;
t99 = qJDD(3) * pkin(3) + t94 * pkin(8) - t72 * t119 + t121;
t135 = (-t41 + t129) * qJ(5) - t99;
t133 = 2 * qJD(5);
t70 = t88 * qJD(3) + t130 * t119;
t40 = qJD(4) * t70 - t130 * qJDD(3) + t74 * t88;
t87 = sin(qJ(6));
t91 = cos(qJ(6));
t43 = t69 * t87 + t70 * t91;
t23 = -t43 * qJD(6) + t40 * t91 - t41 * t87;
t42 = t69 * t91 - t70 * t87;
t24 = t42 * qJD(6) + t40 * t87 + t41 * t91;
t78 = qJD(6) - t80;
t34 = -mrSges(7,2) * t78 + t42 * mrSges(7,3);
t35 = mrSges(7,1) * t78 - t43 * mrSges(7,3);
t55 = -pkin(5) * t80 - pkin(9) * t70;
t67 = t69 ^ 2;
t110 = m(7) * (-t67 * pkin(9) + (-pkin(4) - pkin(5)) * t40 + (-pkin(4) * t80 + t133 + t55) * t70 - t135) + t24 * mrSges(7,2) - t23 * mrSges(7,1) + t43 * t35 - t42 * t34;
t51 = -t69 * mrSges(6,2) + mrSges(6,3) * t80;
t101 = m(6) * (-0.2e1 * qJD(5) * t70 + (t70 * t80 + t40) * pkin(4) + t135) + t40 * mrSges(6,1) + t69 * t51 - t110;
t48 = -mrSges(5,2) * t80 - t69 * mrSges(5,3);
t49 = mrSges(5,1) * t80 - mrSges(5,3) * t70;
t50 = -mrSges(6,1) * t80 + mrSges(6,2) * t70;
t134 = -m(5) * t99 + t40 * mrSges(5,1) + (t49 - t50) * t70 + (mrSges(5,2) - mrSges(6,3)) * t41 + t69 * t48 + t101;
t132 = -m(2) - m(3);
t128 = (mrSges(2,1) - mrSges(3,2));
t127 = -mrSges(2,2) + mrSges(3,3);
t125 = -mrSges(5,3) - mrSges(6,2);
t112 = t92 * t118;
t73 = -qJDD(1) * t89 - t112;
t97 = t131 * t95 - t136;
t29 = (-t74 + t113) * pkin(8) + (-t73 + t112) * pkin(3) + t97;
t114 = -g(3) * t92 + t89 * t57;
t33 = -pkin(3) * t94 + qJDD(3) * pkin(8) - t72 * t120 + t114;
t124 = t130 * t33 + t88 * t29;
t46 = t69 * mrSges(6,1) - mrSges(6,3) * t70;
t123 = -t69 * mrSges(5,1) - mrSges(5,2) * t70 - t46;
t45 = t69 * pkin(4) - qJ(5) * t70;
t68 = qJDD(4) - t73;
t79 = t80 ^ 2;
t104 = -pkin(4) * t79 + t68 * qJ(5) + t80 * t133 - t69 * t45 + t124;
t108 = t130 * t29 - t88 * t33;
t18 = -t68 * pkin(4) - t79 * qJ(5) + t70 * t45 + qJDD(5) - t108;
t13 = (-t41 - t129) * pkin(9) + (t69 * t70 - t68) * pkin(5) + t18;
t14 = -t67 * pkin(5) + t40 * pkin(9) + t55 * t80 + t104;
t27 = -mrSges(7,1) * t42 + mrSges(7,2) * t43;
t63 = qJDD(6) - t68;
t11 = m(7) * (t13 * t91 - t14 * t87) - t24 * mrSges(7,3) + t63 * mrSges(7,1) - t43 * t27 + t78 * t34;
t12 = m(7) * (t13 * t87 + t14 * t91) + t23 * mrSges(7,3) - t63 * mrSges(7,2) + t42 * t27 - t78 * t35;
t106 = m(6) * t104 + t68 * mrSges(6,3) - t87 * t11 + t91 * t12 + t80 * t50;
t6 = m(5) * t124 - t68 * mrSges(5,2) + t123 * t69 + t125 * t40 - t80 * t49 + t106;
t71 = (mrSges(4,1) * t89 + mrSges(4,2) * t92) * qJD(1);
t76 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t119;
t103 = -m(6) * t18 - t91 * t11 - t87 * t12;
t8 = m(5) * t108 + (t48 + t51) * t80 + t123 * t70 + (mrSges(5,1) + mrSges(6,1)) * t68 + t125 * t41 + t103;
t4 = m(4) * t114 - qJDD(3) * mrSges(4,2) + t73 * mrSges(4,3) - qJD(3) * t76 - t71 * t120 + t130 * t6 - t88 * t8;
t75 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t120;
t9 = m(4) * t121 + qJDD(3) * mrSges(4,1) - t74 * mrSges(4,3) + qJD(3) * t75 - t71 * t119 - t134;
t116 = t92 * t4 - t89 * t9;
t105 = -m(3) * (-qJDD(1) * pkin(1) + t102) - t89 * t4 - t92 * t9;
t100 = m(4) * t97 - t73 * mrSges(4,1) + t74 * mrSges(4,2) + t76 * t119 + t75 * t120 + t130 * t8 + t88 * t6;
t98 = -m(3) * (t95 * pkin(1) + t136) + t100;
t2 = m(2) * t109 + t127 * qJDD(1) - (t128 * t95) + t98;
t1 = m(2) * t115 + t128 * qJDD(1) + t127 * t95 + t105;
t3 = [-m(1) * g(1) - t1 * t90 + t2 * t93, t2, -m(3) * g(3) + t116, t4, t6, -t40 * mrSges(6,2) - t69 * t46 + t106, t12; -m(1) * g(2) + t1 * t93 + t2 * t90, t1, -(t95 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t98, t9, t8, -t41 * mrSges(6,3) - t70 * t50 + t101, t11; (-m(1) + t132) * g(3) + t116, t132 * g(3) + t116, qJDD(1) * mrSges(3,2) - t95 * mrSges(3,3) - t105, t100, t134, -t68 * mrSges(6,1) + t41 * mrSges(6,2) + t70 * t46 - t80 * t51 - t103, t110;];
f_new  = t3;
