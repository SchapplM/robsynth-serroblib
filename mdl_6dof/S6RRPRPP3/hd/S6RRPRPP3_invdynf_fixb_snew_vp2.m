% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-05-06 12:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:29:01
% EndTime: 2019-05-06 12:29:06
% DurationCPUTime: 2.00s
% Computational Cost: add. (19615->202), mult. (42076->240), div. (0->0), fcn. (28079->8), ass. (0->92)
t102 = qJD(2) ^ 2;
t98 = sin(qJ(2));
t128 = qJD(1) * t98;
t100 = cos(qJ(2));
t103 = qJD(1) ^ 2;
t101 = cos(qJ(1));
t99 = sin(qJ(1));
t115 = -g(1) * t101 - g(2) * t99;
t77 = -pkin(1) * t103 + qJDD(1) * pkin(7) + t115;
t129 = -t100 * g(3) - t98 * t77;
t83 = (-pkin(2) * t100 - qJ(3) * t98) * qJD(1);
t48 = -qJDD(2) * pkin(2) - t102 * qJ(3) + t83 * t128 + qJDD(3) - t129;
t126 = qJD(1) * qJD(2);
t118 = t100 * t126;
t85 = qJDD(1) * t98 + t118;
t95 = sin(pkin(9));
t96 = cos(pkin(9));
t67 = qJDD(2) * t96 - t85 * t95;
t127 = qJD(1) * t100;
t81 = qJD(2) * t95 + t96 * t128;
t69 = -pkin(3) * t127 - pkin(8) * t81;
t80 = qJD(2) * t96 - t95 * t128;
t79 = t80 ^ 2;
t106 = -t67 * pkin(3) - t79 * pkin(8) + t81 * t69 + t48;
t137 = cos(qJ(4));
t97 = sin(qJ(4));
t61 = -t137 * t80 + t97 * t81;
t90 = -qJD(4) + t127;
t136 = t61 * t90;
t140 = -2 * qJD(5);
t68 = qJDD(2) * t95 + t85 * t96;
t35 = -t61 * qJD(4) + t137 * t68 + t97 * t67;
t62 = t137 * t81 + t97 * t80;
t105 = (-t35 - t136) * qJ(5) + t106 + (-pkin(4) * t90 + t140) * t62;
t139 = 2 * qJD(6);
t34 = t62 * qJD(4) - t137 * t67 + t97 * t68;
t52 = t62 * pkin(5) + qJ(6) * t90;
t53 = t62 * mrSges(7,1) + mrSges(7,3) * t90;
t55 = -t61 * mrSges(7,1) - mrSges(7,2) * t90;
t60 = t61 ^ 2;
t111 = t35 * mrSges(7,2) + t62 * t53 - m(7) * ((pkin(4) + qJ(6)) * t34 + t105 - t62 * t52 - t60 * pkin(5) + t61 * t139) - t34 * mrSges(7,3) - t61 * t55;
t56 = t62 * mrSges(6,1) - mrSges(6,2) * t90;
t109 = t111 - m(6) * (t34 * pkin(4) + t105) + t35 * mrSges(6,3) + t62 * t56;
t54 = t61 * mrSges(6,1) + mrSges(6,3) * t90;
t130 = -mrSges(5,2) * t90 + t61 * mrSges(5,3) + t54;
t51 = -mrSges(5,1) * t90 - t62 * mrSges(5,3);
t144 = m(5) * t106 + t35 * mrSges(5,2) + t62 * t51 - t109 - t130 * t61 - (mrSges(6,2) - mrSges(5,1)) * t34;
t65 = mrSges(4,2) * t127 + mrSges(4,3) * t80;
t66 = -mrSges(4,1) * t127 - mrSges(4,3) * t81;
t143 = m(4) * t48 - t67 * mrSges(4,1) + t68 * mrSges(4,2) - t80 * t65 + t81 * t66 + t144;
t122 = t99 * g(1) - t101 * g(2);
t76 = -qJDD(1) * pkin(1) - t103 * pkin(7) - t122;
t92 = t98 * t126;
t86 = qJDD(1) * t100 - t92;
t45 = (-t85 - t118) * qJ(3) + (-t86 + t92) * pkin(2) + t76;
t121 = -g(3) * t98 + t100 * t77;
t49 = -pkin(2) * t102 + qJDD(2) * qJ(3) + t83 * t127 + t121;
t117 = -0.2e1 * qJD(3) * t81 + t96 * t45 - t95 * t49;
t63 = -mrSges(4,1) * t80 + mrSges(4,2) * t81;
t22 = (-t80 * t127 - t68) * pkin(8) + (t80 * t81 - t86) * pkin(3) + t117;
t123 = 0.2e1 * qJD(3) * t80 + t95 * t45 + t96 * t49;
t25 = -pkin(3) * t79 + t67 * pkin(8) + t69 * t127 + t123;
t116 = t137 * t22 - t97 * t25;
t40 = pkin(4) * t61 - qJ(5) * t62;
t82 = qJDD(4) - t86;
t89 = t90 ^ 2;
t18 = -t82 * pkin(4) - t89 * qJ(5) + t62 * t40 + qJDD(5) - t116;
t39 = -mrSges(7,2) * t62 + mrSges(7,3) * t61;
t125 = m(7) * (t90 * t139 + (t61 * t62 - t82) * qJ(6) + (t35 - t136) * pkin(5) + t18) + t62 * t39 + t35 * mrSges(7,1);
t113 = m(6) * t18 + t125;
t42 = -mrSges(6,2) * t61 - mrSges(6,3) * t62;
t131 = -mrSges(5,1) * t61 - mrSges(5,2) * t62 - t42;
t133 = -mrSges(5,3) - mrSges(6,1);
t134 = mrSges(6,2) - mrSges(7,3);
t7 = m(5) * t116 + t131 * t62 + t133 * t35 + (-t55 + t130) * t90 + (mrSges(5,1) - t134) * t82 - t113;
t132 = t137 * t25 + t97 * t22;
t110 = -t89 * pkin(4) + t82 * qJ(5) - t61 * t40 + t132;
t124 = t90 * t53 - m(7) * (-t34 * pkin(5) - t60 * qJ(6) + qJDD(6) + (t140 - t52) * t90 + t110) - t82 * mrSges(7,2);
t112 = m(6) * (0.2e1 * qJD(5) * t90 - t110) + t124;
t8 = m(5) * t132 + (t51 - t56) * t90 + (-mrSges(5,2) + mrSges(6,3)) * t82 + (-t39 + t131) * t61 + (-mrSges(7,1) + t133) * t34 - t112;
t5 = m(4) * t117 - t86 * mrSges(4,1) - t68 * mrSges(4,3) - t65 * t127 + t137 * t7 - t81 * t63 + t97 * t8;
t6 = m(4) * t123 + t86 * mrSges(4,2) + t67 * mrSges(4,3) + t66 * t127 + t137 * t8 + t80 * t63 - t97 * t7;
t87 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t128;
t88 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t127;
t141 = m(3) * t76 - t86 * mrSges(3,1) + t85 * mrSges(3,2) - (t100 * t88 - t87 * t98) * qJD(1) + t96 * t5 + t95 * t6;
t84 = (-mrSges(3,1) * t100 + mrSges(3,2) * t98) * qJD(1);
t10 = m(3) * t129 + qJDD(2) * mrSges(3,1) - t85 * mrSges(3,3) + qJD(2) * t88 - t84 * t128 - t143;
t4 = m(3) * t121 - qJDD(2) * mrSges(3,2) + t86 * mrSges(3,3) - qJD(2) * t87 + t84 * t127 - t95 * t5 + t96 * t6;
t138 = t100 * t10 + t98 * t4;
t2 = m(2) * t122 + qJDD(1) * mrSges(2,1) - t103 * mrSges(2,2) - t141;
t1 = m(2) * t115 - t103 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t98 * t10 + t100 * t4;
t3 = [-m(1) * g(1) + t1 * t101 - t2 * t99, t1, t4, t6, t8, -t34 * mrSges(6,2) - t61 * t54 - t109, -t111; -m(1) * g(2) + t1 * t99 + t101 * t2, t2, t10, t5, t7, -t82 * mrSges(6,3) + t90 * t56 + (t39 + t42) * t61 + (mrSges(6,1) + mrSges(7,1)) * t34 + t112, -t82 * mrSges(7,3) + t90 * t55 + t125; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t141, t143, t144, t35 * mrSges(6,1) + t62 * t42 + (-t54 + t55) * t90 + t134 * t82 + t113, -t34 * mrSges(7,1) - t61 * t39 - t124;];
f_new  = t3;
