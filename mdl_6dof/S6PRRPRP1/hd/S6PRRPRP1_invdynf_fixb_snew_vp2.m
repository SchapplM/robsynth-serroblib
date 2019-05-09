% Calculate vector of cutting forces with Newton-Euler
% S6PRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-05-05 03:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:36:26
% EndTime: 2019-05-05 03:36:32
% DurationCPUTime: 2.49s
% Computational Cost: add. (30417->175), mult. (65168->230), div. (0->0), fcn. (45912->12), ass. (0->90)
t136 = -2 * qJD(4);
t102 = qJD(2) ^ 2;
t100 = cos(qJ(2));
t91 = sin(pkin(6));
t97 = sin(qJ(2));
t130 = t91 * t97;
t90 = sin(pkin(10));
t93 = cos(pkin(10));
t80 = t90 * g(1) - t93 * g(2);
t94 = cos(pkin(6));
t131 = t80 * t94;
t81 = -t93 * g(1) - t90 * g(2);
t88 = -g(3) + qJDD(1);
t116 = t100 * t81 + t88 * t130 + t97 * t131;
t51 = -t102 * pkin(2) + qJDD(2) * pkin(8) + t116;
t67 = -t91 * t80 + t94 * t88;
t96 = sin(qJ(3));
t99 = cos(qJ(3));
t111 = -t96 * t51 + t99 * t67;
t121 = qJD(2) * qJD(3);
t114 = t99 * t121;
t78 = t96 * qJDD(2) + t114;
t30 = (-t78 + t114) * qJ(4) + (t102 * t96 * t99 + qJDD(3)) * pkin(3) + t111;
t127 = t99 * t51 + t96 * t67;
t79 = t99 * qJDD(2) - t96 * t121;
t124 = qJD(2) * t96;
t82 = qJD(3) * pkin(3) - qJ(4) * t124;
t87 = t99 ^ 2;
t31 = -t87 * t102 * pkin(3) + t79 * qJ(4) - qJD(3) * t82 + t127;
t89 = sin(pkin(11));
t92 = cos(pkin(11));
t72 = (t89 * t99 + t92 * t96) * qJD(2);
t135 = t136 * t72 + t92 * t30 - t89 * t31;
t71 = (t89 * t96 - t92 * t99) * qJD(2);
t134 = (t88 * t91 + t131) * t100 - t97 * t81;
t105 = -qJDD(2) * pkin(2) - t134;
t103 = -t79 * pkin(3) + qJDD(4) + t82 * t124 + (-qJ(4) * t87 - pkin(8)) * t102 + t105;
t101 = qJD(3) ^ 2;
t117 = t136 * t71 + t89 * t30 + t92 * t31;
t54 = t71 * pkin(4) - t72 * pkin(9);
t23 = -t101 * pkin(4) + qJDD(3) * pkin(9) - t71 * t54 + t117;
t58 = -t89 * t78 + t92 * t79;
t59 = t92 * t78 + t89 * t79;
t26 = (qJD(3) * t71 - t59) * pkin(9) + (qJD(3) * t72 - t58) * pkin(4) + t103;
t95 = sin(qJ(5));
t98 = cos(qJ(5));
t113 = -t95 * t23 + t98 * t26;
t61 = t98 * qJD(3) - t95 * t72;
t39 = t61 * qJD(5) + t95 * qJDD(3) + t98 * t59;
t70 = qJD(5) + t71;
t45 = -t70 * mrSges(7,2) + t61 * mrSges(7,3);
t57 = qJDD(5) - t58;
t62 = t95 * qJD(3) + t98 * t72;
t120 = m(7) * (-0.2e1 * qJD(6) * t62 + (t61 * t70 - t39) * qJ(6) + (t61 * t62 + t57) * pkin(5) + t113) + t70 * t45 + t57 * mrSges(7,1);
t41 = -t61 * mrSges(7,1) + t62 * mrSges(7,2);
t42 = -t61 * mrSges(6,1) + t62 * mrSges(6,2);
t46 = -t70 * mrSges(6,2) + t61 * mrSges(6,3);
t13 = m(6) * t113 + t57 * mrSges(6,1) + t70 * t46 + (-t42 - t41) * t62 + (-mrSges(6,3) - mrSges(7,3)) * t39 + t120;
t128 = t98 * t23 + t95 * t26;
t38 = -t62 * qJD(5) + t98 * qJDD(3) - t95 * t59;
t47 = t70 * pkin(5) - t62 * qJ(6);
t60 = t61 ^ 2;
t119 = m(7) * (-t60 * pkin(5) + t38 * qJ(6) + 0.2e1 * qJD(6) * t61 - t70 * t47 + t128) + t61 * t41 + t38 * mrSges(7,3);
t48 = t70 * mrSges(7,1) - t62 * mrSges(7,3);
t49 = t70 * mrSges(6,1) - t62 * mrSges(6,3);
t16 = m(6) * t128 + t38 * mrSges(6,3) + t61 * t42 + (-t49 - t48) * t70 + (-mrSges(6,2) - mrSges(7,2)) * t57 + t119;
t65 = -qJD(3) * mrSges(5,2) - t71 * mrSges(5,3);
t66 = qJD(3) * mrSges(5,1) - t72 * mrSges(5,3);
t107 = -m(5) * t103 + t58 * mrSges(5,1) - t59 * mrSges(5,2) - t98 * t13 - t95 * t16 - t71 * t65 - t72 * t66;
t83 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t124;
t123 = qJD(2) * t99;
t84 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t123;
t133 = (t96 * t83 - t99 * t84) * qJD(2) + m(4) * (-t102 * pkin(8) + t105) - t79 * mrSges(4,1) + t78 * mrSges(4,2) - t107;
t22 = -qJDD(3) * pkin(4) - t101 * pkin(9) + t72 * t54 - t135;
t118 = m(7) * (-t38 * pkin(5) - t60 * qJ(6) + t62 * t47 + qJDD(6) + t22) + t62 * t48 + t39 * mrSges(7,2);
t132 = -m(6) * t22 - t39 * mrSges(6,2) + (t46 + t45) * t61 + (mrSges(6,1) + mrSges(7,1)) * t38 - t62 * t49 - t118;
t10 = m(3) * t134 + qJDD(2) * mrSges(3,1) - t102 * mrSges(3,2) - t133;
t125 = t10 * t100;
t53 = t71 * mrSges(5,1) + t72 * mrSges(5,2);
t11 = m(5) * t117 - qJDD(3) * mrSges(5,2) + t58 * mrSges(5,3) - qJD(3) * t66 - t95 * t13 + t98 * t16 - t71 * t53;
t14 = m(5) * t135 + qJDD(3) * mrSges(5,1) - t59 * mrSges(5,3) + qJD(3) * t65 - t72 * t53 + t132;
t77 = (-mrSges(4,1) * t99 + mrSges(4,2) * t96) * qJD(2);
t7 = m(4) * t111 + qJDD(3) * mrSges(4,1) - t78 * mrSges(4,3) + qJD(3) * t84 + t89 * t11 - t77 * t124 + t92 * t14;
t8 = m(4) * t127 - qJDD(3) * mrSges(4,2) + t79 * mrSges(4,3) - qJD(3) * t83 + t92 * t11 + t77 * t123 - t89 * t14;
t4 = m(3) * t116 - t102 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t96 * t7 + t99 * t8;
t6 = m(3) * t67 + t99 * t7 + t96 * t8;
t115 = m(2) * t88 + t91 * t125 + t4 * t130 + t94 * t6;
t2 = m(2) * t81 - t97 * t10 + t100 * t4;
t1 = m(2) * t80 - t91 * t6 + (t4 * t97 + t125) * t94;
t3 = [-m(1) * g(1) - t90 * t1 + t93 * t2, t2, t4, t8, t11, t16, -t57 * mrSges(7,2) - t70 * t48 + t119; -m(1) * g(2) + t93 * t1 + t90 * t2, t1, t10, t7, t14, t13, -t39 * mrSges(7,3) - t62 * t41 + t120; -m(1) * g(3) + t115, t115, t6, t133, -t107, -t132, -t38 * mrSges(7,1) - t61 * t45 + t118;];
f_new  = t3;
