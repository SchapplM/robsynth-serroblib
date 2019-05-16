% Calculate vector of cutting forces with Newton-Euler
% S6PRRPRP2
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
% Datum: 2019-05-05 03:50
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:43:29
% EndTime: 2019-05-05 03:43:35
% DurationCPUTime: 2.43s
% Computational Cost: add. (30180->174), mult. (64409->230), div. (0->0), fcn. (45254->12), ass. (0->91)
t137 = -2 * qJD(4);
t100 = qJD(2) ^ 2;
t90 = sin(pkin(6));
t96 = sin(qJ(2));
t128 = t90 * t96;
t89 = sin(pkin(10));
t92 = cos(pkin(10));
t77 = t89 * g(1) - t92 * g(2);
t93 = cos(pkin(6));
t129 = t77 * t93;
t78 = -t92 * g(1) - t89 * g(2);
t87 = -g(3) + qJDD(1);
t98 = cos(qJ(2));
t114 = t87 * t128 + t96 * t129 + t98 * t78;
t48 = -t100 * pkin(2) + qJDD(2) * pkin(8) + t114;
t63 = -t90 * t77 + t93 * t87;
t95 = sin(qJ(3));
t97 = cos(qJ(3));
t110 = -t95 * t48 + t97 * t63;
t118 = qJD(2) * qJD(3);
t112 = t97 * t118;
t75 = t95 * qJDD(2) + t112;
t29 = (-t75 + t112) * qJ(4) + (t100 * t95 * t97 + qJDD(3)) * pkin(3) + t110;
t123 = t97 * t48 + t95 * t63;
t76 = t97 * qJDD(2) - t95 * t118;
t121 = qJD(2) * t95;
t79 = qJD(3) * pkin(3) - qJ(4) * t121;
t86 = t97 ^ 2;
t30 = -t86 * t100 * pkin(3) + t76 * qJ(4) - qJD(3) * t79 + t123;
t88 = sin(pkin(11));
t91 = cos(pkin(11));
t69 = (t88 * t97 + t91 * t95) * qJD(2);
t136 = t137 * t69 + t91 * t29 - t88 * t30;
t68 = (t88 * t95 - t91 * t97) * qJD(2);
t135 = (t87 * t90 + t129) * t98 - t96 * t78;
t103 = -qJDD(2) * pkin(2) - t135;
t101 = -t76 * pkin(3) + qJDD(4) + t79 * t121 + (-qJ(4) * t86 - pkin(8)) * t100 + t103;
t131 = cos(qJ(5));
t115 = t137 * t68 + t88 * t29 + t91 * t30;
t51 = t68 * pkin(4) - t69 * pkin(9);
t99 = qJD(3) ^ 2;
t23 = -t99 * pkin(4) + qJDD(3) * pkin(9) - t68 * t51 + t115;
t55 = -t88 * t75 + t91 * t76;
t56 = t91 * t75 + t88 * t76;
t25 = (qJD(3) * t68 - t56) * pkin(9) + (qJD(3) * t69 - t55) * pkin(4) + t101;
t94 = sin(qJ(5));
t125 = t131 * t23 + t94 * t25;
t57 = -t131 * qJD(3) + t94 * t69;
t58 = t94 * qJD(3) + t131 * t69;
t38 = t57 * pkin(5) - t58 * qJ(6);
t67 = qJD(5) + t68;
t46 = -t67 * mrSges(7,1) + t58 * mrSges(7,2);
t54 = qJDD(5) - t55;
t66 = t67 ^ 2;
t117 = m(7) * (-t66 * pkin(5) + t54 * qJ(6) + 0.2e1 * qJD(6) * t67 - t57 * t38 + t125) + t67 * t46 + t54 * mrSges(7,3);
t39 = t57 * mrSges(7,1) - t58 * mrSges(7,3);
t124 = -t57 * mrSges(6,1) - t58 * mrSges(6,2) - t39;
t126 = -mrSges(6,3) - mrSges(7,2);
t35 = t58 * qJD(5) - t131 * qJDD(3) + t94 * t56;
t45 = t67 * mrSges(6,1) - t58 * mrSges(6,3);
t14 = m(6) * t125 - t54 * mrSges(6,2) + t124 * t57 + t126 * t35 - t67 * t45 + t117;
t106 = t131 * t25 - t94 * t23;
t132 = m(7) * (-t54 * pkin(5) - t66 * qJ(6) + t58 * t38 + qJDD(6) - t106);
t36 = -t57 * qJD(5) + t94 * qJDD(3) + t131 * t56;
t43 = -t57 * mrSges(7,2) + t67 * mrSges(7,3);
t44 = -t67 * mrSges(6,2) - t57 * mrSges(6,3);
t16 = m(6) * t106 - t132 + (t44 + t43) * t67 + t124 * t58 + (mrSges(6,1) + mrSges(7,1)) * t54 + t126 * t36;
t61 = -qJD(3) * mrSges(5,2) - t68 * mrSges(5,3);
t62 = qJD(3) * mrSges(5,1) - t69 * mrSges(5,3);
t105 = -m(5) * t101 + t55 * mrSges(5,1) - t56 * mrSges(5,2) - t131 * t16 - t94 * t14 - t68 * t61 - t69 * t62;
t80 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t121;
t120 = qJD(2) * t97;
t81 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t120;
t134 = (t95 * t80 - t97 * t81) * qJD(2) + m(4) * (-t100 * pkin(8) + t103) - t76 * mrSges(4,1) + t75 * mrSges(4,2) - t105;
t22 = -qJDD(3) * pkin(4) - t99 * pkin(9) + t69 * t51 - t136;
t116 = m(7) * (-0.2e1 * qJD(6) * t58 + (t57 * t67 - t36) * qJ(6) + (t58 * t67 + t35) * pkin(5) + t22) + t57 * t43 + t35 * mrSges(7,1);
t133 = m(6) * t22 + t35 * mrSges(6,1) + (t45 - t46) * t58 + (mrSges(6,2) - mrSges(7,3)) * t36 + t57 * t44 + t116;
t10 = m(3) * t135 + qJDD(2) * mrSges(3,1) - t100 * mrSges(3,2) - t134;
t130 = t10 * t98;
t50 = t68 * mrSges(5,1) + t69 * mrSges(5,2);
t11 = m(5) * t115 - qJDD(3) * mrSges(5,2) + t55 * mrSges(5,3) - qJD(3) * t62 + t131 * t14 - t94 * t16 - t68 * t50;
t12 = m(5) * t136 + qJDD(3) * mrSges(5,1) - t56 * mrSges(5,3) + qJD(3) * t61 - t69 * t50 - t133;
t74 = (-mrSges(4,1) * t97 + mrSges(4,2) * t95) * qJD(2);
t7 = m(4) * t110 + qJDD(3) * mrSges(4,1) - t75 * mrSges(4,3) + qJD(3) * t81 + t88 * t11 + t91 * t12 - t74 * t121;
t8 = m(4) * t123 - qJDD(3) * mrSges(4,2) + t76 * mrSges(4,3) - qJD(3) * t80 + t91 * t11 - t88 * t12 + t74 * t120;
t4 = m(3) * t114 - t100 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t95 * t7 + t97 * t8;
t6 = m(3) * t63 + t97 * t7 + t95 * t8;
t113 = m(2) * t87 + t4 * t128 + t90 * t130 + t93 * t6;
t2 = m(2) * t78 - t96 * t10 + t98 * t4;
t1 = m(2) * t77 - t90 * t6 + (t4 * t96 + t130) * t93;
t3 = [-m(1) * g(1) - t89 * t1 + t92 * t2, t2, t4, t8, t11, t14, -t35 * mrSges(7,2) - t57 * t39 + t117; -m(1) * g(2) + t92 * t1 + t89 * t2, t1, t10, t7, t12, t16, -t36 * mrSges(7,3) - t58 * t46 + t116; -m(1) * g(3) + t113, t113, t6, t134, -t105, t133, -t54 * mrSges(7,1) + t36 * mrSges(7,2) + t58 * t39 - t67 * t43 + t132;];
f_new  = t3;
