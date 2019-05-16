% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 18:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:23:50
% EndTime: 2019-05-07 18:23:54
% DurationCPUTime: 1.84s
% Computational Cost: add. (22504->201), mult. (44623->240), div. (0->0), fcn. (30273->8), ass. (0->93)
t137 = cos(qJ(4));
t100 = sin(qJ(3));
t103 = cos(qJ(3));
t104 = cos(qJ(2));
t128 = qJD(1) * qJD(2);
t122 = t104 * t128;
t107 = qJD(1) ^ 2;
t102 = sin(qJ(1));
t105 = cos(qJ(1));
t125 = t102 * g(1) - t105 * g(2);
t78 = -qJDD(1) * pkin(1) - t107 * pkin(7) - t125;
t101 = sin(qJ(2));
t87 = t101 * qJDD(1) + t122;
t96 = t101 * t128;
t88 = t104 * qJDD(1) - t96;
t47 = (-t87 - t122) * pkin(8) + (-t88 + t96) * pkin(2) + t78;
t106 = qJD(2) ^ 2;
t118 = -t105 * g(1) - t102 * g(2);
t79 = -t107 * pkin(1) + qJDD(1) * pkin(7) + t118;
t124 = -t101 * g(3) + t104 * t79;
t129 = t104 * qJD(1);
t86 = (-pkin(2) * t104 - pkin(8) * t101) * qJD(1);
t52 = -t106 * pkin(2) + qJDD(2) * pkin(8) + t129 * t86 + t124;
t121 = -t100 * t52 + t103 * t47;
t130 = qJD(1) * t101;
t83 = t103 * qJD(2) - t100 * t130;
t63 = t83 * qJD(3) + t100 * qJDD(2) + t103 * t87;
t82 = qJDD(3) - t88;
t84 = t100 * qJD(2) + t103 * t130;
t95 = qJD(3) - t129;
t21 = (t83 * t95 - t63) * pkin(9) + (t83 * t84 + t82) * pkin(3) + t121;
t132 = t100 * t47 + t103 * t52;
t62 = -t84 * qJD(3) + t103 * qJDD(2) - t100 * t87;
t70 = t95 * pkin(3) - t84 * pkin(9);
t81 = t83 ^ 2;
t24 = -t81 * pkin(3) + t62 * pkin(9) - t95 * t70 + t132;
t99 = sin(qJ(4));
t134 = t137 * t24 + t99 * t21;
t139 = -2 * qJD(5);
t65 = -t137 * t83 + t99 * t84;
t66 = t137 * t84 + t99 * t83;
t41 = t65 * pkin(4) - t66 * qJ(5);
t80 = -qJDD(4) - t82;
t93 = -qJD(4) - t95;
t92 = t93 ^ 2;
t115 = -t92 * pkin(4) - t80 * qJ(5) + t139 * t93 - t65 * t41 + t134;
t59 = t93 * mrSges(6,1) + t66 * mrSges(6,2);
t142 = m(6) * t115 - t80 * mrSges(6,3) - t93 * t59;
t131 = -t104 * g(3) - t101 * t79;
t113 = qJDD(2) * pkin(2) + t106 * pkin(8) - t86 * t130 + t131;
t110 = t62 * pkin(3) + t81 * pkin(9) - t84 * t70 + t113;
t136 = t65 * t93;
t33 = -t65 * qJD(4) + t137 * t63 + t99 * t62;
t141 = -(t33 + t136) * qJ(5) - t110;
t32 = t66 * qJD(4) - t137 * t62 + t99 * t63;
t43 = -t65 * mrSges(7,1) + t66 * mrSges(7,2);
t56 = t93 * pkin(5) - t66 * qJ(6);
t64 = t65 ^ 2;
t126 = m(7) * (-t64 * pkin(5) + t32 * qJ(6) + 0.2e1 * qJD(6) * t65 - t93 * t56 + t115) + t32 * mrSges(7,3) + t65 * t43;
t42 = t65 * mrSges(6,1) - t66 * mrSges(6,3);
t133 = -t65 * mrSges(5,1) - t66 * mrSges(5,2) - t42;
t135 = -mrSges(5,3) - mrSges(6,2);
t57 = t93 * mrSges(7,1) - t66 * mrSges(7,3);
t58 = -t93 * mrSges(5,1) - t66 * mrSges(5,3);
t10 = m(5) * t134 + (t58 - t57) * t93 + (mrSges(5,2) - mrSges(7,2)) * t80 + t133 * t65 + t135 * t32 + t126 + t142;
t67 = -t83 * mrSges(4,1) + t84 * mrSges(4,2);
t68 = -t95 * mrSges(4,2) + t83 * mrSges(4,3);
t119 = t137 * t21 - t99 * t24;
t17 = t80 * pkin(4) - t92 * qJ(5) + t66 * t41 + qJDD(5) - t119;
t54 = -t93 * mrSges(7,2) + t65 * mrSges(7,3);
t127 = m(7) * (-0.2e1 * qJD(6) * t66 + (-t33 + t136) * qJ(6) + (t65 * t66 + t80) * pkin(5) + t17) + t93 * t54 + t80 * mrSges(7,1);
t117 = m(6) * t17 + t127;
t53 = -t65 * mrSges(6,2) - t93 * mrSges(6,3);
t55 = t93 * mrSges(5,2) - t65 * mrSges(5,3);
t9 = m(5) * t119 + (-t55 - t53) * t93 + (-mrSges(5,1) - mrSges(6,1)) * t80 + (t43 + t133) * t66 + (mrSges(7,3) + t135) * t33 - t117;
t5 = m(4) * t121 + t82 * mrSges(4,1) - t63 * mrSges(4,3) + t99 * t10 + t137 * t9 - t84 * t67 + t95 * t68;
t69 = t95 * mrSges(4,1) - t84 * mrSges(4,3);
t6 = m(4) * t132 - t82 * mrSges(4,2) + t62 * mrSges(4,3) + t10 * t137 + t83 * t67 - t95 * t69 - t99 * t9;
t90 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t130;
t91 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t129;
t140 = m(3) * t78 - t88 * mrSges(3,1) + t87 * mrSges(3,2) + t100 * t6 + t103 * t5 + (t101 * t90 - t104 * t91) * qJD(1);
t85 = (-mrSges(3,1) * t104 + mrSges(3,2) * t101) * qJD(1);
t4 = m(3) * t124 - qJDD(2) * mrSges(3,2) + t88 * mrSges(3,3) - qJD(2) * t90 - t100 * t5 + t103 * t6 + t129 * t85;
t120 = m(7) * (-t64 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t32 + (pkin(4) * t93 + (2 * qJD(5)) + t56) * t66 - t141) + t33 * mrSges(7,2) - t32 * mrSges(7,1) + t66 * t57 - t65 * t54;
t111 = t33 * mrSges(6,3) + t66 * t59 + t120 - m(6) * (t66 * t139 + (-t66 * t93 + t32) * pkin(4) + t141) - t32 * mrSges(6,1) - t65 * t53;
t109 = -m(5) * t110 + t32 * mrSges(5,1) + t33 * mrSges(5,2) + t65 * t55 + t66 * t58 - t111;
t108 = -m(4) * t113 - t62 * mrSges(4,1) + t63 * mrSges(4,2) - t83 * t68 + t84 * t69 + t109;
t8 = m(3) * t131 + qJDD(2) * mrSges(3,1) - t87 * mrSges(3,3) + qJD(2) * t91 - t130 * t85 - t108;
t138 = t101 * t4 + t104 * t8;
t114 = -t80 * mrSges(7,2) - t93 * t57 + t126;
t2 = m(2) * t125 + qJDD(1) * mrSges(2,1) - t107 * mrSges(2,2) - t140;
t1 = m(2) * t118 - t107 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t101 * t8 + t104 * t4;
t3 = [-m(1) * g(1) + t105 * t1 - t102 * t2, t1, t4, t6, t10, -t32 * mrSges(6,2) - t65 * t42 + t114 + t142, t114; -m(1) * g(2) + t102 * t1 + t105 * t2, t2, t8, t5, t9, -t111, -t33 * mrSges(7,3) - t66 * t43 + t127; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t140, t108, t109, t80 * mrSges(6,1) + t93 * t53 + (t42 - t43) * t66 + (mrSges(6,2) - mrSges(7,3)) * t33 + t117, t120;];
f_new  = t3;
