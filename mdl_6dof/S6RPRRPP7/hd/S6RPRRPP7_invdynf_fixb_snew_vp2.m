% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-05-05 21:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPP7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:46:21
% EndTime: 2019-05-05 21:46:23
% DurationCPUTime: 0.77s
% Computational Cost: add. (6695->176), mult. (12661->197), div. (0->0), fcn. (7059->6), ass. (0->81)
t85 = sin(qJ(1));
t87 = cos(qJ(1));
t103 = -t87 * g(1) - t85 * g(2);
t134 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t103;
t86 = cos(qJ(3));
t116 = qJD(1) * t86;
t126 = cos(qJ(4));
t83 = sin(qJ(4));
t67 = t83 * qJD(3) + t116 * t126;
t84 = sin(qJ(3));
t115 = t84 * qJD(1);
t76 = qJD(4) + t115;
t48 = -t76 * mrSges(6,1) + t67 * mrSges(6,2);
t114 = qJD(1) * qJD(3);
t106 = t86 * t114;
t70 = -t84 * qJDD(1) - t106;
t65 = qJDD(4) - t70;
t107 = t84 * t114;
t71 = t86 * qJDD(1) - t107;
t127 = (-pkin(1) - pkin(7));
t89 = qJD(1) ^ 2;
t91 = (t127 * t89) - t134;
t20 = (-t71 + t107) * pkin(8) + (-t70 + t106) * pkin(3) + t91;
t109 = t85 * g(1) - t87 * g(2);
t97 = -t89 * qJ(2) + qJDD(2) - t109;
t54 = qJDD(1) * t127 + t97;
t108 = -t86 * g(3) + t84 * t54;
t69 = (pkin(3) * t84 - pkin(8) * t86) * qJD(1);
t88 = qJD(3) ^ 2;
t24 = -t88 * pkin(3) + qJDD(3) * pkin(8) - t115 * t69 + t108;
t120 = t126 * t24 + t83 * t20;
t129 = 2 * qJD(5);
t66 = -qJD(3) * t126 + t116 * t83;
t39 = t66 * pkin(4) - t67 * qJ(5);
t75 = t76 ^ 2;
t98 = -t75 * pkin(4) + t65 * qJ(5) + t76 * t129 - t66 * t39 + t120;
t133 = m(6) * t98 + t65 * mrSges(6,3) + t76 * t48;
t125 = t66 * t76;
t36 = -t66 * qJD(4) + t83 * qJDD(3) + t126 * t71;
t117 = t84 * g(3) + t86 * t54;
t93 = qJDD(3) * pkin(3) + t88 * pkin(8) - t69 * t116 + t117;
t132 = (-t36 + t125) * qJ(5) - t93;
t35 = t67 * qJD(4) - qJDD(3) * t126 + t83 * t71;
t44 = -t76 * mrSges(5,2) - t66 * mrSges(5,3);
t47 = t76 * mrSges(5,1) - t67 * mrSges(5,3);
t43 = t76 * mrSges(7,2) + t66 * mrSges(7,3);
t45 = -t76 * pkin(5) - t67 * qJ(6);
t46 = -t76 * mrSges(7,1) - t67 * mrSges(7,3);
t64 = t66 ^ 2;
t104 = m(7) * (-t64 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t35 + (-pkin(4) * t76 + t129 + t45) * t67 - t132) + t36 * mrSges(7,2) - t35 * mrSges(7,1) + t67 * t46 - t66 * t43;
t130 = -0.2e1 * t67;
t49 = -t66 * mrSges(6,2) + t76 * mrSges(6,3);
t96 = m(6) * (qJD(5) * t130 + (t67 * t76 + t35) * pkin(4) + t132) + t35 * mrSges(6,1) + t66 * t49 - t104;
t131 = -m(5) * t93 + t35 * mrSges(5,1) + (t47 - t48) * t67 + (mrSges(5,2) - mrSges(6,3)) * t36 + t66 * t44 + t96;
t128 = -m(2) - m(3);
t124 = (mrSges(2,1) - mrSges(3,2));
t123 = -mrSges(2,2) + mrSges(3,3);
t121 = -mrSges(5,3) - mrSges(6,2);
t40 = t66 * mrSges(6,1) - t67 * mrSges(6,3);
t119 = -t66 * mrSges(5,1) - t67 * mrSges(5,2) - t40;
t41 = -t66 * mrSges(7,1) + t67 * mrSges(7,2);
t112 = m(7) * (-t64 * pkin(5) + t35 * qJ(6) + 0.2e1 * qJD(6) * t66 + t76 * t45 + t98) + t66 * t41 + t35 * mrSges(7,3);
t102 = t126 * t20 - t83 * t24;
t16 = -t65 * pkin(4) - t75 * qJ(5) + t67 * t39 + qJDD(5) - t102;
t111 = t76 * t43 + t65 * mrSges(7,1) - m(7) * (qJD(6) * t130 + (-t36 - t125) * qJ(6) + (t66 * t67 - t65) * pkin(5) + t16);
t68 = (mrSges(4,1) * t84 + mrSges(4,2) * t86) * qJD(1);
t73 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t116;
t101 = m(6) * t16 - t111;
t8 = m(5) * t102 + (t44 + t49) * t76 + (mrSges(5,1) + mrSges(6,1)) * t65 + (t41 + t119) * t67 + (mrSges(7,3) + t121) * t36 - t101;
t9 = m(5) * t120 + (-t47 + t46) * t76 + t119 * t66 + (-mrSges(5,2) + mrSges(7,2)) * t65 + t121 * t35 + t112 + t133;
t4 = m(4) * t108 - qJDD(3) * mrSges(4,2) + t70 * mrSges(4,3) - qJD(3) * t73 - t115 * t68 + t126 * t9 - t83 * t8;
t72 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t115;
t5 = m(4) * t117 + qJDD(3) * mrSges(4,1) - t71 * mrSges(4,3) + qJD(3) * t72 - t116 * t68 - t131;
t110 = t86 * t4 - t84 * t5;
t99 = -m(3) * (-qJDD(1) * pkin(1) + t97) - t84 * t4 - t86 * t5;
t95 = t65 * mrSges(7,2) + t76 * t46 + t112;
t94 = m(4) * t91 - t70 * mrSges(4,1) + t71 * mrSges(4,2) + t72 * t115 + t73 * t116 + t126 * t8 + t83 * t9;
t92 = -m(3) * (t89 * pkin(1) + t134) + t94;
t2 = m(2) * t103 + qJDD(1) * t123 - (t124 * t89) + t92;
t1 = m(2) * t109 + qJDD(1) * t124 + t123 * t89 + t99;
t3 = [-m(1) * g(1) - t85 * t1 + t87 * t2, t2, -m(3) * g(3) + t110, t4, t9, -t35 * mrSges(6,2) - t66 * t40 + t133 + t95, t95; -m(1) * g(2) + t87 * t1 + t85 * t2, t1, -(t89 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t92, t5, t8, -t36 * mrSges(6,3) - t67 * t48 + t96, -t36 * mrSges(7,3) - t67 * t41 - t111; (-m(1) + t128) * g(3) + t110, g(3) * t128 + t110, qJDD(1) * mrSges(3,2) - t89 * mrSges(3,3) - t99, t94, t131, -t65 * mrSges(6,1) - t76 * t49 + (t40 - t41) * t67 + (mrSges(6,2) - mrSges(7,3)) * t36 + t101, t104;];
f_new  = t3;
