% Calculate vector of cutting forces with Newton-Euler
% S6RPRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-05-05 17:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPPR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:22:14
% EndTime: 2019-05-05 17:22:16
% DurationCPUTime: 0.69s
% Computational Cost: add. (4320->186), mult. (8608->211), div. (0->0), fcn. (3819->6), ass. (0->89)
t94 = sin(qJ(1));
t97 = cos(qJ(1));
t120 = t94 * g(1) - t97 * g(2);
t99 = qJD(1) ^ 2;
t110 = -t99 * qJ(2) + qJDD(2) - t120;
t34 = (-pkin(1) - pkin(7)) * qJDD(1) + t110;
t93 = sin(qJ(3));
t96 = cos(qJ(3));
t119 = -t96 * g(3) + t93 * t34;
t144 = -qJDD(3) * qJ(4) - t119;
t126 = qJD(1) * t96;
t57 = (pkin(3) * t93 - qJ(4) * t96) * qJD(1);
t143 = t57 * t126 + qJDD(4);
t142 = -2 * qJD(5);
t141 = -m(2) - m(3);
t140 = -pkin(3) - pkin(8);
t139 = -pkin(4) - pkin(8);
t138 = t93 * g(3);
t98 = qJD(3) ^ 2;
t137 = t98 * pkin(3);
t136 = t93 ^ 2 * t99;
t135 = t96 * t34;
t134 = mrSges(2,1) - mrSges(3,2);
t133 = -mrSges(2,2) + mrSges(3,3);
t132 = -mrSges(4,3) - mrSges(5,2);
t131 = -pkin(5) - qJ(4);
t58 = (mrSges(5,1) * t93 - mrSges(5,3) * t96) * qJD(1);
t130 = -t58 - (mrSges(4,1) * t93 + mrSges(4,2) * t96) * qJD(1);
t129 = -t97 * g(1) - t94 * g(2);
t128 = t98 * qJ(4);
t127 = qJD(1) * t93;
t125 = qJD(4) * t96;
t124 = t142 + t57;
t123 = qJD(1) * qJD(3);
t122 = qJD(2) * qJD(1);
t118 = t99 * pkin(1) - qJDD(1) * qJ(2) - t129;
t117 = -t99 * pkin(7) - t118;
t61 = qJDD(1) * t93 + t96 * t123;
t77 = t93 * t123;
t62 = qJDD(1) * t96 - t77;
t113 = t61 * pkin(3) - t62 * qJ(4) + t117;
t68 = -qJD(3) * pkin(4) - qJ(5) * t126;
t80 = -0.2e1 * t122;
t104 = -qJ(5) * t136 + 0.2e1 * qJD(1) * t125 + t68 * t126 + qJDD(5) - t113 + t80;
t11 = t104 + (t131 * t93 + t140 * t96) * t123 + pkin(5) * t62 + t139 * t61;
t106 = t96 * t99 * t93 * pkin(4) + t126 * t142 + (-t62 - t77) * qJ(5) - t138 + t143;
t60 = (pkin(5) * t96 - pkin(8) * t93) * qJD(1);
t14 = t131 * t98 + (-qJD(1) * t60 - t34) * t96 + (-pkin(3) + t139) * qJDD(3) + t106;
t92 = sin(qJ(6));
t95 = cos(qJ(6));
t55 = -qJD(3) * t92 + t95 * t127;
t26 = -qJD(6) * t55 - qJDD(3) * t95 - t61 * t92;
t54 = -qJD(3) * t95 - t92 * t127;
t28 = -mrSges(7,1) * t54 + mrSges(7,2) * t55;
t74 = qJD(6) + t126;
t30 = mrSges(7,1) * t74 - mrSges(7,3) * t55;
t52 = qJDD(6) + t62;
t10 = m(7) * (t11 * t92 + t14 * t95) + t26 * mrSges(7,3) - t52 * mrSges(7,2) + t54 * t28 - t74 * t30;
t66 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t127;
t27 = qJD(6) * t54 - qJDD(3) * t92 + t61 * t95;
t29 = -mrSges(7,2) * t74 + mrSges(7,3) * t54;
t9 = m(7) * (t11 * t95 - t14 * t92) - t27 * mrSges(7,3) + t52 * mrSges(7,1) - t55 * t28 + t74 * t29;
t112 = -t95 * t10 + t92 * t9 - m(6) * (-t128 - t135 + (-pkin(3) - pkin(4)) * qJDD(3) + t106) + t62 * mrSges(6,3) - qJD(3) * t66 - qJDD(3) * mrSges(6,2);
t116 = t135 + t138;
t108 = m(5) * (-qJDD(3) * pkin(3) - t116 - t128 + t143) - t112;
t56 = (mrSges(6,1) * t96 + mrSges(6,2) * t93) * qJD(1);
t67 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t127;
t72 = -mrSges(5,2) * t127 + qJD(3) * mrSges(5,3);
t3 = m(4) * t116 + t132 * t62 + (mrSges(4,1) + mrSges(5,1)) * qJDD(3) + (t67 + t72) * qJD(3) + (t56 + t130) * t126 - t108;
t107 = pkin(4) * t136 - t61 * qJ(5) + t144;
t78 = 0.2e1 * qJD(4) * qJD(3);
t109 = -t26 * mrSges(7,1) - t54 * t29 + m(7) * (qJDD(3) * pkin(5) + qJD(3) * t68 + t78 + t140 * t98 + (t60 - t124) * t127 - t107) + t27 * mrSges(7,2) + t55 * t30;
t69 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t126;
t105 = m(6) * (t137 + (-0.2e1 * qJD(4) - t68) * qJD(3) + t124 * t127 + t107) - t109 - t56 * t127 - t61 * mrSges(6,3) - qJD(3) * t69 - qJDD(3) * mrSges(6,1);
t71 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t126;
t102 = -t105 + m(5) * (-t57 * t127 - t137 - t144 + t78) + qJD(3) * t71 + qJDD(3) * mrSges(5,3);
t70 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t126;
t5 = m(4) * t119 - qJDD(3) * mrSges(4,2) - qJD(3) * t70 + t130 * t127 + t132 * t61 + t102;
t121 = -t93 * t3 + t96 * t5;
t115 = pkin(3) * t96 + qJ(4) * t93;
t114 = t92 * t10 + t95 * t9 + m(6) * (-t61 * pkin(4) - t115 * t123 + t104) + t66 * t127 + t69 * t126 + t61 * mrSges(6,2) + t62 * mrSges(6,1);
t111 = -m(3) * (-qJDD(1) * pkin(1) + t110) - t96 * t3 - t93 * t5;
t79 = 0.2e1 * t122;
t103 = -t62 * mrSges(5,3) - t71 * t126 - t114 + m(5) * (t79 + (t115 * qJD(3) - 0.2e1 * t125) * qJD(1) + t113) + t72 * t127 + t61 * mrSges(5,1);
t101 = t61 * mrSges(4,1) + t103 + m(4) * (t79 + t117) + t67 * t127 + t70 * t126 + t62 * mrSges(4,2);
t100 = -m(3) * (t80 + t118) + t101;
t2 = m(2) * t129 + t133 * qJDD(1) - t134 * t99 + t100;
t1 = m(2) * t120 + t134 * qJDD(1) + t133 * t99 + t111;
t4 = [-m(1) * g(1) - t1 * t94 + t2 * t97, t2, -m(3) * g(3) + t121, t5, -t61 * mrSges(5,2) - t58 * t127 + t102, -t56 * t126 - t112, t10; -m(1) * g(2) + t1 * t97 + t2 * t94, t1, -t99 * mrSges(3,2) - qJDD(1) * mrSges(3,3) - t100, t3, t103, t105, t9; (-m(1) + t141) * g(3) + t121, t141 * g(3) + t121, qJDD(1) * mrSges(3,2) - t99 * mrSges(3,3) - t111, t101, -qJDD(3) * mrSges(5,1) + t62 * mrSges(5,2) - qJD(3) * t72 + (-t56 + t58) * t126 + t108, t114, t109;];
f_new  = t4;
