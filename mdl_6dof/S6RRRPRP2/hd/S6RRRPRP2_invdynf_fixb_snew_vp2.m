% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRP2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:30:26
% EndTime: 2019-05-07 07:30:36
% DurationCPUTime: 3.90s
% Computational Cost: add. (49579->204), mult. (111331->262), div. (0->0), fcn. (81299->10), ass. (0->97)
t101 = sin(qJ(3));
t104 = cos(qJ(3));
t102 = sin(qJ(2));
t105 = cos(qJ(2));
t125 = qJD(1) * qJD(2);
t107 = qJD(1) ^ 2;
t103 = sin(qJ(1));
t106 = cos(qJ(1));
t117 = -g(1) * t106 - g(2) * t103;
t83 = -pkin(1) * t107 + qJDD(1) * pkin(7) + t117;
t129 = t102 * t83;
t136 = pkin(2) * t107;
t86 = qJDD(1) * t102 + t105 * t125;
t55 = qJDD(2) * pkin(2) - t86 * pkin(8) - t129 + (pkin(8) * t125 + t102 * t136 - g(3)) * t105;
t120 = -g(3) * t102 + t105 * t83;
t87 = qJDD(1) * t105 - t102 * t125;
t127 = qJD(1) * t102;
t90 = qJD(2) * pkin(2) - pkin(8) * t127;
t97 = t105 ^ 2;
t56 = pkin(8) * t87 - qJD(2) * t90 - t97 * t136 + t120;
t118 = -t101 * t56 + t104 * t55;
t80 = (-t101 * t102 + t104 * t105) * qJD(1);
t62 = qJD(3) * t80 + t101 * t87 + t104 * t86;
t81 = (t101 * t105 + t102 * t104) * qJD(1);
t95 = qJDD(2) + qJDD(3);
t96 = qJD(2) + qJD(3);
t25 = (t80 * t96 - t62) * qJ(4) + (t80 * t81 + t95) * pkin(3) + t118;
t130 = t101 * t55 + t104 * t56;
t61 = -qJD(3) * t81 - t101 * t86 + t104 * t87;
t76 = pkin(3) * t96 - qJ(4) * t81;
t79 = t80 ^ 2;
t28 = -pkin(3) * t79 + t61 * qJ(4) - t76 * t96 + t130;
t98 = sin(pkin(10));
t99 = cos(pkin(10));
t73 = t80 * t98 + t81 * t99;
t142 = -0.2e1 * qJD(4) * t73 + t99 * t25 - t98 * t28;
t121 = t103 * g(1) - t106 * g(2);
t115 = -qJDD(1) * pkin(1) - t121;
t111 = -t87 * pkin(2) + t90 * t127 + (-pkin(8) * t97 - pkin(7)) * t107 + t115;
t100 = sin(qJ(5));
t109 = -t61 * pkin(3) - t79 * qJ(4) + t81 * t76 + qJDD(4) + t111;
t137 = cos(qJ(5));
t72 = t80 * t99 - t81 * t98;
t122 = 0.2e1 * qJD(4) * t72 + t98 * t25 + t99 * t28;
t51 = -pkin(4) * t72 - pkin(9) * t73;
t94 = t96 ^ 2;
t21 = -pkin(4) * t94 + pkin(9) * t95 + t51 * t72 + t122;
t40 = t61 * t99 - t62 * t98;
t41 = t61 * t98 + t62 * t99;
t23 = (-t72 * t96 - t41) * pkin(9) + (t73 * t96 - t40) * pkin(4) + t109;
t133 = t100 * t23 + t137 * t21;
t39 = qJDD(5) - t40;
t58 = t100 * t73 - t137 * t96;
t59 = t100 * t96 + t137 * t73;
t42 = pkin(5) * t58 - qJ(6) * t59;
t69 = qJD(5) - t72;
t48 = -mrSges(7,1) * t69 + mrSges(7,2) * t59;
t68 = t69 ^ 2;
t124 = m(7) * (-pkin(5) * t68 + qJ(6) * t39 + 0.2e1 * qJD(6) * t69 - t42 * t58 + t133) + t69 * t48 + t39 * mrSges(7,3);
t43 = mrSges(7,1) * t58 - mrSges(7,3) * t59;
t132 = -mrSges(6,1) * t58 - mrSges(6,2) * t59 - t43;
t134 = -mrSges(6,3) - mrSges(7,2);
t30 = t59 * qJD(5) + t100 * t41 - t137 * t95;
t47 = mrSges(6,1) * t69 - mrSges(6,3) * t59;
t12 = m(6) * t133 - t39 * mrSges(6,2) + t132 * t58 + t134 * t30 - t69 * t47 + t124;
t114 = -t100 * t21 + t137 * t23;
t138 = m(7) * (-t39 * pkin(5) - t68 * qJ(6) + t59 * t42 + qJDD(6) - t114);
t31 = -t58 * qJD(5) + t100 * t95 + t137 * t41;
t45 = -mrSges(7,2) * t58 + mrSges(7,3) * t69;
t46 = -mrSges(6,2) * t69 - mrSges(6,3) * t58;
t14 = m(6) * t114 - t138 + (t46 + t45) * t69 + t132 * t59 + (mrSges(6,1) + mrSges(7,1)) * t39 + t134 * t31;
t64 = -mrSges(5,2) * t96 + mrSges(5,3) * t72;
t65 = mrSges(5,1) * t96 - mrSges(5,3) * t73;
t113 = -m(5) * t109 + t40 * mrSges(5,1) - t41 * mrSges(5,2) - t100 * t12 - t137 * t14 + t72 * t64 - t73 * t65;
t75 = -mrSges(4,2) * t96 + mrSges(4,3) * t80;
t77 = mrSges(4,1) * t96 - mrSges(4,3) * t81;
t110 = -m(4) * t111 + t61 * mrSges(4,1) - t62 * mrSges(4,2) + t80 * t75 - t81 * t77 + t113;
t88 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t127;
t126 = qJD(1) * t105;
t89 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t126;
t141 = (t102 * t88 - t105 * t89) * qJD(1) + m(3) * (-t107 * pkin(7) + t115) - t87 * mrSges(3,1) + t86 * mrSges(3,2) - t110;
t20 = -t95 * pkin(4) - t94 * pkin(9) + t73 * t51 - t142;
t123 = m(7) * (-0.2e1 * qJD(6) * t59 + (t58 * t69 - t31) * qJ(6) + (t59 * t69 + t30) * pkin(5) + t20) + t30 * mrSges(7,1) + t58 * t45;
t140 = m(6) * t20 + t30 * mrSges(6,1) + (t47 - t48) * t59 + (mrSges(6,2) - mrSges(7,3)) * t31 + t58 * t46 + t123;
t50 = -mrSges(5,1) * t72 + mrSges(5,2) * t73;
t10 = m(5) * t142 + t95 * mrSges(5,1) - t41 * mrSges(5,3) - t73 * t50 + t96 * t64 - t140;
t74 = -mrSges(4,1) * t80 + mrSges(4,2) * t81;
t9 = m(5) * t122 - t95 * mrSges(5,2) + t40 * mrSges(5,3) - t100 * t14 + t137 * t12 + t72 * t50 - t96 * t65;
t6 = m(4) * t118 + t95 * mrSges(4,1) - t62 * mrSges(4,3) + t99 * t10 - t81 * t74 + t96 * t75 + t98 * t9;
t7 = m(4) * t130 - t95 * mrSges(4,2) + t61 * mrSges(4,3) - t98 * t10 + t80 * t74 - t96 * t77 + t99 * t9;
t85 = (-mrSges(3,1) * t105 + mrSges(3,2) * t102) * qJD(1);
t4 = m(3) * (-t105 * g(3) - t129) - t86 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t85 * t127 + qJD(2) * t89 + t101 * t7 + t104 * t6;
t5 = m(3) * t120 - qJDD(2) * mrSges(3,2) + t87 * mrSges(3,3) - qJD(2) * t88 - t101 * t6 + t104 * t7 + t85 * t126;
t139 = t102 * t5 + t105 * t4;
t8 = m(2) * t121 + qJDD(1) * mrSges(2,1) - t107 * mrSges(2,2) - t141;
t1 = m(2) * t117 - t107 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t102 * t4 + t105 * t5;
t2 = [-m(1) * g(1) + t1 * t106 - t103 * t8, t1, t5, t7, t9, t12, -t30 * mrSges(7,2) - t58 * t43 + t124; -m(1) * g(2) + t1 * t103 + t106 * t8, t8, t4, t6, t10, t14, -t31 * mrSges(7,3) - t59 * t48 + t123; (-m(1) - m(2)) * g(3) + t139, -m(2) * g(3) + t139, t141, -t110, -t113, t140, -t39 * mrSges(7,1) + t31 * mrSges(7,2) + t59 * t43 - t69 * t45 + t138;];
f_new  = t2;
