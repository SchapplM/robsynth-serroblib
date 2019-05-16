% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-05-06 12:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:35:29
% EndTime: 2019-05-06 12:35:34
% DurationCPUTime: 1.72s
% Computational Cost: add. (18191->201), mult. (38178->245), div. (0->0), fcn. (22793->8), ass. (0->93)
t147 = -2 * qJD(3);
t100 = cos(qJ(2));
t101 = cos(qJ(1));
t98 = sin(qJ(1));
t122 = t98 * g(1) - t101 * g(2);
t115 = -qJDD(1) * pkin(1) - t122;
t126 = qJD(1) * qJD(2);
t118 = t100 * t126;
t97 = sin(qJ(2));
t120 = t97 * t126;
t128 = t97 * qJD(1);
t77 = t97 * qJDD(1) + t118;
t108 = pkin(2) * t120 + t128 * t147 + (-t77 - t118) * qJ(3) + t115;
t127 = qJD(1) * t100;
t103 = qJD(1) ^ 2;
t139 = t103 * pkin(7);
t78 = t100 * qJDD(1) - t120;
t84 = pkin(3) * t128 - qJD(2) * pkin(8);
t94 = t100 ^ 2;
t23 = -t84 * t128 + (-pkin(2) - pkin(8)) * t78 + (-pkin(3) * t94 - pkin(7)) * t103 + t108;
t102 = qJD(2) ^ 2;
t117 = -t101 * g(1) - t98 * g(2);
t66 = -t103 * pkin(1) + qJDD(1) * pkin(7) + t117;
t133 = -t100 * g(3) - t97 * t66;
t74 = (-pkin(2) * t100 - qJ(3) * t97) * qJD(1);
t42 = -qJDD(2) * pkin(2) - t102 * qJ(3) + t74 * t128 + qJDD(3) - t133;
t30 = (-t100 * t103 * t97 - qJDD(2)) * pkin(8) + (t77 - t118) * pkin(3) + t42;
t96 = sin(qJ(4));
t99 = cos(qJ(4));
t119 = -t96 * t23 + t99 * t30;
t129 = cos(pkin(9));
t72 = -t96 * qJD(2) - t127 * t99;
t53 = t72 * qJD(4) + t99 * qJDD(2) - t96 * t78;
t73 = t99 * qJD(2) - t127 * t96;
t56 = -t72 * mrSges(5,1) + t73 * mrSges(5,2);
t87 = qJD(4) + t128;
t57 = -t87 * mrSges(5,2) + t72 * mrSges(5,3);
t71 = qJDD(4) + t77;
t142 = -2 * qJD(5);
t17 = (t72 * t87 - t53) * qJ(5) + (t72 * t73 + t71) * pkin(4) + t119;
t135 = t99 * t23 + t96 * t30;
t52 = -t73 * qJD(4) - t96 * qJDD(2) - t99 * t78;
t58 = t87 * pkin(4) - t73 * qJ(5);
t70 = t72 ^ 2;
t19 = -t70 * pkin(4) + t52 * qJ(5) - t87 * t58 + t135;
t95 = sin(pkin(9));
t54 = -t129 * t72 + t95 * t73;
t123 = t129 * t19 + t54 * t142 + t95 * t17;
t55 = t129 * t73 + t95 * t72;
t35 = t54 * pkin(5) - t55 * qJ(6);
t46 = -t87 * mrSges(7,1) + t55 * mrSges(7,2);
t85 = t87 ^ 2;
t124 = m(7) * (-t85 * pkin(5) + t71 * qJ(6) + 0.2e1 * qJD(6) * t87 - t54 * t35 + t123) + t87 * t46 + t71 * mrSges(7,3);
t36 = t54 * mrSges(7,1) - t55 * mrSges(7,3);
t134 = -t54 * mrSges(6,1) - t55 * mrSges(6,2) - t36;
t136 = -mrSges(6,3) - mrSges(7,2);
t31 = -t129 * t52 + t95 * t53;
t45 = t87 * mrSges(6,1) - t55 * mrSges(6,3);
t8 = m(6) * t123 - t71 * mrSges(6,2) + t134 * t54 + t136 * t31 - t87 * t45 + t124;
t114 = t129 * t17 - t95 * t19;
t141 = m(7) * (-t71 * pkin(5) - t85 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t35) * t55 - t114);
t32 = t129 * t53 + t95 * t52;
t43 = -t87 * mrSges(6,2) - t54 * mrSges(6,3);
t44 = -t54 * mrSges(7,2) + t87 * mrSges(7,3);
t9 = m(6) * t114 - t141 + (t43 + t44) * t87 + (mrSges(6,1) + mrSges(7,1)) * t71 + (m(6) * t142 + t134) * t55 + t136 * t32;
t6 = m(5) * t119 + t71 * mrSges(5,1) - t53 * mrSges(5,3) + t129 * t9 - t73 * t56 + t87 * t57 + t95 * t8;
t59 = t87 * mrSges(5,1) - t73 * mrSges(5,3);
t7 = m(5) * t135 - t71 * mrSges(5,2) + t52 * mrSges(5,3) + t129 * t8 + t72 * t56 - t87 * t59 - t95 * t9;
t82 = -mrSges(4,1) * t127 - qJD(2) * mrSges(4,3);
t116 = t96 * t6 - t99 * t7 - m(4) * (-t78 * pkin(2) + t108 - t139) - t82 * t127 + t77 * mrSges(4,3);
t83 = mrSges(4,1) * t128 + qJD(2) * mrSges(4,2);
t131 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t128 - t83;
t138 = mrSges(3,1) - mrSges(4,2);
t81 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t127;
t146 = (-t100 * t81 + t131 * t97) * qJD(1) - t138 * t78 + m(3) * (t115 - t139) + t77 * mrSges(3,2) - t116;
t121 = -t97 * g(3) + t100 * t66;
t145 = t102 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t147 - t74 * t127 - t121;
t109 = -t94 * t103 * pkin(8) + t78 * pkin(3) + qJD(2) * t84 - t145;
t106 = -t52 * pkin(4) - t70 * qJ(5) + t73 * t58 + qJDD(5) + t109;
t111 = t32 * mrSges(7,3) + t55 * t46 - m(7) * (t106 + (t54 * t87 - t32) * qJ(6) + (t55 * t87 + t31) * pkin(5) - 0.2e1 * qJD(6) * t55) - t31 * mrSges(7,1) - t54 * t44;
t107 = m(6) * t106 + t31 * mrSges(6,1) + t32 * mrSges(6,2) + t54 * t43 + t55 * t45 - t111;
t105 = -m(5) * t109 + t52 * mrSges(5,1) - t53 * mrSges(5,2) + t72 * t57 - t73 * t59 - t107;
t104 = m(4) * t145 + t105;
t75 = (mrSges(4,2) * t100 - mrSges(4,3) * t97) * qJD(1);
t132 = t75 + (-mrSges(3,1) * t100 + mrSges(3,2) * t97) * qJD(1);
t137 = mrSges(3,3) + mrSges(4,1);
t11 = -t104 + t132 * t127 + t137 * t78 + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) - t131 * qJD(2) + m(3) * t121;
t113 = -m(4) * t42 - t99 * t6 - t96 * t7;
t4 = m(3) * t133 - t137 * t77 + t138 * qJDD(2) + (t81 - t82) * qJD(2) - t132 * t128 + t113;
t140 = t100 * t4 + t97 * t11;
t2 = m(2) * t122 + qJDD(1) * mrSges(2,1) - t103 * mrSges(2,2) - t146;
t1 = m(2) * t117 - t103 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t100 * t11 - t97 * t4;
t3 = [-m(1) * g(1) + t101 * t1 - t98 * t2, t1, t11, t78 * mrSges(4,2) - t128 * t83 - t116, t7, t8, -t31 * mrSges(7,2) - t54 * t36 + t124; -m(1) * g(2) + t98 * t1 + t101 * t2, t2, t4, -t78 * mrSges(4,1) - qJDD(2) * mrSges(4,3) - qJD(2) * t83 - t127 * t75 + t104, t6, t9, -t111; (-m(1) - m(2)) * g(3) + t140, -m(2) * g(3) + t140, t146, t77 * mrSges(4,1) + qJDD(2) * mrSges(4,2) + qJD(2) * t82 + t128 * t75 - t113, -t105, t107, -t71 * mrSges(7,1) + t32 * mrSges(7,2) + t55 * t36 - t87 * t44 + t141;];
f_new  = t3;
