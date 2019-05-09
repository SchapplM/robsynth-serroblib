% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-05-06 12:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:41:45
% EndTime: 2019-05-06 12:41:48
% DurationCPUTime: 1.02s
% Computational Cost: add. (7183->203), mult. (14492->230), div. (0->0), fcn. (7784->6), ass. (0->90)
t100 = sin(qJ(2));
t103 = cos(qJ(2));
t101 = sin(qJ(1));
t104 = cos(qJ(1));
t127 = t101 * g(1) - t104 * g(2);
t117 = -qJDD(1) * pkin(1) - t127;
t102 = cos(qJ(4));
t131 = qJD(1) * qJD(2);
t124 = t103 * t131;
t125 = t100 * t131;
t132 = t100 * qJD(1);
t76 = t100 * qJDD(1) + t124;
t109 = pkin(2) * t125 - 0.2e1 * qJD(3) * t132 + (-t76 - t124) * qJ(3) + t117;
t133 = qJD(1) * t103;
t106 = qJD(1) ^ 2;
t145 = t106 * pkin(7);
t77 = t103 * qJDD(1) - t125;
t99 = sin(qJ(4));
t71 = t99 * qJD(2) + t102 * t133;
t88 = qJD(4) + t132;
t144 = t71 * t88;
t72 = t102 * qJD(2) - t99 * t133;
t148 = -0.2e1 * t72;
t84 = pkin(3) * t132 - qJD(2) * pkin(8);
t97 = t103 ^ 2;
t21 = -t84 * t132 + (-pkin(2) - pkin(8)) * t77 + (-pkin(3) * t97 - pkin(7)) * t106 + t109;
t105 = qJD(2) ^ 2;
t120 = -t104 * g(1) - t101 * g(2);
t62 = -t106 * pkin(1) + qJDD(1) * pkin(7) + t120;
t138 = -t103 * g(3) - t100 * t62;
t73 = (-pkin(2) * t103 - qJ(3) * t100) * qJD(1);
t29 = -qJDD(2) * pkin(2) - t105 * qJ(3) + t73 * t132 + qJDD(3) - t138;
t25 = (-t100 * t103 * t106 - qJDD(2)) * pkin(8) + (t76 - t124) * pkin(3) + t29;
t123 = t102 * t25 - t99 * t21;
t46 = t71 * pkin(4) - t72 * qJ(5);
t70 = qJDD(4) + t76;
t85 = t88 ^ 2;
t16 = -t70 * pkin(4) - t85 * qJ(5) + t72 * t46 + qJDD(5) - t123;
t43 = -t71 * qJD(4) + t102 * qJDD(2) - t99 * t77;
t50 = t88 * mrSges(7,2) + t71 * mrSges(7,3);
t128 = t88 * t50 + t70 * mrSges(7,1) - m(7) * (qJD(6) * t148 + (-t43 - t144) * qJ(6) + (t71 * t72 - t70) * pkin(5) + t16);
t119 = m(6) * t16 - t128;
t47 = t71 * mrSges(6,1) - t72 * mrSges(6,3);
t139 = -t71 * mrSges(5,1) - t72 * mrSges(5,2) - t47;
t141 = -mrSges(5,3) - mrSges(6,2);
t48 = -t71 * mrSges(7,1) + t72 * mrSges(7,2);
t51 = -t88 * mrSges(5,2) - t71 * mrSges(5,3);
t52 = -t71 * mrSges(6,2) + t88 * mrSges(6,3);
t8 = m(5) * t123 + (t51 + t52) * t88 + (mrSges(5,1) + mrSges(6,1)) * t70 + (t48 + t139) * t72 + (mrSges(7,3) + t141) * t43 - t119;
t81 = -mrSges(4,1) * t133 - qJD(2) * mrSges(4,3);
t140 = t102 * t21 + t99 * t25;
t147 = 2 * qJD(5);
t116 = -t85 * pkin(4) + t70 * qJ(5) + t88 * t147 - t71 * t46 + t140;
t42 = t72 * qJD(4) + t99 * qJDD(2) + t102 * t77;
t53 = -t88 * pkin(5) - t72 * qJ(6);
t69 = t71 ^ 2;
t129 = m(7) * (-t69 * pkin(5) + t42 * qJ(6) + 0.2e1 * qJD(6) * t71 + t88 * t53 + t116) + t71 * t48 + t42 * mrSges(7,3);
t56 = -t88 * mrSges(6,1) + t72 * mrSges(6,2);
t151 = m(6) * t116 + t70 * mrSges(6,3) + t88 * t56;
t54 = -t88 * mrSges(7,1) - t72 * mrSges(7,3);
t55 = t88 * mrSges(5,1) - t72 * mrSges(5,3);
t9 = m(5) * t140 + (-t55 + t54) * t88 + t139 * t71 + (-mrSges(5,2) + mrSges(7,2)) * t70 + t141 * t42 + t129 + t151;
t118 = -t102 * t9 + t99 * t8 - m(4) * (-t77 * pkin(2) + t109 - t145) - t81 * t133 + t76 * mrSges(4,3);
t82 = mrSges(4,1) * t132 + qJD(2) * mrSges(4,2);
t135 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t132 - t82;
t143 = mrSges(3,1) - mrSges(4,2);
t80 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t133;
t153 = (t135 * t100 - t103 * t80) * qJD(1) - t143 * t77 + m(3) * (t117 - t145) + t76 * mrSges(3,2) - t118;
t152 = (-t43 + t144) * qJ(5);
t115 = -m(4) * t29 - t102 * t8 - t99 * t9;
t74 = (mrSges(4,2) * t103 - mrSges(4,3) * t100) * qJD(1);
t136 = t74 + (-mrSges(3,1) * t103 + mrSges(3,2) * t100) * qJD(1);
t142 = mrSges(3,3) + mrSges(4,1);
t4 = m(3) * t138 - t142 * t76 + t143 * qJDD(2) + (t80 - t81) * qJD(2) - t136 * t132 + t115;
t137 = -t100 * g(3) + t103 * t62;
t121 = t105 * pkin(2) - qJDD(2) * qJ(3) - t73 * t133 - t137;
t112 = t97 * t106 * pkin(8) - t77 * pkin(3) - qJD(2) * t84 + t121;
t130 = qJD(3) * qJD(2);
t111 = -t112 + 0.2e1 * t130;
t92 = -0.2e1 * t130;
t122 = m(7) * (-t69 * qJ(6) + qJDD(6) + t92 + (-pkin(4) - pkin(5)) * t42 - t152 + (-pkin(4) * t88 + t147 + t53) * t72 + t112) + t43 * mrSges(7,2) - t42 * mrSges(7,1) + t72 * t54 - t71 * t50;
t110 = -t43 * mrSges(6,3) - t72 * t56 - t122 + m(6) * (qJD(5) * t148 + t152 + (t72 * t88 + t42) * pkin(4) + t111) + t42 * mrSges(6,1) + t71 * t52;
t108 = m(5) * t111 + t42 * mrSges(5,1) + t43 * mrSges(5,2) + t71 * t51 + t72 * t55 + t110;
t107 = -m(4) * (t92 + t121) + t108;
t6 = t107 + t142 * t77 + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) - t135 * qJD(2) + m(3) * t137 + t136 * t133;
t146 = t100 * t6 + t103 * t4;
t114 = t70 * mrSges(7,2) + t88 * t54 + t129;
t2 = m(2) * t127 + qJDD(1) * mrSges(2,1) - t106 * mrSges(2,2) - t153;
t1 = m(2) * t120 - t106 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t100 * t4 + t103 * t6;
t3 = [-m(1) * g(1) + t104 * t1 - t101 * t2, t1, t6, t77 * mrSges(4,2) - t82 * t132 - t118, t9, -t42 * mrSges(6,2) - t71 * t47 + t114 + t151, t114; -m(1) * g(2) + t101 * t1 + t104 * t2, t2, t4, -t77 * mrSges(4,1) - qJDD(2) * mrSges(4,3) - qJD(2) * t82 - t74 * t133 - t107, t8, t110, -t43 * mrSges(7,3) - t72 * t48 - t128; (-m(1) - m(2)) * g(3) + t146, -m(2) * g(3) + t146, t153, t76 * mrSges(4,1) + qJDD(2) * mrSges(4,2) + qJD(2) * t81 + t74 * t132 - t115, t108, -t70 * mrSges(6,1) - t88 * t52 + (t47 - t48) * t72 + (mrSges(6,2) - mrSges(7,3)) * t43 + t119, t122;];
f_new  = t3;
