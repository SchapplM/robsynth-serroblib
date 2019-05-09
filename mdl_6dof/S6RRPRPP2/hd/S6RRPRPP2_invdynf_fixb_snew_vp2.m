% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPP2
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
% Datum: 2019-05-06 12:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:22:57
% EndTime: 2019-05-06 12:23:01
% DurationCPUTime: 1.68s
% Computational Cost: add. (17288->207), mult. (39069->247), div. (0->0), fcn. (26344->8), ass. (0->94)
t102 = sin(qJ(4));
t143 = cos(qJ(4));
t101 = sin(pkin(9));
t103 = sin(qJ(2));
t105 = cos(qJ(2));
t133 = cos(pkin(9));
t85 = (t101 * t105 + t103 * t133) * qJD(1);
t73 = -t143 * qJD(2) + t102 * t85;
t130 = qJD(1) * t105;
t131 = qJD(1) * t103;
t84 = -t101 * t131 + t133 * t130;
t83 = qJD(4) - t84;
t141 = t73 * t83;
t129 = qJD(1) * qJD(2);
t91 = t103 * qJDD(1) + t105 * t129;
t92 = t105 * qJDD(1) - t103 * t129;
t71 = t101 * t92 + t133 * t91;
t39 = -t73 * qJD(4) + t102 * qJDD(2) + t143 * t71;
t150 = (-t39 + t141) * qJ(5);
t107 = qJD(2) ^ 2;
t108 = qJD(1) ^ 2;
t104 = sin(qJ(1));
t106 = cos(qJ(1));
t120 = -t106 * g(1) - t104 * g(2);
t88 = -t108 * pkin(1) + qJDD(1) * pkin(7) + t120;
t134 = t103 * t88;
t142 = pkin(2) * t108;
t46 = qJDD(2) * pkin(2) - t91 * qJ(3) - t134 + (qJ(3) * t129 + t103 * t142 - g(3)) * t105;
t100 = t105 ^ 2;
t124 = -t103 * g(3) + t105 * t88;
t93 = qJD(2) * pkin(2) - qJ(3) * t131;
t47 = t92 * qJ(3) - qJD(2) * t93 - t100 * t142 + t124;
t126 = 0.2e1 * qJD(3) * t84 + t101 * t46 + t133 * t47;
t64 = -t84 * pkin(3) - t85 * pkin(8);
t25 = -t107 * pkin(3) + qJDD(2) * pkin(8) + t84 * t64 + t126;
t125 = t104 * g(1) - t106 * g(2);
t117 = -qJDD(1) * pkin(1) - t125;
t110 = -t92 * pkin(2) + qJDD(3) + t93 * t131 + (-qJ(3) * t100 - pkin(7)) * t108 + t117;
t70 = -t101 * t91 + t133 * t92;
t27 = (-qJD(2) * t84 - t71) * pkin(8) + (qJD(2) * t85 - t70) * pkin(3) + t110;
t138 = t102 * t27 + t143 * t25;
t145 = 2 * qJD(5);
t74 = t102 * qJD(2) + t143 * t85;
t48 = t73 * pkin(4) - t74 * qJ(5);
t69 = qJDD(4) - t70;
t82 = t83 ^ 2;
t116 = -t82 * pkin(4) + t69 * qJ(5) + t83 * t145 - t73 * t48 + t138;
t61 = -t83 * mrSges(6,1) + t74 * mrSges(6,2);
t149 = m(6) * t116 + t69 * mrSges(6,3) + t83 * t61;
t146 = -0.2e1 * t74;
t121 = -t102 * t25 + t143 * t27;
t19 = -t69 * pkin(4) - t82 * qJ(5) + t74 * t48 + qJDD(5) - t121;
t56 = t83 * mrSges(7,2) + t73 * mrSges(7,3);
t127 = t83 * t56 + t69 * mrSges(7,1) - m(7) * (qJD(6) * t146 + (-t39 - t141) * qJ(6) + (t73 * t74 - t69) * pkin(5) + t19);
t119 = m(6) * t19 - t127;
t49 = t73 * mrSges(6,1) - t74 * mrSges(6,3);
t136 = -t73 * mrSges(5,1) - t74 * mrSges(5,2) - t49;
t139 = -mrSges(5,3) - mrSges(6,2);
t50 = -t73 * mrSges(7,1) + t74 * mrSges(7,2);
t55 = -t73 * mrSges(6,2) + t83 * mrSges(6,3);
t57 = -t83 * mrSges(5,2) - t73 * mrSges(5,3);
t11 = m(5) * t121 + (t57 + t55) * t83 + (mrSges(5,1) + mrSges(6,1)) * t69 + (t50 + t136) * t74 + (mrSges(7,3) + t139) * t39 - t119;
t38 = t74 * qJD(4) - t143 * qJDD(2) + t102 * t71;
t58 = -t83 * pkin(5) - t74 * qJ(6);
t72 = t73 ^ 2;
t128 = m(7) * (-t72 * pkin(5) + t38 * qJ(6) + 0.2e1 * qJD(6) * t73 + t83 * t58 + t116) + t73 * t50 + t38 * mrSges(7,3);
t59 = -t83 * mrSges(7,1) - t74 * mrSges(7,3);
t60 = t83 * mrSges(5,1) - t74 * mrSges(5,3);
t12 = m(5) * t138 + (-t60 + t59) * t83 + t136 * t73 + (-mrSges(5,2) + mrSges(7,2)) * t69 + t139 * t38 + t128 + t149;
t75 = -qJD(2) * mrSges(4,2) + t84 * mrSges(4,3);
t76 = qJD(2) * mrSges(4,1) - t85 * mrSges(4,3);
t112 = -m(4) * t110 + t70 * mrSges(4,1) - t71 * mrSges(4,2) - t102 * t12 - t143 * t11 + t84 * t75 - t85 * t76;
t94 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t131;
t95 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t130;
t148 = (t103 * t94 - t105 * t95) * qJD(1) + m(3) * (-t108 * pkin(7) + t117) - t92 * mrSges(3,1) + t91 * mrSges(3,2) - t112;
t137 = -t101 * t47 + t133 * t46;
t113 = qJDD(2) * pkin(3) + t107 * pkin(8) - t85 * t64 + t137;
t132 = qJD(3) * t85;
t79 = -0.2e1 * t132;
t122 = m(7) * (-t72 * qJ(6) + qJDD(6) + t79 + (-pkin(4) - pkin(5)) * t38 - t150 + (-pkin(4) * t83 + t145 + t58) * t74 + t113) + t39 * mrSges(7,2) - t38 * mrSges(7,1) + t74 * t59 - t73 * t56;
t24 = -t113 + 0.2e1 * t132;
t115 = m(6) * (qJD(5) * t146 + t150 + (t74 * t83 + t38) * pkin(4) + t24) + t73 * t55 + t38 * mrSges(6,1) - t122;
t147 = m(5) * t24 + t38 * mrSges(5,1) + (t60 - t61) * t74 + (mrSges(5,2) - mrSges(6,3)) * t39 + t73 * t57 + t115;
t63 = -t84 * mrSges(4,1) + t85 * mrSges(4,2);
t7 = m(4) * t126 - qJDD(2) * mrSges(4,2) + t70 * mrSges(4,3) - qJD(2) * t76 - t102 * t11 + t143 * t12 + t84 * t63;
t8 = qJDD(2) * mrSges(4,1) - t71 * mrSges(4,3) + qJD(2) * t75 + m(4) * (t79 + t137) - t85 * t63 - t147;
t90 = (-mrSges(3,1) * t105 + mrSges(3,2) * t103) * qJD(1);
t4 = m(3) * (-t105 * g(3) - t134) - t91 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t90 * t131 + qJD(2) * t95 + t101 * t7 + t133 * t8;
t5 = m(3) * t124 - qJDD(2) * mrSges(3,2) + t92 * mrSges(3,3) - qJD(2) * t94 - t101 * t8 + t130 * t90 + t133 * t7;
t144 = t103 * t5 + t105 * t4;
t114 = t69 * mrSges(7,2) + t83 * t59 + t128;
t6 = m(2) * t125 + qJDD(1) * mrSges(2,1) - t108 * mrSges(2,2) - t148;
t1 = m(2) * t120 - t108 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t103 * t4 + t105 * t5;
t2 = [-m(1) * g(1) + t106 * t1 - t104 * t6, t1, t5, t7, t12, -t38 * mrSges(6,2) - t73 * t49 + t114 + t149, t114; -m(1) * g(2) + t104 * t1 + t106 * t6, t6, t4, t8, t11, -t39 * mrSges(6,3) - t74 * t61 + t115, -t39 * mrSges(7,3) - t74 * t50 - t127; (-m(1) - m(2)) * g(3) + t144, -m(2) * g(3) + t144, t148, -t112, t147, -t69 * mrSges(6,1) - t83 * t55 + (t49 - t50) * t74 + (mrSges(6,2) - mrSges(7,3)) * t39 + t119, t122;];
f_new  = t2;
