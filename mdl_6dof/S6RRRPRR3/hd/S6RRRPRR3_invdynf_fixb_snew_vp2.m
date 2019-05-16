% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 10:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:10:39
% EndTime: 2019-05-07 10:10:50
% DurationCPUTime: 2.89s
% Computational Cost: add. (36801->203), mult. (78600->257), div. (0->0), fcn. (55550->10), ass. (0->101)
t113 = sin(qJ(2));
t117 = cos(qJ(2));
t138 = qJD(1) * t113;
t119 = qJD(1) ^ 2;
t140 = t117 ^ 2 * t119;
t114 = sin(qJ(1));
t118 = cos(qJ(1));
t139 = t114 * g(1) - t118 * g(2);
t83 = -qJDD(1) * pkin(1) - t119 * pkin(7) - t139;
t136 = qJD(1) * qJD(2);
t90 = t117 * qJDD(1) - t113 * t136;
t93 = qJD(2) * pkin(2) - pkin(8) * t138;
t128 = -t90 * pkin(2) - pkin(8) * t140 + t93 * t138 + t83;
t106 = qJD(2) + qJD(3);
t112 = sin(qJ(3));
t137 = qJD(1) * t117;
t148 = cos(qJ(3));
t81 = t112 * t138 - t148 * t137;
t142 = t106 * t81;
t82 = (t112 * t117 + t148 * t113) * qJD(1);
t89 = t113 * qJDD(1) + t117 * t136;
t57 = t82 * qJD(3) + t112 * t89 - t148 * t90;
t58 = -t81 * qJD(3) + t112 * t90 + t148 * t89;
t125 = t57 * pkin(3) + t128 + (t142 - t58) * qJ(4);
t110 = sin(qJ(6));
t115 = cos(qJ(6));
t147 = pkin(3) * t106;
t150 = 2 * qJD(4);
t76 = -t106 * pkin(4) - t82 * pkin(9);
t80 = t81 ^ 2;
t122 = -t57 * pkin(4) - t80 * pkin(9) - t125 + (-t147 + t150 + t76) * t82;
t105 = qJDD(2) + qJDD(3);
t100 = qJDD(5) - t105;
t111 = sin(qJ(5));
t116 = cos(qJ(5));
t104 = t106 ^ 2;
t133 = -t118 * g(1) - t114 * g(2);
t84 = -t119 * pkin(1) + qJDD(1) * pkin(7) + t133;
t141 = t113 * t84;
t46 = qJDD(2) * pkin(2) - t89 * pkin(8) - t141 + (pkin(2) * t113 * t119 + pkin(8) * t136 - g(3)) * t117;
t135 = -t113 * g(3) + t117 * t84;
t47 = -pkin(2) * t140 + t90 * pkin(8) - qJD(2) * t93 + t135;
t134 = -t112 * t47 + t148 * t46;
t68 = t81 * pkin(3) - t82 * qJ(4);
t30 = -t105 * pkin(3) - t104 * qJ(4) + t82 * t68 + qJDD(4) - t134;
t22 = (-t58 - t142) * pkin(9) + (t81 * t82 - t105) * pkin(4) + t30;
t144 = t112 * t46 + t148 * t47;
t127 = -t104 * pkin(3) + t105 * qJ(4) + t106 * t150 - t81 * t68 + t144;
t24 = -t80 * pkin(4) + t57 * pkin(9) + t106 * t76 + t127;
t145 = t111 * t22 + t116 * t24;
t66 = -t111 * t82 + t116 * t81;
t67 = t111 * t81 + t116 * t82;
t42 = -t66 * pkin(5) - t67 * pkin(10);
t101 = qJD(5) - t106;
t99 = t101 ^ 2;
t17 = -t99 * pkin(5) + t100 * pkin(10) + t66 * t42 + t145;
t34 = -t67 * qJD(5) - t111 * t58 + t116 * t57;
t35 = t66 * qJD(5) + t111 * t57 + t116 * t58;
t18 = (-t101 * t66 - t35) * pkin(10) + (t101 * t67 - t34) * pkin(5) + t122;
t50 = t115 * t101 - t110 * t67;
t26 = t50 * qJD(6) + t110 * t100 + t115 * t35;
t33 = qJDD(6) - t34;
t51 = t110 * t101 + t115 * t67;
t36 = -t50 * mrSges(7,1) + t51 * mrSges(7,2);
t62 = qJD(6) - t66;
t37 = -t62 * mrSges(7,2) + t50 * mrSges(7,3);
t14 = m(7) * (-t110 * t17 + t115 * t18) - t26 * mrSges(7,3) + t33 * mrSges(7,1) - t51 * t36 + t62 * t37;
t25 = -t51 * qJD(6) + t115 * t100 - t110 * t35;
t38 = t62 * mrSges(7,1) - t51 * mrSges(7,3);
t15 = m(7) * (t110 * t18 + t115 * t17) + t25 * mrSges(7,3) - t33 * mrSges(7,2) + t50 * t36 - t62 * t38;
t60 = -t101 * mrSges(6,2) + t66 * mrSges(6,3);
t61 = t101 * mrSges(6,1) - t67 * mrSges(6,3);
t129 = m(6) * t122 - t34 * mrSges(6,1) + t35 * mrSges(6,2) + t110 * t15 + t115 * t14 - t66 * t60 + t67 * t61;
t74 = -t106 * mrSges(5,1) + t82 * mrSges(5,2);
t75 = -t81 * mrSges(5,2) + t106 * mrSges(5,3);
t123 = t58 * mrSges(5,3) + t82 * t74 + t129 - m(5) * ((-(2 * qJD(4)) + t147) * t82 + t125) - t57 * mrSges(5,1) - t81 * t75;
t72 = -t106 * mrSges(4,2) - t81 * mrSges(4,3);
t73 = t106 * mrSges(4,1) - t82 * mrSges(4,3);
t121 = m(4) * t128 + t57 * mrSges(4,1) + t58 * mrSges(4,2) + t81 * t72 + t82 * t73 - t123;
t91 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t138;
t92 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t137;
t151 = m(3) * t83 - t90 * mrSges(3,1) + t89 * mrSges(3,2) + (t113 * t91 - t117 * t92) * qJD(1) + t121;
t41 = -t66 * mrSges(6,1) + t67 * mrSges(6,2);
t10 = m(6) * t145 - t100 * mrSges(6,2) + t34 * mrSges(6,3) - t101 * t61 - t110 * t14 + t115 * t15 + t66 * t41;
t132 = -t111 * t24 + t116 * t22;
t124 = m(7) * (-t100 * pkin(5) - t99 * pkin(10) + t67 * t42 - t132) - t25 * mrSges(7,1) + t26 * mrSges(7,2) - t50 * t37 + t51 * t38;
t11 = m(6) * t132 + t100 * mrSges(6,1) - t35 * mrSges(6,3) + t101 * t60 - t67 * t41 - t124;
t130 = m(5) * t127 + t105 * mrSges(5,3) + t116 * t10 + t106 * t74 - t111 * t11;
t69 = t81 * mrSges(5,1) - t82 * mrSges(5,3);
t143 = -t81 * mrSges(4,1) - t82 * mrSges(4,2) - t69;
t146 = -mrSges(4,3) - mrSges(5,2);
t6 = m(4) * t144 - t105 * mrSges(4,2) - t106 * t73 + t143 * t81 + t146 * t57 + t130;
t126 = -m(5) * t30 - t111 * t10 - t116 * t11;
t7 = m(4) * t134 + t143 * t82 + t146 * t58 + (t72 + t75) * t106 + (mrSges(4,1) + mrSges(5,1)) * t105 + t126;
t88 = (-mrSges(3,1) * t117 + mrSges(3,2) * t113) * qJD(1);
t4 = m(3) * (-t117 * g(3) - t141) - t89 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t88 * t138 + qJD(2) * t92 + t112 * t6 + t148 * t7;
t5 = m(3) * t135 - qJDD(2) * mrSges(3,2) + t90 * mrSges(3,3) - qJD(2) * t91 - t112 * t7 + t88 * t137 + t148 * t6;
t149 = t113 * t5 + t117 * t4;
t8 = m(2) * t139 + qJDD(1) * mrSges(2,1) - t119 * mrSges(2,2) - t151;
t1 = m(2) * t133 - t119 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t113 * t4 + t117 * t5;
t2 = [-m(1) * g(1) + t118 * t1 - t114 * t8, t1, t5, t6, -t57 * mrSges(5,2) - t81 * t69 + t130, t10, t15; -m(1) * g(2) + t114 * t1 + t118 * t8, t8, t4, t7, -t123, t11, t14; (-m(1) - m(2)) * g(3) + t149, -m(2) * g(3) + t149, t151, t121, -t105 * mrSges(5,1) + t58 * mrSges(5,2) - t106 * t75 + t82 * t69 - t126, t129, t124;];
f_new  = t2;
