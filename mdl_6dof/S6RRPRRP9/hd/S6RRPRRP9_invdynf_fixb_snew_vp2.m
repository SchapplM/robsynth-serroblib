% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 18:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRP9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:27:30
% EndTime: 2019-05-06 18:27:41
% DurationCPUTime: 4.88s
% Computational Cost: add. (82270->211), mult. (186760->281), div. (0->0), fcn. (149617->12), ass. (0->105)
t109 = sin(qJ(4));
t113 = cos(qJ(4));
t104 = sin(pkin(11));
t106 = cos(pkin(11));
t107 = cos(pkin(6));
t100 = t107 * qJDD(1) + qJDD(2);
t105 = sin(pkin(6));
t110 = sin(qJ(2));
t114 = cos(qJ(2));
t135 = qJD(1) * t114;
t116 = qJD(1) ^ 2;
t111 = sin(qJ(1));
t115 = cos(qJ(1));
t126 = t111 * g(1) - t115 * g(2);
t90 = t116 * t105 * pkin(8) + qJDD(1) * pkin(1) + t126;
t139 = t107 * t90;
t122 = -t115 * g(1) - t111 * g(2);
t133 = qJDD(1) * t105;
t91 = -t116 * pkin(1) + pkin(8) * t133 + t122;
t140 = t110 * t139 + t114 * t91;
t136 = qJD(1) * t105;
t92 = (-pkin(2) * t114 - qJ(3) * t110) * t136;
t101 = t107 * qJD(1) + qJD(2);
t99 = t101 ^ 2;
t61 = -t99 * pkin(2) + t100 * qJ(3) + (-g(3) * t110 + t92 * t135) * t105 + t140;
t145 = t107 * g(3);
t94 = (qJD(2) * t135 + qJDD(1) * t110) * t105;
t128 = t110 * t136;
t95 = -qJD(2) * t128 + t114 * t133;
t62 = -t95 * pkin(2) - t145 - t94 * qJ(3) + (-t90 + (pkin(2) * t110 - qJ(3) * t114) * t101 * qJD(1)) * t105;
t84 = t104 * t101 + t106 * t128;
t123 = -0.2e1 * qJD(3) * t84 - t104 * t61 + t106 * t62;
t127 = t105 * t135;
t76 = t104 * t100 + t106 * t94;
t83 = t106 * t101 - t104 * t128;
t28 = (-t83 * t127 - t76) * pkin(9) + (t83 * t84 - t95) * pkin(3) + t123;
t129 = 0.2e1 * qJD(3) * t83 + t104 * t62 + t106 * t61;
t75 = t106 * t100 - t104 * t94;
t77 = -pkin(3) * t127 - t84 * pkin(9);
t82 = t83 ^ 2;
t31 = -t82 * pkin(3) + t75 * pkin(9) + t77 * t127 + t129;
t124 = -t109 * t31 + t113 * t28;
t70 = -t109 * t84 + t113 * t83;
t71 = t109 * t83 + t113 * t84;
t56 = -t70 * pkin(4) - t71 * pkin(10);
t87 = qJDD(4) - t95;
t97 = qJD(4) - t127;
t96 = t97 ^ 2;
t22 = -t87 * pkin(4) - t96 * pkin(10) + t71 * t56 - t124;
t108 = sin(qJ(5));
t112 = cos(qJ(5));
t48 = t70 * qJD(4) + t109 * t75 + t113 * t76;
t65 = t108 * t97 + t112 * t71;
t34 = -t65 * qJD(5) - t108 * t48 + t112 * t87;
t64 = -t108 * t71 + t112 * t97;
t35 = t64 * qJD(5) + t108 * t87 + t112 * t48;
t69 = qJD(5) - t70;
t51 = t69 * pkin(5) - t65 * qJ(6);
t52 = t69 * mrSges(7,1) - t65 * mrSges(7,3);
t63 = t64 ^ 2;
t130 = m(7) * (-t34 * pkin(5) - t63 * qJ(6) + t65 * t51 + qJDD(6) + t22) + t35 * mrSges(7,2) + t65 * t52;
t49 = -t69 * mrSges(7,2) + t64 * mrSges(7,3);
t50 = -t69 * mrSges(6,2) + t64 * mrSges(6,3);
t53 = t69 * mrSges(6,1) - t65 * mrSges(6,3);
t146 = m(6) * t22 + t35 * mrSges(6,2) - (t50 + t49) * t64 - (mrSges(6,1) + mrSges(7,1)) * t34 + t65 * t53 + t130;
t142 = t109 * t28 + t113 * t31;
t23 = -t96 * pkin(4) + t87 * pkin(10) + t70 * t56 + t142;
t137 = t105 * t114;
t121 = -g(3) * t137 - t110 * t91 + t114 * t139;
t60 = -t100 * pkin(2) - t99 * qJ(3) + t92 * t128 + qJDD(3) - t121;
t118 = -t75 * pkin(3) - t82 * pkin(9) + t84 * t77 + t60;
t47 = -t71 * qJD(4) - t109 * t76 + t113 * t75;
t26 = (-t70 * t97 - t48) * pkin(10) + (t71 * t97 - t47) * pkin(4) + t118;
t143 = t108 * t26 + t112 * t23;
t138 = t105 * t110;
t125 = -t108 * t23 + t112 * t26;
t46 = qJDD(5) - t47;
t132 = m(7) * (-0.2e1 * qJD(6) * t65 + (t64 * t69 - t35) * qJ(6) + (t64 * t65 + t46) * pkin(5) + t125) + t69 * t49 + t46 * mrSges(7,1);
t42 = -t64 * mrSges(7,1) + t65 * mrSges(7,2);
t43 = -t64 * mrSges(6,1) + t65 * mrSges(6,2);
t13 = m(6) * t125 + t46 * mrSges(6,1) + t69 * t50 + (-t43 - t42) * t65 + (-mrSges(6,3) - mrSges(7,3)) * t35 + t132;
t131 = m(7) * (-t63 * pkin(5) + t34 * qJ(6) + 0.2e1 * qJD(6) * t64 - t69 * t51 + t143) + t34 * mrSges(7,3) + t64 * t42;
t16 = m(6) * t143 + t34 * mrSges(6,3) + t64 * t43 + (-t53 - t52) * t69 + (-mrSges(6,2) - mrSges(7,2)) * t46 + t131;
t66 = -t97 * mrSges(5,2) + t70 * mrSges(5,3);
t67 = t97 * mrSges(5,1) - t71 * mrSges(5,3);
t120 = -m(5) * t118 + t47 * mrSges(5,1) - t48 * mrSges(5,2) - t108 * t16 - t112 * t13 + t70 * t66 - t71 * t67;
t73 = mrSges(4,2) * t127 + t83 * mrSges(4,3);
t74 = -mrSges(4,1) * t127 - t84 * mrSges(4,3);
t117 = m(4) * t60 - t75 * mrSges(4,1) + t76 * mrSges(4,2) - t83 * t73 + t84 * t74 - t120;
t89 = -t101 * mrSges(3,2) + mrSges(3,3) * t127;
t93 = (-mrSges(3,1) * t114 + mrSges(3,2) * t110) * t136;
t10 = m(3) * t121 + t100 * mrSges(3,1) - t94 * mrSges(3,3) + t101 * t89 - t93 * t128 - t117;
t55 = -t70 * mrSges(5,1) + t71 * mrSges(5,2);
t11 = m(5) * t142 - t87 * mrSges(5,2) + t47 * mrSges(5,3) - t108 * t13 + t112 * t16 + t70 * t55 - t97 * t67;
t14 = m(5) * t124 + t87 * mrSges(5,1) - t48 * mrSges(5,3) - t71 * t55 + t97 * t66 - t146;
t72 = -t83 * mrSges(4,1) + t84 * mrSges(4,2);
t7 = m(4) * t123 - t95 * mrSges(4,1) - t76 * mrSges(4,3) + t109 * t11 + t113 * t14 - t73 * t127 - t84 * t72;
t8 = m(4) * t129 + t95 * mrSges(4,2) + t75 * mrSges(4,3) - t109 * t14 + t113 * t11 + t74 * t127 + t83 * t72;
t88 = t101 * mrSges(3,1) - mrSges(3,3) * t128;
t4 = m(3) * (-g(3) * t138 + t140) + t95 * mrSges(3,3) - t100 * mrSges(3,2) + t93 * t127 - t101 * t88 + t106 * t8 - t104 * t7;
t6 = m(3) * (-t105 * t90 - t145) + t94 * mrSges(3,2) - t95 * mrSges(3,1) + t104 * t8 + t106 * t7 + (t110 * t88 - t114 * t89) * t136;
t134 = t10 * t137 + t107 * t6 + t4 * t138;
t2 = m(2) * t122 - t116 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t110 * t10 + t114 * t4;
t1 = m(2) * t126 + qJDD(1) * mrSges(2,1) - t116 * mrSges(2,2) - t105 * t6 + (t114 * t10 + t110 * t4) * t107;
t3 = [-m(1) * g(1) - t111 * t1 + t115 * t2, t2, t4, t8, t11, t16, -t46 * mrSges(7,2) - t69 * t52 + t131; -m(1) * g(2) + t115 * t1 + t111 * t2, t1, t10, t7, t14, t13, -t35 * mrSges(7,3) - t65 * t42 + t132; (-m(1) - m(2)) * g(3) + t134, -m(2) * g(3) + t134, t6, t117, -t120, t146, -t34 * mrSges(7,1) - t64 * t49 + t130;];
f_new  = t3;
