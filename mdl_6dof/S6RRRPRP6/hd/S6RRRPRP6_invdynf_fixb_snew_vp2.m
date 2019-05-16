% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 08:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:59:52
% EndTime: 2019-05-07 08:00:04
% DurationCPUTime: 4.91s
% Computational Cost: add. (86378->212), mult. (190683->282), div. (0->0), fcn. (152232->12), ass. (0->106)
t104 = sin(pkin(6));
t109 = sin(qJ(2));
t113 = cos(qJ(2));
t133 = qJD(1) * qJD(2);
t95 = (-qJDD(1) * t113 + t109 * t133) * t104;
t103 = sin(pkin(11));
t105 = cos(pkin(11));
t108 = sin(qJ(3));
t112 = cos(qJ(3));
t135 = qJD(1) * t113;
t106 = cos(pkin(6));
t115 = qJD(1) ^ 2;
t110 = sin(qJ(1));
t114 = cos(qJ(1));
t126 = t110 * g(1) - g(2) * t114;
t147 = pkin(8) * t104;
t90 = qJDD(1) * pkin(1) + t115 * t147 + t126;
t140 = t106 * t90;
t122 = -g(1) * t114 - g(2) * t110;
t91 = -pkin(1) * t115 + qJDD(1) * t147 + t122;
t141 = t109 * t140 + t113 * t91;
t136 = qJD(1) * t104;
t93 = (-pkin(2) * t113 - pkin(9) * t109) * t136;
t100 = qJD(1) * t106 + qJD(2);
t98 = t100 ^ 2;
t99 = qJDD(1) * t106 + qJDD(2);
t61 = -pkin(2) * t98 + pkin(9) * t99 + (-g(3) * t109 + t93 * t135) * t104 + t141;
t146 = g(3) * t106;
t94 = (qJDD(1) * t109 + t113 * t133) * t104;
t62 = pkin(2) * t95 - pkin(9) * t94 - t146 + (-t90 + (pkin(2) * t109 - pkin(9) * t113) * t100 * qJD(1)) * t104;
t123 = -t108 * t61 + t112 * t62;
t128 = t109 * t136;
t83 = t100 * t112 - t108 * t128;
t70 = qJD(3) * t83 + t108 * t99 + t112 * t94;
t84 = t100 * t108 + t112 * t128;
t87 = qJDD(3) + t95;
t127 = t104 * t135;
t97 = qJD(3) - t127;
t28 = (t83 * t97 - t70) * qJ(4) + (t83 * t84 + t87) * pkin(3) + t123;
t142 = t108 * t62 + t112 * t61;
t69 = -qJD(3) * t84 - t108 * t94 + t112 * t99;
t78 = pkin(3) * t97 - qJ(4) * t84;
t82 = t83 ^ 2;
t31 = -pkin(3) * t82 + qJ(4) * t69 - t78 * t97 + t142;
t75 = t103 * t83 + t105 * t84;
t149 = -0.2e1 * qJD(4) * t75 - t103 * t31 + t105 * t28;
t74 = -t103 * t84 + t105 * t83;
t56 = -pkin(4) * t74 - pkin(10) * t75;
t96 = t97 ^ 2;
t22 = -pkin(4) * t87 - pkin(10) * t96 + t75 * t56 - t149;
t107 = sin(qJ(5));
t111 = cos(qJ(5));
t53 = t103 * t69 + t105 * t70;
t65 = t107 * t97 + t111 * t75;
t36 = -t65 * qJD(5) - t107 * t53 + t111 * t87;
t64 = -t107 * t75 + t111 * t97;
t37 = t64 * qJD(5) + t107 * t87 + t111 * t53;
t73 = qJD(5) - t74;
t46 = pkin(5) * t73 - t65 * qJ(6);
t47 = mrSges(7,1) * t73 - t65 * mrSges(7,3);
t63 = t64 ^ 2;
t130 = m(7) * (-t36 * pkin(5) - t63 * qJ(6) + t65 * t46 + qJDD(6) + t22) + t37 * mrSges(7,2) + t65 * t47;
t44 = -mrSges(7,2) * t73 + t64 * mrSges(7,3);
t45 = -mrSges(6,2) * t73 + t64 * mrSges(6,3);
t48 = mrSges(6,1) * t73 - t65 * mrSges(6,3);
t148 = -m(6) * t22 - t37 * mrSges(6,2) + (t45 + t44) * t64 + (mrSges(6,1) + mrSges(7,1)) * t36 - t65 * t48 - t130;
t129 = 0.2e1 * qJD(4) * t74 + t103 * t28 + t105 * t31;
t23 = -pkin(4) * t96 + pkin(10) * t87 + t56 * t74 + t129;
t137 = t104 * t113;
t121 = -g(3) * t137 - t109 * t91 + t113 * t140;
t60 = -pkin(2) * t99 - pkin(9) * t98 + t93 * t128 - t121;
t117 = -pkin(3) * t69 - qJ(4) * t82 + t84 * t78 + qJDD(4) + t60;
t52 = -t103 * t70 + t105 * t69;
t26 = (-t74 * t97 - t53) * pkin(10) + (t75 * t97 - t52) * pkin(4) + t117;
t144 = t107 * t26 + t111 * t23;
t138 = t104 * t109;
t124 = -t107 * t23 + t111 * t26;
t51 = qJDD(5) - t52;
t132 = m(7) * (-0.2e1 * qJD(6) * t65 + (t64 * t73 - t37) * qJ(6) + (t64 * t65 + t51) * pkin(5) + t124) + t73 * t44 + t51 * mrSges(7,1);
t42 = -mrSges(7,1) * t64 + mrSges(7,2) * t65;
t43 = -mrSges(6,1) * t64 + mrSges(6,2) * t65;
t13 = m(6) * t124 + t51 * mrSges(6,1) + t73 * t45 + (-t43 - t42) * t65 + (-mrSges(6,3) - mrSges(7,3)) * t37 + t132;
t131 = m(7) * (-t63 * pkin(5) + t36 * qJ(6) + 0.2e1 * qJD(6) * t64 - t46 * t73 + t144) + t36 * mrSges(7,3) + t64 * t42;
t16 = m(6) * t144 + t36 * mrSges(6,3) + t64 * t43 + (-t48 - t47) * t73 + (-mrSges(6,2) - mrSges(7,2)) * t51 + t131;
t66 = -mrSges(5,2) * t97 + mrSges(5,3) * t74;
t67 = mrSges(5,1) * t97 - mrSges(5,3) * t75;
t119 = -m(5) * t117 + t52 * mrSges(5,1) - t53 * mrSges(5,2) - t107 * t16 - t111 * t13 + t74 * t66 - t75 * t67;
t77 = -mrSges(4,2) * t97 + mrSges(4,3) * t83;
t79 = mrSges(4,1) * t97 - mrSges(4,3) * t84;
t116 = m(4) * t60 - t69 * mrSges(4,1) + t70 * mrSges(4,2) - t83 * t77 + t84 * t79 - t119;
t89 = -mrSges(3,2) * t100 + mrSges(3,3) * t127;
t92 = (-mrSges(3,1) * t113 + mrSges(3,2) * t109) * t136;
t10 = m(3) * t121 + t99 * mrSges(3,1) - t94 * mrSges(3,3) + t100 * t89 - t92 * t128 - t116;
t55 = -mrSges(5,1) * t74 + mrSges(5,2) * t75;
t11 = m(5) * t129 - t87 * mrSges(5,2) + t52 * mrSges(5,3) - t107 * t13 + t111 * t16 + t74 * t55 - t97 * t67;
t14 = m(5) * t149 + t87 * mrSges(5,1) - t53 * mrSges(5,3) - t75 * t55 + t97 * t66 + t148;
t76 = -mrSges(4,1) * t83 + mrSges(4,2) * t84;
t7 = m(4) * t123 + t87 * mrSges(4,1) - t70 * mrSges(4,3) + t103 * t11 + t105 * t14 - t84 * t76 + t97 * t77;
t8 = m(4) * t142 - t87 * mrSges(4,2) + t69 * mrSges(4,3) - t103 * t14 + t105 * t11 + t83 * t76 - t97 * t79;
t88 = mrSges(3,1) * t100 - mrSges(3,3) * t128;
t4 = m(3) * (-g(3) * t138 + t141) - t95 * mrSges(3,3) - t99 * mrSges(3,2) + t92 * t127 - t100 * t88 + t112 * t8 - t108 * t7;
t6 = m(3) * (-t104 * t90 - t146) + t94 * mrSges(3,2) + t95 * mrSges(3,1) + t108 * t8 + t112 * t7 + (t109 * t88 - t113 * t89) * t136;
t134 = t10 * t137 + t106 * t6 + t4 * t138;
t2 = m(2) * t122 - t115 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t109 * t10 + t113 * t4;
t1 = m(2) * t126 + qJDD(1) * mrSges(2,1) - t115 * mrSges(2,2) - t104 * t6 + (t10 * t113 + t109 * t4) * t106;
t3 = [-m(1) * g(1) - t1 * t110 + t114 * t2, t2, t4, t8, t11, t16, -t51 * mrSges(7,2) - t73 * t47 + t131; -m(1) * g(2) + t1 * t114 + t110 * t2, t1, t10, t7, t14, t13, -t37 * mrSges(7,3) - t65 * t42 + t132; (-m(1) - m(2)) * g(3) + t134, -m(2) * g(3) + t134, t6, t116, -t119, -t148, -t36 * mrSges(7,1) - t64 * t44 + t130;];
f_new  = t3;
