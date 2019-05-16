% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 05:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRP7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:20:24
% EndTime: 2019-05-08 05:20:41
% DurationCPUTime: 5.12s
% Computational Cost: add. (94456->212), mult. (200896->279), div. (0->0), fcn. (161952->12), ass. (0->107)
t107 = sin(qJ(4));
t112 = cos(qJ(4));
t108 = sin(qJ(3));
t113 = cos(qJ(3));
t105 = cos(pkin(6));
t100 = t105 * qJDD(1) + qJDD(2);
t104 = sin(pkin(6));
t109 = sin(qJ(2));
t114 = cos(qJ(2));
t134 = qJD(1) * t114;
t116 = qJD(1) ^ 2;
t110 = sin(qJ(1));
t115 = cos(qJ(1));
t126 = t110 * g(1) - t115 * g(2);
t89 = t116 * t104 * pkin(8) + qJDD(1) * pkin(1) + t126;
t138 = t105 * t89;
t122 = -t115 * g(1) - t110 * g(2);
t132 = qJDD(1) * t104;
t90 = -t116 * pkin(1) + pkin(8) * t132 + t122;
t139 = t109 * t138 + t114 * t90;
t135 = qJD(1) * t104;
t92 = (-pkin(2) * t114 - pkin(9) * t109) * t135;
t101 = t105 * qJD(1) + qJD(2);
t99 = t101 ^ 2;
t61 = -t99 * pkin(2) + t100 * pkin(9) + (-g(3) * t109 + t92 * t134) * t104 + t139;
t145 = t105 * g(3);
t93 = (qJD(2) * t134 + qJDD(1) * t109) * t104;
t128 = t109 * t135;
t94 = -qJD(2) * t128 + t114 * t132;
t62 = -t94 * pkin(2) - t93 * pkin(9) - t145 + (-t89 + (pkin(2) * t109 - pkin(9) * t114) * t101 * qJD(1)) * t104;
t124 = -t108 * t61 + t113 * t62;
t81 = t113 * t101 - t108 * t128;
t70 = t81 * qJD(3) + t108 * t100 + t113 * t93;
t82 = t108 * t101 + t113 * t128;
t86 = qJDD(3) - t94;
t127 = t104 * t134;
t97 = qJD(3) - t127;
t28 = (t81 * t97 - t70) * pkin(10) + (t81 * t82 + t86) * pkin(3) + t124;
t140 = t108 * t62 + t113 * t61;
t69 = -t82 * qJD(3) + t113 * t100 - t108 * t93;
t77 = t97 * pkin(3) - t82 * pkin(10);
t80 = t81 ^ 2;
t31 = -t80 * pkin(3) + t69 * pkin(10) - t97 * t77 + t140;
t123 = -t107 * t31 + t112 * t28;
t72 = -t107 * t82 + t112 * t81;
t73 = t107 * t81 + t112 * t82;
t56 = -pkin(4) * t72 - pkin(11) * t73;
t85 = qJDD(4) + t86;
t96 = qJD(4) + t97;
t95 = t96 ^ 2;
t22 = -t85 * pkin(4) - t95 * pkin(11) + t73 * t56 - t123;
t106 = sin(qJ(5));
t111 = cos(qJ(5));
t45 = t72 * qJD(4) + t107 * t69 + t112 * t70;
t65 = t106 * t96 + t111 * t73;
t34 = -t65 * qJD(5) - t106 * t45 + t111 * t85;
t64 = -t106 * t73 + t111 * t96;
t35 = t64 * qJD(5) + t106 * t85 + t111 * t45;
t71 = qJD(5) - t72;
t51 = pkin(5) * t71 - t65 * qJ(6);
t52 = mrSges(7,1) * t71 - t65 * mrSges(7,3);
t63 = t64 ^ 2;
t129 = m(7) * (-t34 * pkin(5) - t63 * qJ(6) + t65 * t51 + qJDD(6) + t22) + t35 * mrSges(7,2) + t65 * t52;
t49 = -mrSges(7,2) * t71 + t64 * mrSges(7,3);
t50 = -mrSges(6,2) * t71 + t64 * mrSges(6,3);
t53 = mrSges(6,1) * t71 - t65 * mrSges(6,3);
t146 = m(6) * t22 + t35 * mrSges(6,2) - (t50 + t49) * t64 - (mrSges(6,1) + mrSges(7,1)) * t34 + t65 * t53 + t129;
t142 = t107 * t28 + t112 * t31;
t23 = -t95 * pkin(4) + t85 * pkin(11) + t72 * t56 + t142;
t136 = t104 * t114;
t121 = -g(3) * t136 - t109 * t90 + t114 * t138;
t60 = -t100 * pkin(2) - t99 * pkin(9) + t92 * t128 - t121;
t118 = -t69 * pkin(3) - t80 * pkin(10) + t82 * t77 + t60;
t44 = -t73 * qJD(4) - t107 * t70 + t112 * t69;
t26 = (-t72 * t96 - t45) * pkin(11) + (t73 * t96 - t44) * pkin(4) + t118;
t143 = t106 * t26 + t111 * t23;
t137 = t104 * t109;
t125 = -t106 * t23 + t111 * t26;
t43 = qJDD(5) - t44;
t131 = m(7) * (-0.2e1 * qJD(6) * t65 + (t64 * t71 - t35) * qJ(6) + (t64 * t65 + t43) * pkin(5) + t125) + t71 * t49 + t43 * mrSges(7,1);
t47 = -mrSges(7,1) * t64 + mrSges(7,2) * t65;
t48 = -mrSges(6,1) * t64 + mrSges(6,2) * t65;
t13 = m(6) * t125 + t43 * mrSges(6,1) + t71 * t50 + (-t48 - t47) * t65 + (-mrSges(6,3) - mrSges(7,3)) * t35 + t131;
t130 = m(7) * (-t63 * pkin(5) + t34 * qJ(6) + 0.2e1 * qJD(6) * t64 - t51 * t71 + t143) + t34 * mrSges(7,3) + t64 * t47;
t15 = m(6) * t143 + t34 * mrSges(6,3) + t64 * t48 + (-t53 - t52) * t71 + (-mrSges(6,2) - mrSges(7,2)) * t43 + t130;
t66 = -t96 * mrSges(5,2) + t72 * mrSges(5,3);
t67 = t96 * mrSges(5,1) - t73 * mrSges(5,3);
t120 = -m(5) * t118 + t44 * mrSges(5,1) - t45 * mrSges(5,2) - t106 * t15 - t111 * t13 + t72 * t66 - t73 * t67;
t75 = -t97 * mrSges(4,2) + t81 * mrSges(4,3);
t76 = t97 * mrSges(4,1) - t82 * mrSges(4,3);
t117 = m(4) * t60 - t69 * mrSges(4,1) + t70 * mrSges(4,2) - t81 * t75 + t82 * t76 - t120;
t88 = -t101 * mrSges(3,2) + mrSges(3,3) * t127;
t91 = (-mrSges(3,1) * t114 + mrSges(3,2) * t109) * t135;
t10 = m(3) * t121 + t100 * mrSges(3,1) - t93 * mrSges(3,3) + t101 * t88 - t91 * t128 - t117;
t55 = -mrSges(5,1) * t72 + mrSges(5,2) * t73;
t11 = m(5) * t142 - t85 * mrSges(5,2) + t44 * mrSges(5,3) - t106 * t13 + t111 * t15 + t72 * t55 - t96 * t67;
t16 = m(5) * t123 + t85 * mrSges(5,1) - t45 * mrSges(5,3) - t73 * t55 + t96 * t66 - t146;
t74 = -mrSges(4,1) * t81 + mrSges(4,2) * t82;
t7 = m(4) * t124 + t86 * mrSges(4,1) - t70 * mrSges(4,3) + t107 * t11 + t112 * t16 - t82 * t74 + t97 * t75;
t8 = m(4) * t140 - t86 * mrSges(4,2) + t69 * mrSges(4,3) - t107 * t16 + t112 * t11 + t81 * t74 - t97 * t76;
t87 = t101 * mrSges(3,1) - mrSges(3,3) * t128;
t4 = m(3) * (-g(3) * t137 + t139) + t94 * mrSges(3,3) - t100 * mrSges(3,2) + t91 * t127 - t101 * t87 + t113 * t8 - t108 * t7;
t6 = m(3) * (-t104 * t89 - t145) + t93 * mrSges(3,2) - t94 * mrSges(3,1) + t108 * t8 + t113 * t7 + (t109 * t87 - t114 * t88) * t135;
t133 = t10 * t136 + t105 * t6 + t4 * t137;
t2 = m(2) * t122 - t116 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t109 * t10 + t114 * t4;
t1 = m(2) * t126 + qJDD(1) * mrSges(2,1) - t116 * mrSges(2,2) - t104 * t6 + (t114 * t10 + t109 * t4) * t105;
t3 = [-m(1) * g(1) - t110 * t1 + t115 * t2, t2, t4, t8, t11, t15, -t43 * mrSges(7,2) - t71 * t52 + t130; -m(1) * g(2) + t115 * t1 + t110 * t2, t1, t10, t7, t16, t13, -t35 * mrSges(7,3) - t65 * t47 + t131; (-m(1) - m(2)) * g(3) + t133, -m(2) * g(3) + t133, t6, t117, -t120, t146, -t34 * mrSges(7,1) - t64 * t49 + t129;];
f_new  = t3;
