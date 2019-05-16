% Calculate vector of cutting forces with Newton-Euler
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-05-05 06:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:00:37
% EndTime: 2019-05-05 06:00:42
% DurationCPUTime: 2.11s
% Computational Cost: add. (25261->178), mult. (50131->228), div. (0->0), fcn. (31742->12), ass. (0->99)
t144 = -2 * qJD(4);
t101 = cos(qJ(3));
t102 = cos(qJ(2));
t91 = sin(pkin(11));
t93 = cos(pkin(11));
t72 = t91 * g(1) - t93 * g(2);
t94 = cos(pkin(6));
t134 = t72 * t94;
t73 = -t93 * g(1) - t91 * g(2);
t90 = -g(3) + qJDD(1);
t92 = sin(pkin(6));
t98 = sin(qJ(2));
t141 = (t90 * t92 + t134) * t102 - t98 * t73;
t110 = -qJDD(2) * pkin(2) - t141;
t100 = cos(qJ(5));
t123 = qJD(2) * qJD(3);
t117 = t101 * t123;
t97 = sin(qJ(3));
t119 = t97 * t123;
t69 = t97 * qJDD(2) + t117;
t85 = t97 * qJD(2);
t105 = pkin(3) * t119 + t85 * t144 + (-t69 - t117) * qJ(4) + t110;
t103 = qJD(3) ^ 2;
t104 = qJD(2) ^ 2;
t133 = t92 * t98;
t121 = t102 * t73 + t90 * t133 + t98 * t134;
t39 = -t104 * pkin(2) + qJDD(2) * pkin(8) + t121;
t36 = t97 * t39;
t66 = (-pkin(3) * t101 - qJ(4) * t97) * qJD(2);
t115 = -t103 * qJ(4) + t66 * t85 + qJDD(4) + t36;
t137 = pkin(9) * t104;
t138 = -pkin(3) - pkin(9);
t53 = -t92 * t72 + t94 * t90;
t24 = t69 * pkin(4) + t138 * qJDD(3) + (-pkin(4) * t123 - t97 * t137 - t53) * t101 + t115;
t70 = t101 * qJDD(2) - t119;
t79 = pkin(4) * t85 - qJD(3) * pkin(9);
t89 = t101 ^ 2;
t26 = -t79 * t85 + t138 * t70 + (-pkin(4) * t89 - pkin(8)) * t104 + t105;
t96 = sin(qJ(5));
t118 = t100 * t24 - t96 * t26;
t124 = qJD(2) * t101;
t64 = -t96 * qJD(3) - t100 * t124;
t45 = t64 * qJD(5) + t100 * qJDD(3) - t96 * t70;
t61 = qJDD(5) + t69;
t65 = t100 * qJD(3) - t96 * t124;
t82 = t85 + qJD(5);
t16 = (t64 * t82 - t45) * pkin(10) + (t64 * t65 + t61) * pkin(5) + t118;
t130 = t100 * t26 + t96 * t24;
t44 = -t65 * qJD(5) - t96 * qJDD(3) - t100 * t70;
t52 = t82 * pkin(5) - t65 * pkin(10);
t60 = t64 ^ 2;
t17 = -t60 * pkin(5) + t44 * pkin(10) - t82 * t52 + t130;
t95 = sin(qJ(6));
t99 = cos(qJ(6));
t46 = t99 * t64 - t95 * t65;
t29 = t46 * qJD(6) + t95 * t44 + t99 * t45;
t47 = t95 * t64 + t99 * t65;
t35 = -t46 * mrSges(7,1) + t47 * mrSges(7,2);
t80 = qJD(6) + t82;
t40 = -t80 * mrSges(7,2) + t46 * mrSges(7,3);
t57 = qJDD(6) + t61;
t14 = m(7) * (t99 * t16 - t95 * t17) - t29 * mrSges(7,3) + t57 * mrSges(7,1) - t47 * t35 + t80 * t40;
t28 = -t47 * qJD(6) + t99 * t44 - t95 * t45;
t41 = t80 * mrSges(7,1) - t47 * mrSges(7,3);
t15 = m(7) * (t95 * t16 + t99 * t17) + t28 * mrSges(7,3) - t57 * mrSges(7,2) + t46 * t35 - t80 * t41;
t48 = -t64 * mrSges(6,1) + t65 * mrSges(6,2);
t50 = -t82 * mrSges(6,2) + t64 * mrSges(6,3);
t11 = m(6) * t118 + t61 * mrSges(6,1) - t45 * mrSges(6,3) + t99 * t14 + t95 * t15 - t65 * t48 + t82 * t50;
t51 = t82 * mrSges(6,1) - t65 * mrSges(6,3);
t12 = m(6) * t130 - t61 * mrSges(6,2) + t44 * mrSges(6,3) - t95 * t14 + t99 * t15 + t64 * t48 - t82 * t51;
t135 = t104 * pkin(8);
t76 = -mrSges(5,1) * t124 - qJD(3) * mrSges(5,3);
t113 = -t100 * t12 + t96 * t11 - m(5) * (-t70 * pkin(3) + t105 - t135) - t76 * t124 + t69 * mrSges(5,3);
t77 = mrSges(5,1) * t85 + qJD(3) * mrSges(5,2);
t127 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t85 - t77;
t132 = mrSges(4,1) - mrSges(5,2);
t75 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t124;
t143 = (-t101 * t75 + t127 * t97) * qJD(2) - t132 * t70 + m(4) * (t110 - t135) + t69 * mrSges(4,2) - t113;
t129 = t101 * t39 + t97 * t53;
t142 = t103 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t144 - t66 * t124 - t129;
t8 = m(3) * t141 + qJDD(2) * mrSges(3,1) - t104 * mrSges(3,2) - t143;
t136 = t102 * t8;
t131 = mrSges(4,3) + mrSges(5,1);
t67 = (mrSges(5,2) * t101 - mrSges(5,3) * t97) * qJD(2);
t128 = t67 + (-mrSges(4,1) * t101 + mrSges(4,2) * t97) * qJD(2);
t126 = t101 * t53;
t108 = t70 * pkin(4) + qJD(3) * t79 - t89 * t137 - t142;
t111 = -t28 * mrSges(7,1) - t46 * t40 + m(7) * (-t44 * pkin(5) - t60 * pkin(10) + t65 * t52 + t108) + t29 * mrSges(7,2) + t47 * t41;
t107 = m(6) * t108 - t44 * mrSges(6,1) + t45 * mrSges(6,2) - t64 * t50 + t65 * t51 + t111;
t106 = -m(5) * t142 + t107;
t13 = t106 + t128 * t124 + t131 * t70 - t127 * qJD(3) + (-mrSges(4,2) + mrSges(5,3)) * qJDD(3) + m(4) * t129;
t112 = -m(5) * (-qJDD(3) * pkin(3) + t115 - t126) - t100 * t11 - t96 * t12;
t9 = m(4) * (-t36 + t126) - t131 * t69 + t132 * qJDD(3) + (t75 - t76) * qJD(3) - t128 * t85 + t112;
t4 = m(3) * t121 - t104 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t101 * t13 - t97 * t9;
t6 = m(3) * t53 + t101 * t9 + t97 * t13;
t120 = m(2) * t90 + t4 * t133 + t92 * t136 + t94 * t6;
t2 = m(2) * t73 + t102 * t4 - t98 * t8;
t1 = m(2) * t72 - t92 * t6 + (t4 * t98 + t136) * t94;
t3 = [-m(1) * g(1) - t91 * t1 + t93 * t2, t2, t4, t13, t70 * mrSges(5,2) - t77 * t85 - t113, t12, t15; -m(1) * g(2) + t93 * t1 + t91 * t2, t1, t8, t9, -t70 * mrSges(5,1) - qJDD(3) * mrSges(5,3) - qJD(3) * t77 - t67 * t124 - t106, t11, t14; -m(1) * g(3) + t120, t120, t6, t143, t69 * mrSges(5,1) + qJDD(3) * mrSges(5,2) + qJD(3) * t76 + t67 * t85 - t112, t107, t111;];
f_new  = t3;
