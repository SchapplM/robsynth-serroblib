% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-05-05 22:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:38:26
% EndTime: 2019-05-05 22:38:34
% DurationCPUTime: 3.18s
% Computational Cost: add. (34061->193), mult. (82792->238), div. (0->0), fcn. (63336->10), ass. (0->100)
t102 = qJD(1) ^ 2;
t101 = cos(qJ(1));
t98 = sin(qJ(1));
t124 = t98 * g(1) - t101 * g(2);
t120 = qJDD(2) - t124;
t94 = cos(pkin(10));
t91 = t94 ^ 2;
t93 = sin(pkin(10));
t128 = t93 ^ 2 + t91;
t106 = (-pkin(2) * t94 - pkin(1)) * qJDD(1) + (-t128 * pkin(7) - qJ(2)) * t102 + t120;
t100 = cos(qJ(3));
t97 = sin(qJ(3));
t116 = t100 * t94 - t93 * t97;
t115 = t100 * t93 + t94 * t97;
t82 = t115 * qJD(1);
t72 = -t82 * qJD(3) + qJDD(1) * t116;
t76 = qJD(3) * pkin(3) - t82 * pkin(8);
t81 = t116 * qJD(1);
t80 = t81 ^ 2;
t105 = -t72 * pkin(3) - t80 * pkin(8) + t82 * t76 + t106;
t139 = cos(qJ(4));
t96 = sin(qJ(4));
t63 = -t139 * t81 + t96 * t82;
t92 = qJD(3) + qJD(4);
t137 = t63 * t92;
t142 = -2 * qJD(5);
t126 = t81 * qJD(3);
t73 = qJDD(1) * t115 + t126;
t38 = -t63 * qJD(4) + t139 * t73 + t96 * t72;
t64 = t139 * t82 + t96 * t81;
t103 = (-t38 + t137) * qJ(5) + t105 + (t92 * pkin(4) + t142) * t64;
t125 = qJD(1) * qJD(2);
t123 = -t94 * g(3) - 0.2e1 * t93 * t125;
t127 = pkin(7) * qJDD(1);
t138 = pkin(2) * t102;
t118 = -t101 * g(1) - t98 * g(2);
t83 = -t102 * pkin(1) + qJDD(1) * qJ(2) + t118;
t57 = (t94 * t138 - t127 - t83) * t93 + t123;
t121 = -t93 * g(3) + (0.2e1 * t125 + t83) * t94;
t59 = t94 * t127 - t91 * t138 + t121;
t122 = t100 * t57 - t97 * t59;
t23 = (-t73 + t126) * pkin(8) + (t81 * t82 + qJDD(3)) * pkin(3) + t122;
t129 = t100 * t59 + t97 * t57;
t29 = -t80 * pkin(3) + t72 * pkin(8) - qJD(3) * t76 + t129;
t119 = t139 * t23 - t96 * t29;
t44 = t63 * pkin(4) - t64 * qJ(5);
t88 = t92 ^ 2;
t89 = qJDD(3) + qJDD(4);
t19 = -t89 * pkin(4) - t88 * qJ(5) + t64 * t44 + qJDD(5) - t119;
t14 = (t63 * t64 - t89) * pkin(9) + (t38 + t137) * pkin(5) + t19;
t37 = t64 * qJD(4) - t139 * t72 + t96 * t73;
t55 = t64 * pkin(5) - t92 * pkin(9);
t62 = t63 ^ 2;
t15 = -t62 * pkin(5) - t64 * t55 + t103 + (pkin(4) + pkin(9)) * t37;
t95 = sin(qJ(6));
t99 = cos(qJ(6));
t47 = t99 * t63 - t95 * t92;
t26 = t47 * qJD(6) + t95 * t37 + t99 * t89;
t36 = qJDD(6) + t38;
t48 = t95 * t63 + t99 * t92;
t39 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t61 = qJD(6) + t64;
t40 = -t61 * mrSges(7,2) + t47 * mrSges(7,3);
t12 = m(7) * (t99 * t14 - t95 * t15) - t26 * mrSges(7,3) + t36 * mrSges(7,1) - t48 * t39 + t61 * t40;
t25 = -t48 * qJD(6) + t99 * t37 - t95 * t89;
t41 = t61 * mrSges(7,1) - t48 * mrSges(7,3);
t13 = m(7) * (t95 * t14 + t99 * t15) + t25 * mrSges(7,3) - t36 * mrSges(7,2) + t47 * t39 - t61 * t41;
t54 = t64 * mrSges(6,1) + t92 * mrSges(6,2);
t113 = t95 * t12 - t99 * t13 - m(6) * (t37 * pkin(4) + t103) + t38 * mrSges(6,3) + t64 * t54;
t53 = t63 * mrSges(6,1) - t92 * mrSges(6,3);
t130 = t92 * mrSges(5,2) + t63 * mrSges(5,3) + t53;
t134 = mrSges(5,1) - mrSges(6,2);
t52 = t92 * mrSges(5,1) - t64 * mrSges(5,3);
t143 = -m(5) * t105 - t38 * mrSges(5,2) + t130 * t63 - t134 * t37 - t64 * t52 + t113;
t74 = -qJD(3) * mrSges(4,2) + t81 * mrSges(4,3);
t75 = qJD(3) * mrSges(4,1) - t82 * mrSges(4,3);
t104 = m(4) * t106 - t72 * mrSges(4,1) + t73 * mrSges(4,2) - t81 * t74 + t82 * t75 - t143;
t146 = -t104 - m(3) * (-qJDD(1) * pkin(1) - t102 * qJ(2) + t120);
t145 = t128 * mrSges(3,3);
t117 = -t94 * mrSges(3,1) + t93 * mrSges(3,2);
t114 = qJDD(1) * mrSges(3,3) + t102 * t117;
t132 = t139 * t29 + t96 * t23;
t109 = -t88 * pkin(4) + t89 * qJ(5) - t63 * t44 + t132;
t111 = -t25 * mrSges(7,1) - t47 * t40 + m(7) * (-t37 * pkin(5) - t62 * pkin(9) + ((2 * qJD(5)) + t55) * t92 + t109) + t26 * mrSges(7,2) + t48 * t41;
t108 = -m(6) * (t92 * t142 - t109) + t111;
t46 = -t63 * mrSges(6,2) - t64 * mrSges(6,3);
t131 = -t63 * mrSges(5,1) - t64 * mrSges(5,2) - t46;
t133 = -mrSges(5,3) - mrSges(6,1);
t10 = m(5) * t132 + (-t52 + t54) * t92 + (-mrSges(5,2) + mrSges(6,3)) * t89 + t131 * t63 + t133 * t37 + t108;
t69 = -t81 * mrSges(4,1) + t82 * mrSges(4,2);
t112 = -m(6) * t19 - t99 * t12 - t95 * t13;
t9 = m(5) * t119 - t130 * t92 + t131 * t64 + t133 * t38 + t134 * t89 + t112;
t6 = m(4) * t122 + qJDD(3) * mrSges(4,1) - t73 * mrSges(4,3) + qJD(3) * t74 + t96 * t10 + t139 * t9 - t82 * t69;
t7 = m(4) * t129 - qJDD(3) * mrSges(4,2) + t72 * mrSges(4,3) - qJD(3) * t75 + t139 * t10 + t81 * t69 - t96 * t9;
t4 = m(3) * t123 + t97 * t7 + t100 * t6 + (-m(3) * t83 - t114) * t93;
t5 = m(3) * t121 + t100 * t7 + t114 * t94 - t97 * t6;
t141 = t94 * t4 + t93 * t5;
t8 = m(2) * t124 + (-mrSges(2,2) + t145) * t102 + (mrSges(2,1) - t117) * qJDD(1) + t146;
t1 = m(2) * t118 - t102 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t93 * t4 + t94 * t5;
t2 = [-m(1) * g(1) + t101 * t1 - t98 * t8, t1, t5, t7, t10, -t37 * mrSges(6,2) - t63 * t53 - t113, t13; -m(1) * g(2) + t98 * t1 + t101 * t8, t8, t4, t6, t9, t37 * mrSges(6,1) - t89 * mrSges(6,3) + t63 * t46 - t92 * t54 - t108, t12; (-m(1) - m(2)) * g(3) + t141, -m(2) * g(3) + t141, t117 * qJDD(1) - t102 * t145 - t146, t104, -t143, t38 * mrSges(6,1) + t89 * mrSges(6,2) + t64 * t46 + t92 * t53 - t112, t111;];
f_new  = t2;
