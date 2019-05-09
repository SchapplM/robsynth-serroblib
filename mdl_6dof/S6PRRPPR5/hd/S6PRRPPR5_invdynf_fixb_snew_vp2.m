% Calculate vector of cutting forces with Newton-Euler
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-05-05 03:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:24:35
% EndTime: 2019-05-05 03:24:40
% DurationCPUTime: 2.08s
% Computational Cost: add. (23471->177), mult. (48957->230), div. (0->0), fcn. (30641->12), ass. (0->97)
t144 = -2 * qJD(4);
t100 = cos(qJ(3));
t101 = cos(qJ(2));
t91 = sin(pkin(10));
t94 = cos(pkin(10));
t73 = g(1) * t91 - g(2) * t94;
t95 = cos(pkin(6));
t136 = t73 * t95;
t74 = -g(1) * t94 - g(2) * t91;
t89 = -g(3) + qJDD(1);
t92 = sin(pkin(6));
t98 = sin(qJ(2));
t141 = (t89 * t92 + t136) * t101 - t98 * t74;
t109 = -qJDD(2) * pkin(2) - t141;
t123 = qJD(2) * qJD(3);
t117 = t100 * t123;
t97 = sin(qJ(3));
t118 = t97 * t123;
t125 = qJD(2) * t97;
t70 = qJDD(2) * t97 + t117;
t104 = pkin(3) * t118 + t125 * t144 + (-t70 - t117) * qJ(4) + t109;
t102 = qJD(3) ^ 2;
t103 = qJD(2) ^ 2;
t135 = t92 * t98;
t120 = t101 * t74 + t89 * t135 + t98 * t136;
t39 = -pkin(2) * t103 + qJDD(2) * pkin(8) + t120;
t36 = t97 * t39;
t67 = (-pkin(3) * t100 - qJ(4) * t97) * qJD(2);
t114 = -t102 * qJ(4) + t67 * t125 + qJDD(4) + t36;
t126 = qJ(5) * t103;
t132 = -pkin(3) - qJ(5);
t53 = -t73 * t92 + t89 * t95;
t24 = t70 * pkin(4) + t132 * qJDD(3) + (-pkin(4) * t123 - t97 * t126 - t53) * t100 + t114;
t71 = qJDD(2) * t100 - t118;
t77 = pkin(4) * t125 - qJD(3) * qJ(5);
t88 = t100 ^ 2;
t26 = -t77 * t125 + t132 * t71 + (-pkin(4) * t88 - pkin(8)) * t103 + t104;
t124 = qJD(2) * t100;
t90 = sin(pkin(11));
t93 = cos(pkin(11));
t62 = qJD(3) * t93 - t90 * t124;
t116 = -0.2e1 * qJD(5) * t62 + t93 * t24 - t90 * t26;
t51 = qJDD(3) * t93 - t71 * t90;
t61 = -qJD(3) * t90 - t93 * t124;
t16 = (t61 * t125 - t51) * pkin(9) + (t61 * t62 + t70) * pkin(5) + t116;
t121 = 0.2e1 * qJD(5) * t61 + t90 * t24 + t93 * t26;
t50 = -qJDD(3) * t90 - t71 * t93;
t52 = pkin(5) * t125 - t62 * pkin(9);
t60 = t61 ^ 2;
t17 = -t60 * pkin(5) + t50 * pkin(9) - t52 * t125 + t121;
t96 = sin(qJ(6));
t99 = cos(qJ(6));
t43 = t61 * t99 - t62 * t96;
t33 = t43 * qJD(6) + t50 * t96 + t51 * t99;
t44 = t61 * t96 + t62 * t99;
t35 = -mrSges(7,1) * t43 + mrSges(7,2) * t44;
t82 = qJD(6) + t125;
t40 = -mrSges(7,2) * t82 + t43 * mrSges(7,3);
t64 = qJDD(6) + t70;
t14 = m(7) * (t16 * t99 - t17 * t96) - t33 * mrSges(7,3) + t64 * mrSges(7,1) - t44 * t35 + t82 * t40;
t32 = -t44 * qJD(6) + t50 * t99 - t51 * t96;
t41 = mrSges(7,1) * t82 - t44 * mrSges(7,3);
t15 = m(7) * (t16 * t96 + t17 * t99) + t32 * mrSges(7,3) - t64 * mrSges(7,2) + t43 * t35 - t82 * t41;
t45 = -mrSges(6,1) * t61 + mrSges(6,2) * t62;
t48 = -mrSges(6,2) * t125 + t61 * mrSges(6,3);
t11 = m(6) * t116 + t70 * mrSges(6,1) - t51 * mrSges(6,3) + t48 * t125 + t99 * t14 + t96 * t15 - t62 * t45;
t49 = mrSges(6,1) * t125 - t62 * mrSges(6,3);
t12 = m(6) * t121 - t70 * mrSges(6,2) + t50 * mrSges(6,3) - t49 * t125 - t96 * t14 + t99 * t15 + t61 * t45;
t137 = t103 * pkin(8);
t78 = -mrSges(5,1) * t124 - qJD(3) * mrSges(5,3);
t112 = t90 * t11 - t93 * t12 - m(5) * (-t71 * pkin(3) + t104 - t137) - t78 * t124 + t70 * mrSges(5,3);
t79 = mrSges(5,1) * t125 + qJD(3) * mrSges(5,2);
t129 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t125 - t79;
t134 = mrSges(4,1) - mrSges(5,2);
t76 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t124;
t143 = qJD(2) * (-t100 * t76 + t129 * t97) - t134 * t71 + m(4) * (t109 - t137) + t70 * mrSges(4,2) - t112;
t131 = t100 * t39 + t97 * t53;
t142 = pkin(3) * t102 - qJDD(3) * qJ(4) + qJD(3) * t144 - t67 * t124 - t131;
t8 = m(3) * t141 + qJDD(2) * mrSges(3,1) - t103 * mrSges(3,2) - t143;
t138 = t101 * t8;
t133 = mrSges(4,3) + mrSges(5,1);
t68 = (mrSges(5,2) * t100 - mrSges(5,3) * t97) * qJD(2);
t130 = t68 + (-mrSges(4,1) * t100 + mrSges(4,2) * t97) * qJD(2);
t128 = t100 * t53;
t106 = pkin(4) * t71 + qJD(3) * t77 - t88 * t126 + qJDD(5) - t142;
t110 = -t32 * mrSges(7,1) - t43 * t40 + m(7) * (-t50 * pkin(5) - t60 * pkin(9) + t62 * t52 + t106) + t33 * mrSges(7,2) + t44 * t41;
t107 = m(6) * t106 - t50 * mrSges(6,1) + t51 * mrSges(6,2) - t61 * t48 + t62 * t49 + t110;
t105 = -m(5) * t142 + t107;
t13 = m(4) * t131 + t105 + t130 * t124 + t133 * t71 - t129 * qJD(3) + (-mrSges(4,2) + mrSges(5,3)) * qJDD(3);
t111 = -m(5) * (-qJDD(3) * pkin(3) + t114 - t128) - t93 * t11 - t90 * t12;
t9 = m(4) * (-t36 + t128) - t133 * t70 + t134 * qJDD(3) + (t76 - t78) * qJD(3) - t130 * t125 + t111;
t4 = m(3) * t120 - mrSges(3,1) * t103 - qJDD(2) * mrSges(3,2) + t100 * t13 - t9 * t97;
t6 = m(3) * t53 + t100 * t9 + t13 * t97;
t119 = m(2) * t89 + t135 * t4 + t138 * t92 + t6 * t95;
t2 = m(2) * t74 + t101 * t4 - t8 * t98;
t1 = m(2) * t73 - t92 * t6 + (t4 * t98 + t138) * t95;
t3 = [-m(1) * g(1) - t1 * t91 + t2 * t94, t2, t4, t13, t71 * mrSges(5,2) - t125 * t79 - t112, t12, t15; -m(1) * g(2) + t1 * t94 + t2 * t91, t1, t8, t9, -t71 * mrSges(5,1) - qJDD(3) * mrSges(5,3) - qJD(3) * t79 - t124 * t68 - t105, t11, t14; -m(1) * g(3) + t119, t119, t6, t143, t70 * mrSges(5,1) + qJDD(3) * mrSges(5,2) + qJD(3) * t78 + t125 * t68 - t111, t107, t110;];
f_new  = t3;
