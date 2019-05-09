% Calculate vector of cutting forces with Newton-Euler
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-05-04 23:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:57:02
% EndTime: 2019-05-04 22:57:06
% DurationCPUTime: 1.89s
% Computational Cost: add. (21313->164), mult. (47433->208), div. (0->0), fcn. (34146->12), ass. (0->93)
t84 = sin(pkin(10));
t87 = cos(pkin(10));
t72 = t84 * g(1) - t87 * g(2);
t88 = cos(pkin(6));
t128 = t72 * t88;
t73 = -t87 * g(1) - t84 * g(2);
t82 = -g(3) + qJDD(1);
t85 = sin(pkin(6));
t91 = sin(qJ(2));
t93 = cos(qJ(2));
t134 = (t82 * t85 + t128) * t93 - t91 * t73;
t101 = qJDD(3) - t134;
t95 = qJD(2) ^ 2;
t130 = cos(qJ(4));
t83 = sin(pkin(11));
t86 = cos(pkin(11));
t90 = sin(qJ(4));
t135 = -t86 * t130 + t83 * t90;
t66 = t135 * qJD(2);
t117 = t66 * qJD(4);
t118 = pkin(8) * qJDD(2);
t115 = qJD(2) * qJD(3);
t61 = -t85 * t72 + t88 * t82;
t121 = -0.2e1 * t83 * t115 + t86 * t61;
t131 = pkin(3) * t95;
t126 = t85 * t91;
t113 = t82 * t126 + t91 * t128 + t93 * t73;
t41 = -t95 * pkin(2) + qJDD(2) * qJ(3) + t113;
t27 = (t86 * t131 - t118 - t41) * t83 + t121;
t114 = t83 * t61 + (0.2e1 * t115 + t41) * t86;
t81 = t86 ^ 2;
t28 = t86 * t118 - t81 * t131 + t114;
t110 = t130 * t27 - t90 * t28;
t105 = t130 * t83 + t86 * t90;
t67 = t105 * qJD(2);
t44 = t66 * pkin(4) - t67 * qJ(5);
t94 = qJD(4) ^ 2;
t21 = -qJDD(4) * pkin(4) - t94 * qJ(5) + t67 * t44 + qJDD(5) - t110;
t51 = t105 * qJDD(2) - t117;
t16 = (t66 * t67 - qJDD(4)) * pkin(9) + (t51 + t117) * pkin(5) + t21;
t116 = t67 * qJD(4);
t50 = t135 * qJDD(2) + t116;
t60 = t67 * pkin(5) - qJD(4) * pkin(9);
t65 = t66 ^ 2;
t133 = -2 * qJD(5);
t119 = t83 ^ 2 + t81;
t97 = (-pkin(3) * t86 - pkin(2)) * qJDD(2) + (-t119 * pkin(8) - qJ(3)) * t95 + t101;
t96 = pkin(4) * t116 + t67 * t133 + (-t51 + t117) * qJ(5) + t97;
t19 = -t65 * pkin(5) - t67 * t60 + (pkin(4) + pkin(9)) * t50 + t96;
t89 = sin(qJ(6));
t92 = cos(qJ(6));
t52 = -t89 * qJD(4) + t92 * t66;
t32 = t52 * qJD(6) + t92 * qJDD(4) + t89 * t50;
t53 = t92 * qJD(4) + t89 * t66;
t35 = -t52 * mrSges(7,1) + t53 * mrSges(7,2);
t64 = qJD(6) + t67;
t39 = -t64 * mrSges(7,2) + t52 * mrSges(7,3);
t49 = qJDD(6) + t51;
t14 = m(7) * (t92 * t16 - t89 * t19) - t32 * mrSges(7,3) + t49 * mrSges(7,1) - t53 * t35 + t64 * t39;
t31 = -t53 * qJD(6) - t89 * qJDD(4) + t92 * t50;
t40 = t64 * mrSges(7,1) - t53 * mrSges(7,3);
t15 = m(7) * (t89 * t16 + t92 * t19) + t31 * mrSges(7,3) - t49 * mrSges(7,2) + t52 * t35 - t64 * t40;
t59 = t67 * mrSges(6,1) + qJD(4) * mrSges(6,2);
t106 = t89 * t14 - t92 * t15 - m(6) * (t50 * pkin(4) + t96) + t67 * t59 + t51 * mrSges(6,3);
t58 = t66 * mrSges(6,1) - qJD(4) * mrSges(6,3);
t120 = -qJD(4) * mrSges(5,2) - t66 * mrSges(5,3) - t58;
t125 = mrSges(5,1) - mrSges(6,2);
t57 = qJD(4) * mrSges(5,1) - t67 * mrSges(5,3);
t98 = m(5) * t97 + t51 * mrSges(5,2) + t120 * t66 + t125 * t50 + t67 * t57 - t106;
t137 = -m(4) * (-qJDD(2) * pkin(2) - t95 * qJ(3) + t101) - t98;
t136 = t119 * mrSges(4,3);
t109 = -t86 * mrSges(4,1) + t83 * mrSges(4,2);
t11 = m(3) * t134 + (-mrSges(3,2) + t136) * t95 + (mrSges(3,1) - t109) * qJDD(2) + t137;
t129 = t11 * t93;
t124 = -mrSges(5,3) - mrSges(6,1);
t123 = t130 * t28 + t90 * t27;
t46 = -t66 * mrSges(6,2) - t67 * mrSges(6,3);
t122 = -t66 * mrSges(5,1) - t67 * mrSges(5,2) - t46;
t107 = qJDD(2) * mrSges(4,3) + t95 * t109;
t100 = -t94 * pkin(4) + qJDD(4) * qJ(5) - t66 * t44 + t123;
t103 = -t31 * mrSges(7,1) - t52 * t39 + m(7) * (-t50 * pkin(5) - t65 * pkin(9) + ((2 * qJD(5)) + t60) * qJD(4) + t100) + t32 * mrSges(7,2) + t53 * t40;
t99 = -m(6) * (qJD(4) * t133 - t100) + t103;
t12 = m(5) * t123 + t122 * t66 + t124 * t50 + (-mrSges(5,2) + mrSges(6,3)) * qJDD(4) + (-t57 + t59) * qJD(4) + t99;
t104 = -m(6) * t21 - t92 * t14 - t89 * t15;
t9 = m(5) * t110 + t120 * qJD(4) + t125 * qJDD(4) + t122 * t67 + t124 * t51 + t104;
t7 = m(4) * t121 + t90 * t12 + t130 * t9 + (-m(4) * t41 - t107) * t83;
t8 = m(4) * t114 + t107 * t86 + t130 * t12 - t90 * t9;
t4 = m(3) * t113 - t95 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t83 * t7 + t86 * t8;
t6 = m(3) * t61 + t86 * t7 + t83 * t8;
t111 = m(2) * t82 + t4 * t126 + t85 * t129 + t88 * t6;
t2 = m(2) * t73 - t91 * t11 + t93 * t4;
t1 = m(2) * t72 - t85 * t6 + (t4 * t91 + t129) * t88;
t3 = [-m(1) * g(1) - t84 * t1 + t87 * t2, t2, t4, t8, t12, -t50 * mrSges(6,2) - t66 * t58 - t106, t15; -m(1) * g(2) + t87 * t1 + t84 * t2, t1, t11, t7, t9, t50 * mrSges(6,1) - qJDD(4) * mrSges(6,3) - qJD(4) * t59 + t66 * t46 - t99, t14; -m(1) * g(3) + t111, t111, t6, t109 * qJDD(2) - t95 * t136 - t137, t98, t51 * mrSges(6,1) + qJDD(4) * mrSges(6,2) + qJD(4) * t58 + t67 * t46 - t104, t103;];
f_new  = t3;
