% Calculate vector of inverse dynamics joint torques for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:11
% EndTime: 2019-07-18 13:30:14
% DurationCPUTime: 1.25s
% Computational Cost: add. (1610->201), mult. (2564->280), div. (0->0), fcn. (1272->12), ass. (0->102)
t78 = sin(qJ(5));
t147 = t78 / 0.2e1;
t155 = mrSges(5,2) - mrSges(6,3);
t82 = cos(qJ(5));
t52 = -mrSges(6,1) * t82 + mrSges(6,2) * t78;
t154 = t52 - mrSges(5,1);
t123 = pkin(2) * qJD(2);
t80 = sin(qJ(3));
t116 = t80 * t123;
t79 = sin(qJ(4));
t83 = cos(qJ(4));
t117 = qJD(2) + qJD(3);
t84 = cos(qJ(3));
t96 = pkin(3) * t117 + t123 * t84;
t153 = -t79 * t116 + t83 * t96;
t25 = t116 * t83 + t79 * t96;
t71 = qJD(4) + t117;
t19 = t71 * pkin(6) + t25;
t12 = qJD(1) * t82 - t19 * t78;
t118 = qJD(2) * qJD(3);
t46 = (qJDD(2) * t80 + t118 * t84) * pkin(2);
t45 = (qJDD(2) * t84 - t118 * t80) * pkin(2);
t76 = qJDD(2) + qJDD(3);
t92 = pkin(3) * t76 + t45;
t7 = t153 * qJD(4) + t83 * t46 + t79 * t92;
t70 = qJDD(4) + t76;
t5 = pkin(6) * t70 + t7;
t2 = qJD(5) * t12 + qJDD(1) * t78 + t5 * t82;
t13 = qJD(1) * t78 + t19 * t82;
t3 = -qJD(5) * t13 + qJDD(1) * t82 - t5 * t78;
t152 = t2 * t82 - t3 * t78;
t77 = qJ(2) + qJ(3);
t74 = qJ(4) + t77;
t65 = sin(t74);
t66 = cos(t74);
t151 = g(1) * t66 + g(2) * t65;
t109 = mrSges(6,1) * t78 + mrSges(6,2) * t82;
t149 = -t153 * t109 + qJD(5) * (Ifges(6,5) * t82 - Ifges(6,6) * t78) / 0.2e1 - (t12 * t82 + t13 * t78) * mrSges(6,3);
t148 = t45 * mrSges(4,1) + Ifges(4,3) * t76;
t105 = -t12 * t78 + t13 * t82;
t130 = t71 * t78;
t47 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t130;
t129 = t71 * t82;
t48 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t129;
t97 = t71 * mrSges(5,2) + t78 * t47 - t82 * t48;
t145 = m(5) + m(6);
t72 = sin(t77);
t143 = g(1) * t72;
t81 = sin(qJ(2));
t142 = g(1) * t81;
t136 = Ifges(6,4) * t78;
t135 = Ifges(6,4) * t82;
t134 = Ifges(6,2) * t82;
t126 = t80 * t83;
t102 = t79 * t84 + t126;
t36 = t102 * t123;
t131 = t153 * t36;
t128 = t78 * t83;
t127 = t79 * t80;
t125 = t82 * t83;
t68 = t84 * pkin(2) + pkin(3);
t40 = pkin(2) * t126 + t79 * t68;
t73 = cos(t77);
t64 = pkin(3) * t73;
t85 = cos(qJ(2));
t124 = t85 * pkin(2) + t64;
t122 = qJD(4) * t79;
t121 = qJD(4) * t83;
t120 = qJD(5) * t78;
t119 = qJD(5) * t82;
t115 = m(2) + m(3) + m(4) + m(5);
t32 = t52 * t71;
t113 = t71 * mrSges(5,1) - t32;
t111 = qJD(3) * t117;
t16 = t68 * t122 + (qJD(3) * t102 + t121 * t80) * pkin(2);
t39 = pkin(2) * t127 - t83 * t68;
t8 = qJD(4) * t25 + t79 * t46 - t83 * t92;
t110 = -t153 * t16 + t39 * t8;
t108 = t134 + t136;
t33 = -t120 * t71 + t70 * t82;
t22 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t33;
t34 = t119 * t71 + t70 * t78;
t23 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t34;
t104 = t82 * t22 - t78 * t23;
t101 = t83 * t84 - t127;
t100 = t117 * t123;
t98 = t78 * (Ifges(6,1) * t82 - t136);
t95 = t154 * t66 + t155 * t65;
t93 = -t154 * t65 + t155 * t66;
t91 = -t119 * t12 - t120 * t13 + t152;
t90 = -t73 * mrSges(4,1) + t72 * mrSges(4,2) + t95;
t89 = -m(6) * t66 * pkin(6) + t73 * mrSges(4,2) + t93;
t26 = Ifges(6,6) * qJD(5) + t108 * t71;
t53 = Ifges(6,4) * t129;
t27 = Ifges(6,1) * t130 + Ifges(6,5) * qJD(5) + t53;
t88 = t82 * (Ifges(6,4) * t34 + Ifges(6,2) * t33) / 0.2e1 + t33 * t108 / 0.2e1 + t34 * (Ifges(6,1) * t78 + t135) / 0.2e1 - t26 * t120 / 0.2e1 + Ifges(5,3) * t70 + (Ifges(6,1) * t34 + Ifges(6,4) * t33) * t147 + t154 * t8 + (t27 + t71 * (-Ifges(6,2) * t78 + t135)) * t119 / 0.2e1 + t152 * mrSges(6,3) + (0.2e1 * Ifges(6,5) * t147 + Ifges(6,6) * t82) * qJDD(5) + (t98 * t71 / 0.2e1 + t149) * qJD(5);
t87 = m(6) * t91 - t119 * t47 - t120 * t48 + t104;
t86 = -t7 * mrSges(5,2) + t88;
t59 = t65 * pkin(6);
t37 = t101 * t123;
t11 = -mrSges(6,1) * t33 + mrSges(6,2) * t34;
t1 = [t48 * t119 - t47 * t120 + t78 * t22 + t82 * t23 + m(6) * (qJD(5) * t105 + t2 * t78 + t3 * t82) + t115 * qJDD(1) + (-m(6) - t115) * g(3); m(5) * (t40 * t7 + t110) + (-t16 * t71 - t39 * t70) * mrSges(5,1) + (t145 * t142 + (-t111 * t84 - t80 * t76) * mrSges(4,2) + (-t111 * t80 + t84 * t76) * mrSges(4,1) + (-g(2) * t85 + t45 * t84 + t46 * t80 + t142) * m(4)) * pkin(2) + t87 * (pkin(6) + t40) + t88 + (t81 * mrSges(3,1) + t85 * mrSges(3,2) + (pkin(3) * t145 + mrSges(4,1)) * t72 + t89) * g(1) + (-t85 * mrSges(3,1) + t81 * mrSges(3,2) - m(5) * t124 - m(6) * (t59 + t124) + t90) * g(2) + t39 * t11 - t46 * mrSges(4,2) + t16 * t32 + (-t40 * t70 - t7) * mrSges(5,2) + Ifges(3,3) * qJDD(2) + m(6) * t110 + (m(5) * t25 + m(6) * t105 - t97) * (t68 * t121 + (qJD(3) * t101 - t80 * t122) * pkin(2)) + t148; (-m(6) * (t59 + t64) + t90) * g(2) - m(5) * (t25 * t37 - t131) + t87 * (pkin(3) * t79 + pkin(6)) + (m(6) * t143 - t79 * t70 * mrSges(5,2) + (-m(6) * t8 + t70 * mrSges(5,1) - t11) * t83 + (m(6) * (-t12 * t128 + t125 * t13 - t153 * t79) + t79 * t32 - t47 * t128 + t48 * t125 + (-t79 * mrSges(5,1) - t83 * mrSges(5,2)) * t71) * qJD(4) + (-g(2) * t73 + t121 * t25 - t122 * t153 + t7 * t79 - t8 * t83 + t143) * m(5)) * pkin(3) + t86 - m(6) * (t105 * t37 - t131) + t97 * t37 + t113 * t36 + (t100 * t84 - t46) * mrSges(4,2) + t80 * mrSges(4,1) * t100 + (t72 * mrSges(4,1) + t89) * g(1) + t148; -(-m(6) * (-t105 + t25) - t97) * t153 + ((-t82 * t47 - t78 * t48) * qJD(5) + (t91 - t151) * m(6) + t104) * pkin(6) + t93 * g(1) + t113 * t25 + t95 * g(2) + t86; t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t34 + Ifges(6,6) * t33 + Ifges(6,3) * qJDD(5) + g(3) * t52 - t12 * t48 + t13 * t47 + (t26 * t147 + (-t98 / 0.2e1 + t134 * t147) * t71 - (t27 + t53) * t82 / 0.2e1 - t149) * t71 + t151 * t109;];
tau  = t1;
