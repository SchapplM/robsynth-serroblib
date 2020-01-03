% Calculate vector of inverse dynamics joint torques for
% S5PRPRR9
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:37
% EndTime: 2019-12-31 17:39:40
% DurationCPUTime: 1.44s
% Computational Cost: add. (1272->187), mult. (1817->244), div. (0->0), fcn. (826->6), ass. (0->90)
t62 = sin(qJ(5));
t64 = cos(qJ(5));
t40 = -mrSges(6,1) * t64 + mrSges(6,2) * t62;
t145 = mrSges(5,1) - t40;
t108 = t62 * mrSges(6,3);
t61 = -qJD(2) + qJD(4);
t35 = qJD(5) * mrSges(6,1) - t61 * t108;
t105 = t64 * mrSges(6,3);
t36 = -qJD(5) * mrSges(6,2) + t61 * t105;
t66 = -pkin(2) - pkin(3);
t48 = t66 * qJD(2) + qJD(3);
t63 = sin(qJ(4));
t65 = cos(qJ(4));
t98 = qJ(3) * qJD(2);
t21 = t48 * t63 + t65 * t98;
t13 = pkin(7) * t61 + t21;
t10 = -qJD(1) * t64 - t13 * t62;
t75 = t62 * qJD(1) - t13 * t64;
t77 = t10 * t62 + t64 * t75;
t124 = m(6) * t77 + t61 * mrSges(5,2) + t62 * t35 - t64 * t36;
t144 = -m(5) * t21 + t124;
t143 = m(6) * pkin(4) + t145;
t20 = t48 * t65 - t63 * t98;
t12 = -pkin(4) * t61 - t20;
t85 = t145 * t61;
t142 = -m(6) * t12 + t85;
t141 = qJD(4) * t20;
t122 = t62 / 0.2e1;
t140 = m(4) + m(5);
t139 = m(5) * t20 + t142;
t80 = mrSges(6,1) * t62 + mrSges(6,2) * t64;
t135 = t12 * t80;
t133 = mrSges(4,1) + mrSges(3,1);
t132 = -mrSges(4,3) + mrSges(3,2);
t38 = t65 * qJ(3) + t63 * t66;
t100 = qJD(5) * t62;
t99 = qJD(5) * t64;
t131 = -t10 * t99 + t100 * t75;
t96 = qJD(2) * qJD(3);
t50 = qJDD(2) * qJ(3) + t96;
t101 = qJD(4) * t21;
t44 = t66 * qJDD(2) + qJDD(3);
t6 = t44 * t65 - t50 * t63 - t101;
t60 = -qJDD(2) + qJDD(4);
t25 = -t61 * t100 + t60 * t64;
t14 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t25;
t26 = t60 * t62 + t61 * t99;
t15 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t26;
t130 = t64 * t14 - t62 * t15;
t4 = -pkin(4) * t60 - t6;
t129 = m(6) * t4 - mrSges(6,1) * t25 + mrSges(6,2) * t26;
t94 = pkin(8) + qJ(2);
t57 = cos(t94);
t84 = sin(t94);
t22 = -t57 * t65 - t84 * t63;
t23 = t57 * t63 - t84 * t65;
t128 = t23 * mrSges(5,2) + t143 * t22;
t127 = t22 * mrSges(5,2) - t143 * t23;
t78 = t10 * t64 - t62 * t75;
t5 = t63 * t44 + t65 * t50 + t141;
t3 = pkin(7) * t60 + t5;
t1 = t10 * qJD(5) - qJDD(1) * t62 + t3 * t64;
t2 = t75 * qJD(5) - qJDD(1) * t64 - t3 * t62;
t83 = t1 * t64 - t2 * t62;
t123 = m(6) * (-t78 * qJD(5) + t83) + t130 - t35 * t99 - t36 * t100;
t120 = Ifges(6,4) * t62;
t119 = Ifges(6,4) * t64;
t118 = Ifges(6,2) * t64;
t114 = t60 * mrSges(5,1);
t113 = t60 * mrSges(5,2);
t110 = t61 * t62;
t109 = t61 * t64;
t102 = t57 * pkin(2) + t84 * qJ(3);
t97 = qJDD(2) * mrSges(4,1);
t93 = t57 * pkin(3) + t102;
t88 = m(2) + m(3) + t140;
t87 = -m(6) * pkin(7) - mrSges(6,3);
t79 = Ifges(6,1) * t64 - t120;
t42 = t118 + t120;
t37 = -qJ(3) * t63 + t65 * t66;
t74 = -t84 * pkin(2) + t57 * qJ(3);
t29 = qJD(5) * (Ifges(6,5) * t64 - Ifges(6,6) * t62);
t69 = -t84 * pkin(3) + t74;
t16 = Ifges(6,6) * qJD(5) + t42 * t61;
t49 = Ifges(6,4) * t109;
t17 = Ifges(6,1) * t110 + Ifges(6,5) * qJD(5) + t49;
t68 = Ifges(5,3) * t60 + qJD(5) * t135 + t4 * t40 + t6 * mrSges(5,1) + t1 * t105 + t25 * t42 / 0.2e1 + t26 * (Ifges(6,1) * t62 + t119) / 0.2e1 - t5 * mrSges(5,2) + (Ifges(6,1) * t26 + Ifges(6,4) * t25) * t122 + t64 * (Ifges(6,4) * t26 + Ifges(6,2) * t25) / 0.2e1 - t16 * t100 / 0.2e1 + t17 * t99 / 0.2e1 - t2 * t108 + t131 * mrSges(6,3) + (0.2e1 * Ifges(6,5) * t122 + Ifges(6,6) * t64) * qJDD(5) + (t29 + t79 * t110 + (-Ifges(6,2) * t62 + t119) * t109) * qJD(5) / 0.2e1;
t67 = qJD(2) ^ 2;
t56 = -qJDD(2) * pkin(2) + qJDD(3);
t7 = [-t36 * t99 + t35 * t100 - t62 * t14 - t64 * t15 + m(6) * (t77 * qJD(5) - t1 * t62 - t2 * t64) + t88 * qJDD(1) + (-m(6) - t88) * g(3); pkin(2) * t97 + t37 * t114 + m(5) * (t37 * t6 + t38 * t5) - t38 * t113 + m(4) * (-pkin(2) * t56 + (t50 + t96) * qJ(3)) - t56 * mrSges(4,1) - t68 + t129 * (pkin(4) - t37) - t139 * (qJD(3) * t63 + t38 * qJD(4)) + (Ifges(4,2) + Ifges(3,3)) * qJDD(2) - t144 * (qJD(3) * t65 + t37 * qJD(4)) + 0.2e1 * t50 * mrSges(4,3) + t123 * (-pkin(7) + t38) + (-m(4) * t102 - m(6) * (pkin(7) * t23 + t93) - t23 * mrSges(6,3) - m(5) * t93 + t132 * t84 - t133 * t57 + t128) * g(2) + (-m(4) * t74 - m(5) * t69 - m(6) * (t22 * pkin(7) + t69) - t22 * mrSges(6,3) + t133 * t84 + t132 * t57 + t127) * g(1); -t97 - t67 * mrSges(4,3) + (-t67 * qJ(3) + t56) * m(4) + (t114 + m(5) * (t6 + t101) - t124 * qJD(4) - t129) * t65 + (-t113 + (-t64 * t35 - t62 * t36) * qJD(5) - t85 * qJD(4) + m(5) * (t5 - t141) + m(6) * (qJD(4) * t12 + t131 + t83) + t130) * t63 + (t139 * t63 + t144 * t65) * qJD(2) + (m(6) + t140) * (-g(1) * t84 + g(2) * t57); t68 + (-t87 * t23 - t128) * g(2) + (-t87 * t22 - t127) * g(1) - t129 * pkin(4) + t142 * t21 + t124 * t20 + t123 * pkin(7); t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t26 + Ifges(6,6) * t25 + Ifges(6,3) * qJDD(5) - g(3) * t40 - t10 * t36 - t75 * t35 + (-t135 + t16 * t122 - t29 / 0.2e1 + (-t62 * t79 / 0.2e1 + t118 * t122) * t61 + t78 * mrSges(6,3) - (t17 + t49) * t64 / 0.2e1) * t61 + (-t22 * g(1) - t23 * g(2)) * t80;];
tau = t7;
