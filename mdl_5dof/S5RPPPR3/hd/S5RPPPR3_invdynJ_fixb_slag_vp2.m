% Calculate vector of inverse dynamics joint torques for
% S5RPPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:52
% EndTime: 2019-12-31 17:43:56
% DurationCPUTime: 3.38s
% Computational Cost: add. (1007->226), mult. (2073->304), div. (0->0), fcn. (1276->10), ass. (0->103)
t79 = sin(pkin(7));
t66 = pkin(1) * t79 + qJ(3);
t51 = qJD(1) * qJD(3) + qJDD(1) * t66;
t78 = sin(pkin(8));
t80 = cos(pkin(8));
t137 = t78 ^ 2 + t80 ^ 2;
t126 = m(5) + m(6);
t105 = m(4) + t126;
t136 = Ifges(4,4) - Ifges(5,5);
t118 = t78 * mrSges(5,3);
t93 = t80 * mrSges(5,1) + t118;
t94 = mrSges(4,1) * t80 - mrSges(4,2) * t78;
t135 = -t94 - t93 - mrSges(3,1);
t134 = m(6) * pkin(6) + mrSges(3,2) + mrSges(6,3);
t81 = cos(pkin(7));
t124 = pkin(1) * t81;
t115 = t78 * qJ(4);
t97 = pkin(2) + t115;
t129 = -(pkin(3) + pkin(4)) * t80 - t97;
t38 = t124 - t129;
t21 = qJD(1) * t38 - qJD(3);
t84 = cos(qJ(5));
t110 = qJD(1) * t84;
t82 = sin(qJ(5));
t111 = qJD(1) * t82;
t46 = -t110 * t80 - t111 * t78;
t47 = t110 * t78 - t111 * t80;
t133 = m(6) * t21 - mrSges(6,1) * t46 + mrSges(6,2) * t47 + t93 * qJD(1);
t77 = qJ(1) + pkin(7);
t72 = sin(t77);
t73 = cos(t77);
t132 = -g(1) * t73 - g(2) * t72;
t29 = t78 * qJDD(2) + t80 * t51;
t120 = t29 * t80;
t131 = t51 * t137 + t120 + t132;
t88 = -pkin(3) * t80 - t97;
t48 = t88 - t124;
t130 = qJDD(1) * t48;
t127 = t47 / 0.2e1;
t83 = sin(qJ(1));
t123 = pkin(1) * t83;
t85 = cos(qJ(1));
t74 = t85 * pkin(1);
t122 = -pkin(6) + t66;
t121 = Ifges(6,4) * t47;
t119 = t73 * t80;
t58 = t66 * qJD(1);
t40 = t78 * qJD(2) + t80 * t58;
t117 = t40 * t80 * qJD(3) + t66 * t120;
t113 = qJD(1) * t78;
t112 = qJD(1) * t80;
t109 = qJD(4) * t78;
t107 = qJDD(1) * t78;
t106 = qJDD(1) * t80;
t103 = t73 * pkin(2) + t72 * qJ(3) + t74;
t69 = -pkin(2) - t124;
t39 = qJD(2) * t80 - t78 * t58;
t28 = qJDD(2) * t80 - t78 * t51;
t96 = pkin(3) * t119 + t73 * t115 + t103;
t90 = t78 * t82 + t80 * t84;
t44 = t90 * qJD(5);
t55 = t78 * t84 - t80 * t82;
t15 = -qJD(1) * t44 + qJDD(1) * t55;
t45 = t55 * qJD(5);
t16 = -qJD(1) * t45 - qJDD(1) * t90;
t92 = -t16 * mrSges(6,1) + t15 * mrSges(6,2);
t91 = -mrSges(6,1) * t90 - mrSges(6,2) * t55;
t37 = qJD(4) - t39;
t23 = -pkin(6) * t113 + t37;
t24 = -pkin(6) * t112 + t40;
t3 = t23 * t84 - t24 * t82;
t4 = t23 * t82 + t24 * t84;
t49 = t122 * t78;
t50 = t122 * t80;
t11 = t49 * t84 - t50 * t82;
t12 = t49 * t82 + t50 * t84;
t25 = qJDD(4) - t28;
t89 = -qJD(1) * t109 + qJDD(3);
t68 = mrSges(4,2) * t107;
t57 = qJDD(1) * t69 + qJDD(3);
t41 = Ifges(6,4) * t46;
t36 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t47;
t35 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t46;
t34 = t90 * t73;
t33 = t55 * t73;
t32 = t90 * t72;
t31 = t55 * t72;
t30 = qJD(1) * t48 + qJD(3);
t27 = t40 * t112;
t22 = t29 * t78;
t20 = -pkin(6) * t106 + t29;
t19 = t89 + t130;
t18 = -pkin(6) * t107 + t25;
t13 = qJDD(1) * t38 - t89;
t10 = Ifges(6,1) * t47 + Ifges(6,5) * qJD(5) + t41;
t9 = Ifges(6,2) * t46 + Ifges(6,6) * qJD(5) + t121;
t8 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t16;
t7 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t15;
t6 = qJD(3) * t55 - qJD(5) * t12;
t5 = qJD(3) * t90 + qJD(5) * t11;
t2 = -qJD(5) * t4 + t18 * t84 - t20 * t82;
t1 = qJD(5) * t3 + t18 * t82 + t20 * t84;
t14 = [t21 * (mrSges(6,1) * t45 - mrSges(6,2) * t44) + qJD(5) * (-Ifges(6,5) * t44 - Ifges(6,6) * t45) / 0.2e1 + t46 * (-Ifges(6,4) * t44 - Ifges(6,2) * t45) / 0.2e1 + (-Ifges(6,1) * t44 - Ifges(6,4) * t45) * t127 + (t3 * t44 - t4 * t45) * mrSges(6,3) + (Ifges(2,3) + Ifges(3,3) + (0.2e1 * t81 * mrSges(3,1) - 0.2e1 * t79 * mrSges(3,2) + m(3) * (t79 ^ 2 + t81 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (-mrSges(6,3) * t1 - Ifges(6,4) * t15 - Ifges(6,2) * t16 - Ifges(6,6) * qJDD(5)) * t90 + t133 * t109 + (-m(4) * t103 - m(5) * t96 - m(3) * t74 - m(6) * (pkin(4) * t119 + t96) - t34 * mrSges(6,1) - t33 * mrSges(6,2) - mrSges(2,1) * t85 + mrSges(2,2) * t83 + t134 * t72 + t135 * t73) * g(2) + (-t28 * t78 + t131) * mrSges(4,3) + (t25 * t78 + t131) * mrSges(5,2) + m(4) * (t57 * t69 + (-qJD(3) * t39 - t28 * t66) * t78 + t117) + m(5) * (t19 * t48 + (qJD(3) * t37 - qJD(4) * t30 + t25 * t66) * t78 + t117) + m(6) * (t1 * t12 + t11 * t2 + t13 * t38 + t3 * t6 + t4 * t5) + (m(3) * t123 + mrSges(2,1) * t83 + t32 * mrSges(6,1) + mrSges(2,2) * t85 + t31 * mrSges(6,2) + t134 * t73 - t105 * (t73 * qJ(3) - t123) + (m(4) * pkin(2) - m(5) * t88 - t129 * m(6) - t135) * t72) * g(1) + (-mrSges(6,3) * t2 + Ifges(6,1) * t15 + Ifges(6,4) * t16 + Ifges(6,5) * qJDD(5)) * t55 + (-t19 - t130) * t93 + t69 * t68 - t44 * t10 / 0.2e1 - t45 * t9 / 0.2e1 + t5 * t35 + t6 * t36 + t11 * t7 + t12 * t8 - t13 * t91 + t38 * t92 - t57 * t94 + (t136 * t80 + (Ifges(4,1) + Ifges(5,1)) * t78) * t107 + (-t69 * mrSges(4,1) + (Ifges(4,2) + Ifges(5,3)) * t80 + t136 * t78) * t106; m(3) * qJDD(2) - t44 * t35 - t45 * t36 - t90 * t7 + t55 * t8 + m(4) * (t28 * t80 + t22) + m(5) * (-t25 * t80 + t22) + m(6) * (t1 * t55 - t2 * t90 - t3 * t45 - t4 * t44) + (-m(3) - t105) * g(3); t46 * t35 - t47 * t36 + t68 + (-t118 + (-mrSges(4,1) - mrSges(5,1)) * t80) * qJDD(1) - t92 - (mrSges(5,2) + mrSges(4,3)) * qJD(1) ^ 2 * t137 + (-g(1) * t72 + g(2) * t73) * t105 + (-t3 * t47 + t4 * t46 - t13) * m(6) + (-t113 * t37 + t19 - t27) * m(5) + (t113 * t39 - t27 + t57) * m(4); t84 * t7 + t82 * t8 + (t84 * t35 - t82 * t36) * qJD(5) + t126 * t80 * g(3) + m(5) * t25 + m(6) * (t1 * t82 + t2 * t84 + (-t3 * t82 + t4 * t84) * qJD(5)) + (qJDD(1) * mrSges(5,2) + (m(5) * t30 - t133) * qJD(1) + t126 * t132) * t78; Ifges(6,5) * t15 + Ifges(6,6) * t16 + Ifges(6,3) * qJDD(5) - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t21 * (mrSges(6,1) * t47 + mrSges(6,2) * t46) - t47 * (Ifges(6,1) * t46 - t121) / 0.2e1 + t9 * t127 - qJD(5) * (Ifges(6,5) * t46 - Ifges(6,6) * t47) / 0.2e1 - t3 * t35 + t4 * t36 - g(1) * (mrSges(6,1) * t33 - mrSges(6,2) * t34) - g(2) * (mrSges(6,1) * t31 - mrSges(6,2) * t32) - g(3) * t91 + (t3 * t46 + t4 * t47) * mrSges(6,3) - (-Ifges(6,2) * t47 + t10 + t41) * t46 / 0.2e1;];
tau = t14;
