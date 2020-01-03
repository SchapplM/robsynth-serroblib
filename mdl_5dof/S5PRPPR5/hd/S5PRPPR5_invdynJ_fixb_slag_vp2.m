% Calculate vector of inverse dynamics joint torques for
% S5PRPPR5
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:03
% EndTime: 2019-12-31 17:38:08
% DurationCPUTime: 2.92s
% Computational Cost: add. (1029->230), mult. (1916->310), div. (0->0), fcn. (1107->8), ass. (0->103)
t73 = -pkin(2) - pkin(3);
t72 = cos(qJ(2));
t117 = qJD(1) * t72;
t98 = qJD(3) - t117;
t34 = qJD(2) * t73 + t98;
t70 = sin(qJ(2));
t118 = qJD(1) * t70;
t52 = qJD(2) * qJ(3) + t118;
t65 = sin(pkin(8));
t67 = cos(pkin(8));
t11 = t34 * t67 - t52 * t65;
t12 = t65 * t34 + t67 * t52;
t69 = sin(qJ(5));
t71 = cos(qJ(5));
t91 = mrSges(6,1) * t71 - mrSges(6,2) * t69;
t38 = t91 * qJD(2);
t160 = -m(4) * t52 - m(5) * (-t11 * t65 + t12 * t67) - t65 * t38;
t119 = mrSges(6,3) * qJD(2);
t46 = qJD(5) * mrSges(6,1) + t119 * t69;
t47 = -qJD(5) * mrSges(6,2) - t119 * t71;
t147 = -t69 * t46 + t71 * t47;
t159 = qJD(2) * mrSges(5,2) + t147;
t108 = m(4) + m(5) + m(6);
t78 = m(6) * pkin(4) + t91;
t83 = t65 * t72 - t67 * t70;
t158 = t83 * (-t78 - mrSges(5,1)) + (m(4) * pkin(2) - m(5) * t73 + mrSges(4,1)) * t70 + (-qJ(3) * t108 - mrSges(4,3)) * t72;
t112 = qJD(5) * t71;
t113 = qJD(5) * t69;
t10 = -qJD(2) * pkin(6) + t12;
t7 = qJD(4) * t71 - t10 * t69;
t8 = qJD(4) * t69 + t10 * t71;
t107 = qJD(1) * qJD(2);
t57 = t70 * t107;
t41 = qJDD(1) * t72 - t57;
t82 = qJDD(3) - t41;
t16 = qJDD(2) * t73 + t82;
t106 = qJD(2) * qJD(3);
t149 = qJDD(2) * qJ(3) + t106;
t58 = t72 * t107;
t42 = t70 * qJDD(1) + t58;
t21 = t42 + t149;
t6 = t65 * t16 + t67 * t21;
t4 = -qJDD(2) * pkin(6) + t6;
t1 = qJD(5) * t7 + qJDD(4) * t69 + t4 * t71;
t2 = -qJD(5) * t8 + qJDD(4) * t71 - t4 * t69;
t96 = t1 * t71 - t2 * t69;
t157 = -t7 * t112 - t8 * t113 + t96;
t123 = mrSges(3,1) + mrSges(4,1);
t129 = t67 * t72;
t32 = t65 * t70 + t129;
t156 = -t32 * mrSges(5,1) + mrSges(5,2) * t83 - t123 * t72;
t143 = -t69 / 0.2e1;
t94 = -t69 * t7 + t71 * t8;
t153 = m(6) * t94;
t152 = g(3) * t83;
t151 = qJD(2) * mrSges(5,1) + t38;
t105 = qJD(2) * qJD(5);
t39 = -qJDD(2) * t71 + t105 * t69;
t22 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t39;
t40 = -qJDD(2) * t69 - t105 * t71;
t23 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t40;
t148 = t71 * t22 - t69 * t23;
t66 = sin(pkin(7));
t68 = cos(pkin(7));
t146 = g(1) * t68 + g(2) * t66;
t9 = qJD(2) * pkin(4) - t11;
t90 = mrSges(6,1) * t69 + mrSges(6,2) * t71;
t144 = -t9 * t90 + qJD(5) * (-Ifges(6,5) * t71 + Ifges(6,6) * t69) / 0.2e1;
t138 = t65 * t9;
t135 = Ifges(6,4) * t69;
t134 = Ifges(6,4) * t71;
t131 = t66 * t70;
t128 = t68 * t70;
t122 = -mrSges(3,2) + mrSges(4,3);
t44 = t67 * qJ(3) + t65 * t73;
t121 = t72 * pkin(2) + t70 * qJ(3);
t111 = qJDD(2) * mrSges(4,1);
t110 = qJDD(2) * mrSges(5,1);
t109 = qJDD(2) * mrSges(5,2);
t103 = t72 * pkin(3) + t121;
t97 = g(1) * t66 - g(2) * t68;
t95 = -t8 * t69 - t7 * t71;
t89 = -t69 * Ifges(6,1) - t134;
t88 = -Ifges(6,2) * t71 - t135;
t5 = t16 * t67 - t21 * t65;
t85 = (-qJD(2) * pkin(2) + t98) * t70 + t52 * t72;
t43 = -qJ(3) * t65 + t67 * t73;
t80 = t69 * (-Ifges(6,1) * t71 + t135);
t79 = t71 * (Ifges(6,2) * t69 - t134);
t76 = qJD(5) * t95 + t96;
t75 = (-t71 * t46 - t69 * t47) * qJD(5) + t148;
t74 = qJD(2) ^ 2;
t35 = pkin(4) - t43;
t31 = Ifges(6,5) * qJD(5) + qJD(2) * t89;
t30 = Ifges(6,6) * qJD(5) + qJD(2) * t88;
t29 = t32 * qJD(2);
t28 = t83 * qJD(2);
t24 = -qJDD(2) * pkin(2) + t82;
t19 = -t128 * t65 - t129 * t68;
t17 = -t129 * t66 - t131 * t65;
t13 = -mrSges(6,1) * t39 + mrSges(6,2) * t40;
t3 = qJDD(2) * pkin(4) - t5;
t14 = [m(2) * qJDD(1) + t32 * t13 + t151 * t28 + t159 * t29 + (t122 * t72 - t123 * t70) * t74 - t75 * t83 + (-m(2) - m(3) - t108) * g(3) + m(3) * (t41 * t72 + t42 * t70) + m(4) * (qJD(2) * t85 + t21 * t70 - t24 * t72) + m(6) * (t28 * t9 + t29 * t94 + t3 * t32 - t76 * t83) + m(5) * (-t11 * t28 + t12 * t29 - t32 * t5 - t6 * t83) + (t122 * t70 - t156) * qJDD(2); -(t80 + t79) * t105 / 0.2e1 + (-g(1) * t19 - g(2) * t17 - t152 - t157) * mrSges(6,3) + (-m(6) * (pkin(6) * t83 + t103) - t78 * t32 - m(4) * t121 - m(5) * t103 + t70 * mrSges(3,2) + t156) * g(3) + (m(6) * t35 + t91) * t3 + (t41 + t57) * mrSges(3,1) + (-pkin(2) * t24 + qJ(3) * t21 - t85 * qJD(1)) * m(4) + (-m(5) * t12 - t153 - t159) * t32 * qJD(1) + t144 * qJD(5) + (-m(6) * (pkin(6) * t17 + t131 * t73) + t17 * mrSges(5,2) + t158 * t66) * g(2) + (-m(6) * (pkin(6) * t19 + t128 * t73) + t19 * mrSges(5,2) + t158 * t68) * g(1) - t31 * t112 / 0.2e1 + t30 * t113 / 0.2e1 - t43 * t110 + t146 * (t70 * mrSges(3,1) + mrSges(3,2) * t72) + (m(6) * t76 - t112 * t46 - t113 * t47 + t148) * (-pkin(6) + t44) + (-g(3) * t70 + t149 + t21 - t58) * mrSges(4,3) + (m(5) * t11 - m(6) * t9 - t151) * (t117 * t65 - t118 * t67) + (Ifges(6,1) * t40 + Ifges(6,4) * t39) * t143 + t39 * t88 / 0.2e1 + t40 * t89 / 0.2e1 + (-t42 + t58) * mrSges(3,2) + (m(6) * (t67 * t94 + t138) + t147 * t67 - t160) * qJD(3) + t44 * t109 + pkin(2) * t111 + (t67 * t106 + t6) * mrSges(5,2) + (t65 * t106 - t5) * mrSges(5,1) - t71 * (Ifges(6,4) * t40 + Ifges(6,2) * t39) / 0.2e1 + t35 * t13 + (-t24 + t57) * mrSges(4,1) + (Ifges(5,3) + Ifges(4,2) + Ifges(3,3)) * qJDD(2) + m(5) * (t43 * t5 + t44 * t6) + (0.2e1 * Ifges(6,5) * t143 - Ifges(6,6) * t71) * qJDD(5); -t74 * mrSges(4,3) + m(4) * t24 - t111 + (m(5) * t5 - m(6) * t3 - t74 * mrSges(5,2) - t110 - t13) * t67 + (m(5) * t6 + m(6) * t157 - t74 * mrSges(5,1) + t109 + t75) * t65 + (-m(6) * t138 + (-t147 - t153) * t67 + t160) * qJD(2) + (g(3) * t72 - t146 * t70) * t108; t69 * t22 + t71 * t23 + t147 * qJD(5) + (qJD(5) * t94 + t1 * t69 + t2 * t71 + t97) * m(6) + (qJDD(4) + t97) * m(5); Ifges(6,5) * t40 + Ifges(6,6) * t39 + Ifges(6,3) * qJDD(5) - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t7 * t47 + t8 * t46 - g(1) * ((t19 * t69 - t66 * t71) * mrSges(6,1) + (t19 * t71 + t66 * t69) * mrSges(6,2)) - g(2) * ((t17 * t69 + t68 * t71) * mrSges(6,1) + (t17 * t71 - t68 * t69) * mrSges(6,2)) - t90 * t152 + (t71 * t31 / 0.2e1 + t30 * t143 + (t80 / 0.2e1 + t79 / 0.2e1) * qJD(2) + t95 * mrSges(6,3) - t144) * qJD(2);];
tau = t14;
