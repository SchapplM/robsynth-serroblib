% Calculate vector of inverse dynamics joint torques for
% S5RPPRP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:44
% EndTime: 2019-12-31 17:50:49
% DurationCPUTime: 2.60s
% Computational Cost: add. (904->222), mult. (1588->272), div. (0->0), fcn. (696->8), ass. (0->99)
t146 = Ifges(6,4) + Ifges(5,4);
t154 = Ifges(6,1) + Ifges(5,1);
t153 = -Ifges(6,2) - Ifges(5,2);
t62 = cos(qJ(4));
t156 = t146 * t62;
t60 = sin(qJ(4));
t155 = t146 * t60;
t145 = Ifges(6,5) + Ifges(5,5);
t144 = Ifges(6,6) + Ifges(5,6);
t152 = t153 * t60 + t156;
t151 = t154 * t62 - t155;
t150 = -t60 / 0.2e1;
t128 = t62 / 0.2e1;
t127 = -m(4) - m(5);
t99 = m(6) - t127;
t56 = qJ(1) + pkin(7);
t52 = sin(t56);
t53 = cos(t56);
t149 = -g(1) * t52 + g(2) * t53;
t57 = sin(pkin(7));
t48 = pkin(1) * t57 + qJ(3);
t37 = pkin(4) * t60 + t48;
t21 = qJD(1) * t37 + qJD(5);
t38 = t48 * qJD(1);
t148 = t38 * (mrSges(5,1) * t62 - mrSges(5,2) * t60) + t21 * (mrSges(6,1) * t62 - mrSges(6,2) * t60) + (-t144 * t62 - t145 * t60) * qJD(4) / 0.2e1;
t102 = qJD(4) * t62;
t58 = cos(pkin(7));
t50 = -pkin(1) * t58 - pkin(2);
t45 = -pkin(6) + t50;
t28 = qJDD(1) * t45 + qJDD(3);
t29 = qJD(1) * t45 + qJD(3);
t93 = t62 * qJDD(2) + t29 * t102 + t60 * t28;
t95 = qJD(2) * qJD(4);
t3 = -t60 * t95 + t93;
t12 = qJD(2) * t62 + t29 * t60;
t4 = -qJD(4) * t12 - qJDD(2) * t60 + t62 * t28;
t139 = t3 * t60 + t4 * t62;
t11 = -qJD(2) * t60 + t62 * t29;
t74 = -t11 * t60 + t12 * t62;
t147 = m(5) * (qJD(4) * t74 + t139);
t143 = t152 * qJD(1) + t144 * qJD(4);
t142 = t151 * qJD(1) + t145 * qJD(4);
t97 = qJD(1) * qJD(4);
t34 = qJDD(1) * t62 - t60 * t97;
t96 = qJD(1) * qJD(5);
t1 = qJDD(4) * pkin(4) - qJ(5) * t34 - t62 * t96 + t4;
t35 = -qJDD(1) * t60 - t62 * t97;
t2 = qJ(5) * t35 + (-t95 - t96) * t60 + t93;
t138 = t1 * t62 + t2 * t60;
t136 = mrSges(3,1) + mrSges(5,3) + mrSges(6,3) - mrSges(4,2);
t115 = t62 * mrSges(6,2);
t82 = mrSges(5,1) * t60 + mrSges(5,2) * t62;
t91 = m(6) * pkin(4) + mrSges(6,1);
t135 = -t60 * t91 + mrSges(3,2) - mrSges(4,3) - t115 - t82;
t133 = (t153 * t62 - t155) * t150 + (-t154 * t60 - t156) * t128;
t61 = sin(qJ(1));
t126 = pkin(1) * t61;
t63 = cos(qJ(1));
t55 = t63 * pkin(1);
t16 = qJDD(4) * mrSges(6,1) - mrSges(6,3) * t34;
t17 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t34;
t112 = t16 + t17;
t18 = -qJDD(4) * mrSges(6,2) + mrSges(6,3) * t35;
t19 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t35;
t111 = t18 + t19;
t105 = qJD(1) * t60;
t39 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t105;
t40 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t105;
t110 = t39 + t40;
t104 = qJD(1) * t62;
t41 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t104;
t42 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t104;
t109 = -t41 - t42;
t107 = qJ(5) - t45;
t103 = qJD(4) * t60;
t101 = qJ(5) * qJD(1);
t100 = qJDD(1) * mrSges(4,2);
t98 = qJD(1) * qJD(3);
t81 = mrSges(6,1) * t60 + t115;
t32 = t81 * qJD(1);
t90 = -m(6) * t21 - t32;
t86 = -t35 * mrSges(6,1) + t34 * mrSges(6,2);
t23 = t107 * t62;
t7 = -t101 * t62 + t11;
t5 = qJD(4) * pkin(4) + t7;
t8 = -t101 * t60 + t12;
t84 = -t5 * t62 - t60 * t8;
t83 = -t5 * t60 + t8 * t62;
t30 = qJDD(1) * t48 + t98;
t73 = qJD(3) * t38 + t30 * t48;
t59 = -qJ(5) - pkin(6);
t44 = pkin(4) * t102 + qJD(3);
t36 = qJDD(1) * t50 + qJDD(3);
t33 = t82 * qJD(1);
t22 = t107 * t60;
t10 = -qJD(4) * t23 - qJD(5) * t60;
t9 = -qJD(5) * t62 + t103 * t107;
t6 = -pkin(4) * t35 + qJDD(5) + t30;
t13 = [(m(3) * t126 + mrSges(2,1) * t61 + mrSges(2,2) * t63 - t99 * (t53 * qJ(3) - t126) + t135 * t53 + (-m(6) * (-pkin(2) + t59) - m(5) * (-pkin(2) - pkin(6)) + m(4) * pkin(2) + t136) * t52) * g(1) + (-m(3) * t55 - mrSges(2,1) * t63 + mrSges(2,2) * t61 - t99 * (t53 * pkin(2) + t52 * qJ(3) + t55) + (-m(5) * pkin(6) + m(6) * t59 - t136) * t53 + t135 * t52) * g(2) + t50 * t100 + m(5) * t73 + m(6) * (-t1 * t23 + t10 * t8 - t2 * t22 + t21 * t44 + t37 * t6 + t5 * t9) + (t98 + t30) * mrSges(4,3) + (t40 * t102 - t42 * t103 + t62 * t17 + t60 * t19 + t147) * t45 - t142 * t103 / 0.2e1 - t143 * t102 / 0.2e1 + (-t102 * t8 + t103 * t5 - t138) * mrSges(6,3) + (-t102 * t12 + t103 * t11 - t139) * mrSges(5,3) + t30 * t82 + (t48 * mrSges(4,3) + Ifges(4,1) + Ifges(2,3) + Ifges(3,3) + (0.2e1 * mrSges(3,1) * t58 - 0.2e1 * t57 * mrSges(3,2) + m(3) * (t57 ^ 2 + t58 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t133 * t97 + (0.2e1 * t145 * t128 - t144 * t60) * qJDD(4) + t148 * qJD(4) + t37 * t86 + t44 * t32 + t48 * (-mrSges(5,1) * t35 + mrSges(5,2) * t34) + t36 * mrSges(4,2) + t10 * t39 + t9 * t41 - t22 * t18 - t23 * t16 + qJD(3) * t33 + t6 * t81 + m(4) * (t36 * t50 + t73) + t151 * t34 / 0.2e1 + t152 * t35 / 0.2e1 + (t146 * t34 - t153 * t35) * t150 + (t146 * t35 + t154 * t34) * t128; t111 * t62 - t112 * t60 + (m(3) + m(4)) * qJDD(2) + (t109 * t62 - t110 * t60) * qJD(4) + m(5) * (t3 * t62 - t4 * t60 + (-t11 * t62 - t12 * t60) * qJD(4)) + m(6) * (qJD(4) * t84 - t1 * t60 + t2 * t62) + (-m(3) - t99) * g(3); t100 + t112 * t62 + t111 * t60 + (t109 * t60 + t110 * t62) * qJD(4) + m(4) * t36 + t147 + m(6) * (qJD(4) * t83 + t138) + (-mrSges(4,3) * qJD(1) + t127 * t38 - t33 + t90) * qJD(1) + t149 * t99; t4 * mrSges(5,1) + t1 * mrSges(6,1) - t3 * mrSges(5,2) - t2 * mrSges(6,2) - t11 * t40 + t12 * t42 - t7 * t39 + t144 * t35 + t145 * t34 + (Ifges(6,3) + Ifges(5,3)) * qJDD(4) + (t82 + t81) * g(3) + (t16 + (g(3) * t60 + t1) * m(6)) * pkin(4) + (t41 - m(6) * (-t5 + t7)) * t8 + (-t133 * qJD(1) + t90 * t62 * pkin(4) + t83 * mrSges(6,3) + t74 * mrSges(5,3) + t142 * t60 / 0.2e1 + t143 * t128 - t148) * qJD(1) - t149 * ((mrSges(5,2) + mrSges(6,2)) * t60 - t62 * (mrSges(5,1) + t91)); (t60 * t39 + t62 * t41) * qJD(1) + (-g(1) * t53 - g(2) * t52 - t84 * qJD(1) + t6) * m(6) + t86;];
tau = t13;
