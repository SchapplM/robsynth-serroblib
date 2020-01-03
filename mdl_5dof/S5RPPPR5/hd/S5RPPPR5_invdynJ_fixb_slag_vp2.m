% Calculate vector of inverse dynamics joint torques for
% S5RPPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:21
% EndTime: 2019-12-31 17:46:25
% DurationCPUTime: 2.46s
% Computational Cost: add. (1378->240), mult. (2568->319), div. (0->0), fcn. (1524->10), ass. (0->106)
t100 = sin(pkin(8));
t102 = cos(pkin(8));
t163 = (t100 ^ 2 + t102 ^ 2) * mrSges(5,3);
t162 = mrSges(4,2) - t163;
t101 = sin(pkin(7));
t115 = mrSges(5,1) * t102 - mrSges(5,2) * t100;
t132 = qJD(1) * t102;
t103 = cos(pkin(7));
t108 = -pkin(1) - pkin(2);
t69 = qJD(1) * t108 + qJD(2);
t89 = t101 * qJ(2);
t46 = -qJD(1) * t89 + t103 * t69;
t34 = qJD(1) * pkin(3) + qJD(4) - t46;
t29 = pkin(4) * t132 + t34;
t106 = sin(qJ(5));
t107 = cos(qJ(5));
t55 = t100 * t106 - t107 * t102;
t50 = t55 * qJD(1);
t57 = t100 * t107 + t102 * t106;
t51 = t57 * qJD(1);
t161 = (m(5) * t34 + m(6) * t29 - mrSges(6,1) * t50 - mrSges(6,2) * t51 + qJD(1) * t115) * t101;
t124 = qJDD(1) * t102;
t125 = qJDD(1) * t100;
t54 = mrSges(5,1) * t124 - mrSges(5,2) * t125;
t155 = t55 * qJD(5);
t24 = qJD(1) * t155 - qJDD(1) * t57;
t53 = t57 * qJD(5);
t25 = qJD(1) * t53 + qJDD(1) * t55;
t7 = -t25 * mrSges(6,1) + mrSges(6,2) * t24;
t160 = -t54 - t7;
t159 = mrSges(3,1) + mrSges(2,1);
t158 = -mrSges(3,3) + mrSges(2,2);
t131 = qJD(1) * t103;
t44 = t57 * t101;
t157 = -qJD(5) * t44 + t131 * t55;
t156 = t101 * t155 + t131 * t57;
t127 = qJD(1) * qJD(2);
t70 = qJDD(1) * qJ(2) + t127;
t68 = qJDD(1) * t108 + qJDD(2);
t33 = t101 * t68 + t103 * t70;
t28 = -qJDD(1) * qJ(4) - qJD(1) * qJD(4) + t33;
t84 = t102 * qJDD(3);
t12 = -t100 * t28 + t84;
t13 = qJDD(3) * t100 + t102 * t28;
t154 = -t12 * t100 + t13 * t102;
t47 = qJ(2) * t131 + t101 * t69;
t35 = -qJD(1) * qJ(4) + t47;
t86 = t102 * qJD(3);
t26 = -t100 * t35 + t86;
t27 = qJD(3) * t100 + t102 * t35;
t153 = -t100 * t26 + t102 * t27;
t152 = -m(5) - m(6) - m(4);
t98 = pkin(8) + qJ(5);
t87 = sin(t98);
t88 = cos(t98);
t117 = -mrSges(6,1) * t88 + mrSges(6,2) * t87;
t92 = t102 * pkin(4);
t151 = m(5) * pkin(3) + mrSges(4,1) + m(6) * (t92 + pkin(3)) - t117 + t115;
t147 = -m(5) * qJ(4) + m(6) * (-pkin(6) - qJ(4)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t145 = -t51 / 0.2e1;
t64 = qJ(2) * t103 + t101 * t108;
t59 = -qJ(4) + t64;
t144 = pkin(6) - t59;
t143 = cos(qJ(1));
t142 = sin(qJ(1));
t141 = Ifges(6,4) * t51;
t140 = pkin(1) * t143 + qJ(2) * t142;
t129 = qJDD(1) * mrSges(3,1);
t128 = qJDD(1) * mrSges(4,1);
t32 = -t101 * t70 + t103 * t68;
t63 = t103 * t108 - t89;
t61 = pkin(3) - t63;
t119 = -pkin(1) * t142 + qJ(2) * t143;
t56 = -t101 * t142 - t103 * t143;
t58 = t101 * t143 - t103 * t142;
t118 = -g(1) * t58 + g(2) * t56;
t114 = -t101 * t46 + t103 * t47;
t18 = t86 + (pkin(6) * qJD(1) - t35) * t100;
t19 = -pkin(6) * t132 + t27;
t3 = -t106 * t19 + t107 * t18;
t4 = t106 * t18 + t107 * t19;
t40 = t144 * t100;
t41 = t144 * t102;
t8 = t106 * t41 + t107 * t40;
t9 = t106 * t40 - t107 * t41;
t30 = qJDD(1) * pkin(3) + qJDD(4) - t32;
t109 = qJD(1) ^ 2;
t82 = -pkin(1) * qJDD(1) + qJDD(2);
t72 = qJD(2) * t103 - qJD(4);
t49 = t61 + t92;
t48 = Ifges(6,4) * t50;
t45 = t55 * t101;
t43 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t51;
t42 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t50;
t22 = pkin(4) * t124 + t30;
t17 = -Ifges(6,1) * t51 + Ifges(6,5) * qJD(5) + t48;
t16 = Ifges(6,2) * t50 + Ifges(6,6) * qJD(5) - t141;
t15 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t25;
t14 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t24;
t11 = -pkin(6) * t124 + t13;
t10 = t84 + (pkin(6) * qJDD(1) - t28) * t100;
t6 = -qJD(5) * t9 - t57 * t72;
t5 = qJD(5) * t8 - t55 * t72;
t2 = -qJD(5) * t4 + t10 * t107 - t106 * t11;
t1 = qJD(5) * t3 + t10 * t106 + t107 * t11;
t20 = [m(5) * (t30 * t61 + (t13 * t59 + t27 * t72) * t102 + (-t12 * t59 - t26 * t72) * t100) + (t101 * t127 - t32) * mrSges(4,1) + (-m(3) * t119 + t158 * t143 + t159 * t142 - t151 * t58 + t152 * (-pkin(2) * t142 + t119) + t147 * t56) * g(1) + (-m(3) * t140 - t159 * t143 + t158 * t142 + t152 * (pkin(2) * t143 + t140) + t147 * t58 + t151 * t56) * g(2) - t154 * mrSges(5,3) + (t64 * mrSges(4,2) - t163 * t59 + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1) + (-mrSges(6,2) * t22 + mrSges(6,3) * t2 - Ifges(6,1) * t24 - Ifges(6,4) * t25 - Ifges(6,5) * qJDD(5)) * t57 + m(6) * (t1 * t9 + t2 * t8 + t22 * t49 + t3 * t6 + t4 * t5) - (-Ifges(5,1) * t100 - Ifges(5,4) * t102) * t125 + m(3) * (-pkin(1) * t82 + (t70 + t127) * qJ(2)) - t63 * t128 - (-Ifges(5,4) * t100 - Ifges(5,2) * t102) * t124 - qJD(1) * t72 * t163 + m(4) * (qJD(2) * t114 + t32 * t63 + t33 * t64) + t30 * t115 + 0.2e1 * t70 * mrSges(3,3) - t82 * mrSges(3,1) + t61 * t54 + t53 * t16 / 0.2e1 + t5 * t42 + t6 * t43 + t49 * t7 + t8 * t14 + t9 * t15 + (-t155 * t3 + t4 * t53) * mrSges(6,3) + t155 * t17 / 0.2e1 + t50 * (Ifges(6,4) * t155 + Ifges(6,2) * t53) / 0.2e1 + t29 * (-mrSges(6,1) * t53 + mrSges(6,2) * t155) + qJD(5) * (Ifges(6,5) * t155 + Ifges(6,6) * t53) / 0.2e1 + (Ifges(6,1) * t155 + Ifges(6,4) * t53) * t145 + (-mrSges(6,1) * t22 + mrSges(6,3) * t1 + Ifges(6,4) * t24 + Ifges(6,2) * t25 + Ifges(6,6) * qJDD(5)) * t55 + qJD(2) * t161 + pkin(1) * t129 + (t103 * t127 + t33) * mrSges(4,2); -t129 - t44 * t14 - t45 * t15 + t156 * t43 + t157 * t42 + (-m(3) * qJ(2) - mrSges(3,3)) * t109 + (-t109 * t162 - t128 + t160) * t103 + (-t109 * mrSges(4,1) + qJDD(1) * t162) * t101 + m(4) * (t101 * t33 + t103 * t32) + m(5) * (t101 * t154 - t30 * t103) + m(3) * t82 + (-t1 * t45 - t103 * t22 + t156 * t3 + t157 * t4 - t2 * t44) * m(6) + (-m(5) * t103 * t153 - m(4) * t114 - t161) * qJD(1) + (-g(1) * t142 + g(2) * t143) * (m(3) - t152); m(4) * qJDD(3) - t55 * t14 + t57 * t15 - t155 * t42 - t53 * t43 + m(5) * (t100 * t13 + t102 * t12) + m(6) * (t1 * t57 - t155 * t4 - t2 * t55 - t3 * t53) - t152 * g(3); -t109 * t163 - t50 * t42 - t51 * t43 + (-t3 * t51 - t4 * t50 + t118 + t22) * m(6) + (qJD(1) * t153 + t118 + t30) * m(5) - t160; Ifges(6,5) * t24 + Ifges(6,6) * t25 + Ifges(6,3) * qJDD(5) - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t29 * (-mrSges(6,1) * t51 + mrSges(6,2) * t50) + t51 * (Ifges(6,1) * t50 + t141) / 0.2e1 + t16 * t145 - qJD(5) * (Ifges(6,5) * t50 + Ifges(6,6) * t51) / 0.2e1 - t3 * t42 + t4 * t43 - g(3) * t117 + (t3 * t50 - t4 * t51) * mrSges(6,3) - (Ifges(6,2) * t51 + t17 + t48) * t50 / 0.2e1 + (-g(1) * t56 - g(2) * t58) * (mrSges(6,1) * t87 + mrSges(6,2) * t88);];
tau = t20;
