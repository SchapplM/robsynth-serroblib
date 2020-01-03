% Calculate time derivative of joint inertia matrix for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR9_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR9_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:19
% EndTime: 2019-12-31 18:02:22
% DurationCPUTime: 0.96s
% Computational Cost: add. (710->173), mult. (1501->282), div. (0->0), fcn. (1094->6), ass. (0->87)
t45 = sin(qJ(5));
t46 = sin(qJ(4));
t47 = cos(qJ(5));
t78 = qJD(5) * t47;
t48 = cos(qJ(4));
t81 = qJD(4) * t48;
t53 = t45 * t81 + t46 * t78;
t108 = Ifges(5,1) - Ifges(5,2);
t84 = t45 ^ 2 + t47 ^ 2;
t79 = qJD(5) * t46;
t68 = t45 * t79;
t54 = t47 * t81 - t68;
t57 = mrSges(6,1) * t45 + mrSges(6,2) * t47;
t19 = t57 * t46;
t7 = -mrSges(6,1) * t53 - mrSges(6,2) * t54;
t107 = -t19 * t81 + t46 * t7;
t29 = -mrSges(6,1) * t47 + mrSges(6,2) * t45;
t60 = -m(6) * pkin(4) - mrSges(5,1) + t29;
t106 = 0.2e1 * m(6);
t95 = Ifges(6,4) * t45;
t30 = Ifges(6,2) * t47 + t95;
t103 = t30 / 0.2e1;
t102 = t45 / 0.2e1;
t101 = -t47 / 0.2e1;
t100 = pkin(4) * t46;
t99 = pkin(7) * t48;
t94 = Ifges(6,4) * t47;
t93 = Ifges(6,6) * t45;
t92 = Ifges(6,6) * t48;
t43 = sin(pkin(8));
t44 = cos(pkin(8));
t49 = -pkin(1) - pkin(2);
t85 = t44 * qJ(2) + t43 * t49;
t21 = -pkin(6) + t85;
t40 = t46 ^ 2;
t91 = t21 * t40;
t42 = t48 ^ 2;
t90 = t21 * t42;
t89 = t45 * t46;
t88 = t45 * t48;
t87 = t46 * t47;
t86 = t47 * t48;
t83 = qJD(2) * t44;
t82 = qJD(4) * t46;
t80 = qJD(5) * t45;
t77 = t43 * qJD(2);
t74 = t44 * t77;
t73 = t48 * t83;
t72 = t43 * t82;
t70 = t46 * t81;
t63 = -t43 * qJ(2) + t44 * t49;
t20 = pkin(3) - t63;
t12 = t48 * pkin(4) + t46 * pkin(7) + t20;
t52 = qJD(5) * t12 - t21 * t82 + t73;
t59 = -qJD(5) * t21 * t48 + t77 + (t99 - t100) * qJD(4);
t1 = t45 * t59 + t47 * t52;
t3 = t12 * t47 - t21 * t88;
t64 = -qJD(5) * t3 + t1;
t62 = 0.2e1 * t70;
t58 = -t46 * mrSges(5,1) - t48 * mrSges(5,2);
t56 = Ifges(6,1) * t47 - t95;
t31 = Ifges(6,1) * t45 + t94;
t55 = -Ifges(6,2) * t45 + t94;
t18 = t43 * t86 - t44 * t45;
t17 = -t43 * t88 - t44 * t47;
t51 = Ifges(6,5) * t68 + (-Ifges(6,5) * t86 - Ifges(6,3) * t46) * qJD(4) + t53 * Ifges(6,6);
t8 = qJD(5) * t17 - t47 * t72;
t9 = -qJD(5) * t18 + t45 * t72;
t50 = -t9 * t45 + t8 * t47 + (-t17 * t47 - t18 * t45) * qJD(5);
t36 = Ifges(6,5) * t78;
t28 = t40 * t74;
t27 = mrSges(6,1) * t48 + mrSges(6,3) * t87;
t26 = -mrSges(6,2) * t48 + mrSges(6,3) * t89;
t25 = t56 * qJD(5);
t24 = t55 * qJD(5);
t23 = t58 * qJD(4);
t22 = t57 * qJD(5);
t15 = Ifges(6,5) * t48 - t46 * t56;
t14 = -t46 * t55 + t92;
t13 = t83 * t91;
t11 = mrSges(6,2) * t82 + mrSges(6,3) * t53;
t10 = -mrSges(6,1) * t82 + mrSges(6,3) * t54;
t6 = t31 * t79 + (-Ifges(6,5) * t46 - t48 * t56) * qJD(4);
t5 = t30 * t79 + (-Ifges(6,6) * t46 - t48 * t55) * qJD(4);
t4 = t12 * t45 + t21 * t86;
t2 = -t45 * t52 + t47 * t59;
t16 = [-t6 * t87 + 0.2e1 * m(5) * (t13 + (t20 * t43 + t44 * t90) * qJD(2)) + t48 * t51 + 0.2e1 * t20 * t23 + 0.2e1 * t1 * t26 + 0.2e1 * t2 * t27 + 0.2e1 * t3 * t10 + 0.2e1 * t4 * t11 + (t21 ^ 2 * t70 + t1 * t4 + t2 * t3 + t13) * t106 + t5 * t89 + (0.2e1 * Ifges(5,4) * t48 + t108 * t46) * t81 + 0.2e1 * (-t46 * t19 + mrSges(4,2)) * t83 + 0.2e1 * (t48 * mrSges(5,1) - mrSges(5,2) * t46 + mrSges(4,1)) * t77 - 0.2e1 * (t40 + t42) * mrSges(5,3) * t83 - t54 * t15 + t53 * t14 + 0.2e1 * t107 * t21 + ((-Ifges(6,3) + t108) * t48 + (Ifges(6,5) * t47 - 0.2e1 * Ifges(5,4) - t93) * t46) * t82 + 0.2e1 * (m(4) * (-t43 * t63 + t44 * t85) + m(3) * qJ(2) + mrSges(3,3)) * qJD(2); t17 * t10 + t18 * t11 - t44 * t23 + t8 * t26 + t9 * t27 + t107 * t43 + m(6) * (t21 * t43 * t62 + t1 * t18 + t17 * t2 + t3 * t9 + t4 * t8 + t28) + m(5) * (t28 + (t42 - 0.1e1) * t74); (t43 ^ 2 * t70 + t17 * t9 + t18 * t8) * t106; -t48 * t7 + (m(6) * (-t3 * t88 + t4 * t86 - t90 + t91) + t26 * t86 - t27 * t88) * qJD(4) + (-qJD(4) * t19 + m(6) * (t1 * t47 - t2 * t45 - t3 * t78 - t4 * t80 - t73) - t26 * t80 + t47 * t11 - t27 * t78 - t45 * t10) * t46; m(6) * (t50 * t46 + ((-t17 * t45 + t18 * t47) * t48 + (t40 - t42) * t43) * qJD(4)); m(6) * (-0.1e1 + t84) * t62; -pkin(4) * t7 + (-mrSges(5,2) * t83 + t36 / 0.2e1 + (t21 * t60 - Ifges(5,5)) * qJD(4)) * t48 + (t81 * t103 - t2 * mrSges(6,3) + t6 / 0.2e1 + (-t14 / 0.2e1 - t92 / 0.2e1 - t4 * mrSges(6,3)) * qJD(5) + (m(6) * (-qJD(5) * t4 - t2) - qJD(5) * t26 - t10) * pkin(7)) * t45 + (-t31 * t81 / 0.2e1 + qJD(5) * t15 / 0.2e1 + t5 / 0.2e1 + t64 * mrSges(6,3) + (m(6) * t64 - qJD(5) * t27 + t11) * pkin(7)) * t47 + (t25 * t101 + t24 * t102 + t21 * t22 + (t102 * t31 + t103 * t47) * qJD(5) + (-Ifges(6,5) * t45 / 0.2e1 + Ifges(6,6) * t101 + Ifges(5,6) + t21 * mrSges(5,2)) * qJD(4) + t60 * t83) * t46; (m(6) * pkin(7) + mrSges(6,3)) * t50 + ((qJD(4) * mrSges(5,2) + t22) * t46 + t60 * t81) * t43; -t48 * t22 + (t46 * t29 + m(6) * (t84 * t99 - t100) + t84 * t48 * mrSges(6,3) + t58) * qJD(4); -0.2e1 * pkin(4) * t22 + t47 * t24 + t45 * t25 + (-t30 * t45 + t31 * t47) * qJD(5); mrSges(6,1) * t2 - mrSges(6,2) * t1 + t51; mrSges(6,1) * t9 - mrSges(6,2) * t8; t7; t36 + (pkin(7) * t29 - t93) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t16(1), t16(2), t16(4), t16(7), t16(11); t16(2), t16(3), t16(5), t16(8), t16(12); t16(4), t16(5), t16(6), t16(9), t16(13); t16(7), t16(8), t16(9), t16(10), t16(14); t16(11), t16(12), t16(13), t16(14), t16(15);];
Mq = res;
