% Calculate time derivative of joint inertia matrix for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP7_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP7_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:10
% EndTime: 2019-12-31 17:20:13
% DurationCPUTime: 1.38s
% Computational Cost: add. (556->197), mult. (1462->294), div. (0->0), fcn. (930->4), ass. (0->92)
t93 = Ifges(5,4) + Ifges(4,5);
t61 = cos(qJ(3));
t86 = qJD(3) * t61;
t59 = sin(qJ(3));
t88 = qJD(3) * t59;
t113 = Ifges(5,6) * t88 + t93 * t86;
t112 = -Ifges(5,2) - Ifges(4,3);
t60 = sin(qJ(2));
t62 = cos(qJ(2));
t89 = qJD(2) * t62;
t82 = t59 * t89;
t64 = t60 * t86 + t82;
t101 = Ifges(5,5) * t59;
t73 = Ifges(5,1) * t61 + t101;
t103 = Ifges(4,4) * t59;
t74 = Ifges(4,1) * t61 - t103;
t111 = (t73 + t74) * qJD(3);
t40 = (pkin(2) * t60 - pkin(6) * t62) * qJD(2);
t42 = -pkin(2) * t62 - t60 * pkin(6) - pkin(1);
t85 = qJD(3) * t62;
t90 = qJD(2) * t60;
t5 = pkin(5) * (t59 * t90 - t61 * t85) + t61 * t40 - t42 * t88;
t70 = pkin(3) * t61 + qJ(4) * t59;
t84 = qJD(4) * t61;
t110 = t70 * qJD(3) - t84;
t109 = 2 * m(4);
t108 = -0.2e1 * pkin(1);
t107 = 0.2e1 * pkin(5);
t104 = t62 * pkin(5);
t102 = Ifges(4,4) * t61;
t100 = Ifges(5,5) * t61;
t99 = Ifges(4,6) * t61;
t98 = Ifges(4,6) * t62;
t97 = t59 * t60;
t96 = t60 * t61;
t95 = t61 * t42;
t94 = qJD(2) / 0.2e1;
t22 = -Ifges(5,4) * t62 + t60 * t73;
t23 = -Ifges(4,5) * t62 + t60 * t74;
t92 = t22 + t23;
t19 = t61 * t104 + t59 * t42;
t87 = qJD(3) * t60;
t81 = t61 * t89;
t44 = -Ifges(5,3) * t61 + t101;
t45 = Ifges(4,2) * t61 + t103;
t80 = t44 / 0.2e1 - t45 / 0.2e1;
t46 = Ifges(5,1) * t59 - t100;
t47 = Ifges(4,1) * t59 + t102;
t79 = t47 / 0.2e1 + t46 / 0.2e1;
t71 = Ifges(5,3) * t59 + t100;
t20 = -Ifges(5,6) * t62 + t60 * t71;
t72 = -Ifges(4,2) * t59 + t102;
t21 = t60 * t72 - t98;
t78 = t20 - t21 + t98;
t77 = -t61 * mrSges(4,1) + t59 * mrSges(4,2);
t76 = mrSges(4,1) * t59 + mrSges(4,2) * t61;
t43 = -t61 * mrSges(5,1) - t59 * mrSges(5,3);
t75 = mrSges(5,1) * t59 - mrSges(5,3) * t61;
t69 = pkin(3) * t59 - qJ(4) * t61;
t68 = -t64 * Ifges(5,6) + t112 * t90 - t81 * t93;
t67 = pkin(5) + t69;
t65 = -t59 * t87 + t81;
t15 = mrSges(5,2) * t81 + (-mrSges(5,1) * qJD(2) - mrSges(5,2) * t88) * t60;
t4 = t59 * t40 + t42 * t86 + (-t59 * t85 - t61 * t90) * pkin(5);
t41 = -pkin(2) - t70;
t39 = -mrSges(5,2) * t97 - mrSges(5,3) * t62;
t38 = mrSges(5,1) * t62 + mrSges(5,2) * t96;
t37 = -mrSges(4,1) * t62 - mrSges(4,3) * t96;
t36 = mrSges(4,2) * t62 - mrSges(4,3) * t97;
t33 = t72 * qJD(3);
t32 = t71 * qJD(3);
t31 = t76 * qJD(3);
t30 = t75 * qJD(3);
t27 = t75 * t60;
t25 = qJD(3) * t69 - qJD(4) * t59;
t24 = t67 * t60;
t18 = -t104 * t59 + t95;
t17 = -mrSges(5,2) * t64 + mrSges(5,3) * t90;
t16 = -mrSges(4,2) * t90 - mrSges(4,3) * t64;
t14 = mrSges(4,1) * t90 - mrSges(4,3) * t65;
t13 = -t95 + (pkin(5) * t59 + pkin(3)) * t62;
t12 = -qJ(4) * t62 + t19;
t11 = mrSges(4,1) * t64 + mrSges(4,2) * t65;
t10 = mrSges(5,1) * t64 - mrSges(5,3) * t65;
t9 = -t47 * t87 + (Ifges(4,5) * t60 + t62 * t74) * qJD(2);
t8 = -t46 * t87 + (Ifges(5,4) * t60 + t62 * t73) * qJD(2);
t7 = -t45 * t87 + (Ifges(4,6) * t60 + t62 * t72) * qJD(2);
t6 = -t44 * t87 + (Ifges(5,6) * t60 + t62 * t71) * qJD(2);
t3 = t110 * t60 + t67 * t89;
t2 = -pkin(3) * t90 - t5;
t1 = qJ(4) * t90 - qJD(4) * t62 + t4;
t26 = [0.2e1 * t1 * t39 + 0.2e1 * t24 * t10 + 0.2e1 * t12 * t17 + 0.2e1 * t13 * t15 + 0.2e1 * t18 * t14 + 0.2e1 * t19 * t16 + 0.2e1 * t2 * t38 + 0.2e1 * t3 * t27 + 0.2e1 * t4 * t36 + 0.2e1 * t5 * t37 + (t18 * t5 + t19 * t4) * t109 + 0.2e1 * m(5) * (t1 * t12 + t13 * t2 + t24 * t3) + ((mrSges(3,2) * t108 + 0.2e1 * Ifges(3,4) * t62 + t59 * t78 + t61 * t92) * qJD(2) + t68) * t62 + (t11 * t107 + (t8 + t9) * t61 + (t6 - t7) * t59 + (t78 * t61 + (t62 * t93 - t92) * t59) * qJD(3) + (mrSges(3,1) * t108 + (-0.2e1 * Ifges(3,4) + t93 * t61 + (-Ifges(4,6) + Ifges(5,6)) * t59) * t60 + (pkin(5) ^ 2 * t109 + t107 * t76 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) + t112) * t62) * qJD(2)) * t60; m(5) * (t24 * t25 + t3 * t41) + t24 * t30 + t41 * t10 + t3 * t43 + t25 * t27 - pkin(2) * t11 + (t8 / 0.2e1 + t9 / 0.2e1 - t5 * mrSges(4,3) + t2 * mrSges(5,2) + (-t12 * mrSges(5,2) - t19 * mrSges(4,3) + t98 / 0.2e1 + t20 / 0.2e1 - t21 / 0.2e1) * qJD(3)) * t59 + (-t6 / 0.2e1 + t7 / 0.2e1 + t4 * mrSges(4,3) + t1 * mrSges(5,2) + (t13 * mrSges(5,2) - t18 * mrSges(4,3) + t22 / 0.2e1 + t23 / 0.2e1) * qJD(3)) * t61 + (Ifges(3,5) + t79 * t61 + t80 * t59 + (-m(4) * pkin(2) - mrSges(3,1) + t77) * pkin(5)) * t89 + ((t16 + t17) * t61 + (-t14 + t15) * t59 + ((-t37 + t38) * t61 + (-t36 - t39) * t59) * qJD(3) + m(5) * (t1 * t61 - t12 * t88 + t13 * t86 + t2 * t59) + m(4) * (-t18 * t86 - t19 * t88 + t4 * t61 - t5 * t59)) * pkin(6) + (t99 * t94 - Ifges(3,6) * qJD(2) + (qJD(2) * mrSges(3,2) + t31) * pkin(5) + (t32 / 0.2e1 - t33 / 0.2e1 - t79 * qJD(3) + t93 * t94) * t59 + (-Ifges(5,6) * t94 + t80 * qJD(3) + t111 / 0.2e1) * t61) * t60 - t113 * t62 / 0.2e1; -0.2e1 * pkin(2) * t31 + 0.2e1 * t30 * t41 + (-t32 + t33) * t61 + t111 * t59 + 0.2e1 * (m(5) * t41 + t43) * t25 + ((t46 + t47) * t61 + (t44 - t45) * t59) * qJD(3); -Ifges(4,6) * t82 - pkin(3) * t15 + m(5) * (-pkin(3) * t2 + qJ(4) * t1 + qJD(4) * t12) + qJD(4) * t39 + qJ(4) * t17 + t1 * mrSges(5,3) - t2 * mrSges(5,1) - t4 * mrSges(4,2) + t5 * mrSges(4,1) + (-t93 * t59 - t99) * t87 - t68; -Ifges(4,6) * t88 - t110 * mrSges(5,2) + (m(5) * t84 + (-m(5) * t70 + t43 + t77) * qJD(3)) * pkin(6) + t113; 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4); m(5) * t2 + t15; (m(5) * pkin(6) + mrSges(5,2)) * t86; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t26(1), t26(2), t26(4), t26(7); t26(2), t26(3), t26(5), t26(8); t26(4), t26(5), t26(6), t26(9); t26(7), t26(8), t26(9), t26(10);];
Mq = res;
