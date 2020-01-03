% Calculate time derivative of joint inertia matrix for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR7_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR7_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:58
% EndTime: 2019-12-31 19:35:01
% DurationCPUTime: 1.02s
% Computational Cost: add. (1206->187), mult. (2669->275), div. (0->0), fcn. (2319->6), ass. (0->95)
t56 = sin(pkin(8));
t57 = cos(pkin(8));
t59 = sin(qJ(2));
t61 = cos(qJ(2));
t38 = t56 * t61 + t57 * t59;
t35 = t38 * qJD(2);
t90 = t57 * t61;
t37 = t56 * t59 - t90;
t58 = sin(qJ(5));
t60 = cos(qJ(5));
t80 = qJD(5) * t60;
t65 = t58 * t35 + t37 * t80;
t110 = 2 * mrSges(5,1) + 2 * mrSges(4,3);
t101 = pkin(3) + pkin(7);
t53 = -pkin(2) * t61 - pkin(1);
t66 = -t38 * qJ(4) + t53;
t16 = t101 * t37 + t66;
t83 = -qJ(3) - pkin(6);
t43 = t83 * t59;
t45 = t83 * t61;
t25 = -t57 * t43 - t45 * t56;
t19 = pkin(4) * t38 + t25;
t3 = -t16 * t58 + t19 * t60;
t82 = qJD(2) * t59;
t36 = qJD(2) * t90 - t56 * t82;
t55 = pkin(2) * t82;
t63 = -t36 * qJ(4) - t38 * qJD(4) + t55;
t8 = t101 * t35 + t63;
t76 = qJD(2) * t83;
t34 = qJD(3) * t61 + t59 * t76;
t62 = -t59 * qJD(3) + t61 * t76;
t17 = t34 * t56 - t57 * t62;
t9 = pkin(4) * t36 + t17;
t1 = qJD(5) * t3 + t58 * t9 + t60 * t8;
t4 = t16 * t60 + t19 * t58;
t2 = -qJD(5) * t4 - t58 * t8 + t60 * t9;
t109 = t1 * t58 + t2 * t60;
t13 = pkin(3) * t35 + t63;
t108 = -0.2e1 * t13;
t22 = pkin(3) * t37 + t66;
t107 = -0.2e1 * t22;
t106 = 0.2e1 * t53;
t105 = m(4) * pkin(2);
t104 = t37 / 0.2e1;
t103 = -t58 / 0.2e1;
t102 = t60 / 0.2e1;
t100 = pkin(2) * t56;
t99 = pkin(2) * t57;
t96 = mrSges(6,3) * t37;
t95 = Ifges(6,4) * t58;
t94 = Ifges(6,4) * t60;
t93 = Ifges(6,6) * t38;
t70 = Ifges(6,2) * t60 + t95;
t14 = t37 * t70 + t93;
t92 = t14 * t60;
t71 = Ifges(6,1) * t58 + t94;
t15 = Ifges(6,5) * t38 + t37 * t71;
t91 = t15 * t58;
t47 = Ifges(6,1) * t60 - t95;
t88 = t58 * t47;
t87 = t60 * t35;
t46 = -Ifges(6,2) * t58 + t94;
t86 = t60 * t46;
t84 = mrSges(5,2) - mrSges(4,1);
t81 = qJD(5) * t58;
t79 = 0.2e1 * t61;
t78 = t37 * t81;
t52 = -pkin(3) - t99;
t75 = Ifges(6,5) * t65 + Ifges(6,6) * t87 + Ifges(6,3) * t36;
t72 = t3 * t58 - t4 * t60;
t44 = mrSges(6,1) * t58 + mrSges(6,2) * t60;
t69 = -Ifges(6,5) * t58 - Ifges(6,6) * t60;
t18 = t57 * t34 + t56 * t62;
t26 = t43 * t56 - t45 * t57;
t68 = t17 * t25 + t18 * t26;
t23 = mrSges(6,1) * t38 - t58 * t96;
t24 = -mrSges(6,2) * t38 + t60 * t96;
t67 = -t58 * t23 + t60 * t24;
t64 = t78 - t87;
t50 = qJ(4) + t100;
t49 = -pkin(7) + t52;
t42 = t71 * qJD(5);
t41 = t70 * qJD(5);
t40 = -mrSges(6,1) * t80 + mrSges(6,2) * t81;
t32 = t36 * mrSges(5,3);
t31 = t36 * mrSges(4,2);
t21 = (-mrSges(6,1) * t60 + mrSges(6,2) * t58) * t37;
t20 = -pkin(4) * t37 + t26;
t12 = t36 * mrSges(6,1) - mrSges(6,3) * t65;
t11 = -mrSges(6,2) * t36 - mrSges(6,3) * t64;
t10 = -t35 * pkin(4) + t18;
t7 = mrSges(6,1) * t64 + mrSges(6,2) * t65;
t6 = Ifges(6,1) * t65 - Ifges(6,4) * t64 + Ifges(6,5) * t36;
t5 = Ifges(6,4) * t65 - Ifges(6,2) * t64 + Ifges(6,6) * t36;
t27 = [0.2e1 * t1 * t24 + 0.2e1 * t10 * t21 + 0.2e1 * t4 * t11 + 0.2e1 * t3 * t12 + 0.2e1 * t2 * t23 + 0.2e1 * t20 * t7 + t32 * t107 + t31 * t106 + 0.2e1 * m(4) * t68 + 0.2e1 * m(5) * (t13 * t22 + t68) + 0.2e1 * m(6) * (t1 * t4 + t10 * t20 + t2 * t3) + (mrSges(5,3) * t108 + t110 * t17 + t75) * t38 + (mrSges(5,2) * t108 + t60 * t5 + t58 * t6 - t18 * t110 + (t15 * t60 + (-t14 - t93) * t58) * qJD(5)) * t37 + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t61) * t79 + (-0.2e1 * pkin(1) * mrSges(3,1) + t105 * t106 + 0.2e1 * pkin(2) * (mrSges(4,1) * t37 + mrSges(4,2) * t38) - 0.2e1 * Ifges(3,4) * t59 + (-Ifges(3,2) + Ifges(3,1)) * t79) * t59) * qJD(2) + (mrSges(4,1) * t106 + mrSges(5,2) * t107 - t26 * t110 + t91 + t92 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t38 + 0.2e1 * (Ifges(4,2) + Ifges(5,3)) * t37) * t35 + (t25 * t110 + ((2 * Ifges(4,1)) + (2 * Ifges(5,2)) + Ifges(6,3)) * t38 + (-0.2e1 * Ifges(4,4) - 0.2e1 * Ifges(5,6) - t69) * t37) * t36; t10 * t44 - t20 * t40 + t50 * t7 + (-mrSges(4,2) + mrSges(5,3)) * t18 + t84 * t17 + (-t37 * mrSges(5,1) + t21) * qJD(4) + (-t41 * t104 + t49 * t12 - t2 * mrSges(6,3) + t6 / 0.2e1) * t60 + (-t5 / 0.2e1 - t42 * t104 + t49 * t11 - t1 * mrSges(6,3)) * t58 + m(5) * (qJD(4) * t26 + t17 * t52 + t18 * t50) + m(6) * (qJD(4) * t20 + t50 * t10 + t109 * t49) + (-t17 * t57 + t18 * t56) * t105 + (t52 * mrSges(5,1) - mrSges(4,3) * t99 + Ifges(6,5) * t102 + Ifges(6,6) * t103 - Ifges(5,4) + Ifges(4,5)) * t36 + (Ifges(3,5) * t61 - Ifges(3,6) * t59 + (-mrSges(3,1) * t61 + mrSges(3,2) * t59) * pkin(6)) * qJD(2) + (t86 / 0.2e1 + Ifges(5,5) - Ifges(4,6) + t88 / 0.2e1 - t50 * mrSges(5,1) - mrSges(4,3) * t100) * t35 + (-t92 / 0.2e1 + t38 * t69 / 0.2e1 - t91 / 0.2e1 + (t102 * t47 + t103 * t46) * t37 + t72 * mrSges(6,3) + (-m(6) * t72 + t67) * t49) * qJD(5); -0.2e1 * t50 * t40 + t41 * t58 - t42 * t60 + (-t86 - t88) * qJD(5) + 0.2e1 * (mrSges(5,3) + t44 + (m(5) + m(6)) * t50) * qJD(4); m(4) * t55 + t60 * t11 - t58 * t12 + t31 - t32 - t84 * t35 + (-t60 * t23 - t58 * t24) * qJD(5) + m(6) * (t1 * t60 - t2 * t58 + (-t3 * t60 - t4 * t58) * qJD(5)) + m(5) * t13; 0; 0; t36 * mrSges(5,1) + t58 * t11 + t60 * t12 + t67 * qJD(5) + m(6) * (-qJD(5) * t72 + t109) + m(5) * t17; 0; 0; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 - Ifges(6,6) * t78 + t75; (-t44 * t49 + t69) * qJD(5); t40; -t44 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t27(1), t27(2), t27(4), t27(7), t27(11); t27(2), t27(3), t27(5), t27(8), t27(12); t27(4), t27(5), t27(6), t27(9), t27(13); t27(7), t27(8), t27(9), t27(10), t27(14); t27(11), t27(12), t27(13), t27(14), t27(15);];
Mq = res;
