% Calculate time derivative of joint inertia matrix for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR11_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR11_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:14
% EndTime: 2024-09-27 21:45:15
% DurationCPUTime: 0.92s
% Computational Cost: add. (2186->194), mult. (6590->325), div. (0->0), fcn. (6005->8), ass. (0->93)
t73 = sin(qJ(3));
t76 = cos(qJ(3));
t69 = sin(pkin(5));
t88 = t69 * qJD(3);
t51 = (mrSges(4,1) * t73 + mrSges(4,2) * t76) * t88;
t110 = qJD(3) + qJD(4);
t72 = sin(qJ(4));
t75 = cos(qJ(4));
t49 = (-t72 * t73 + t75 * t76) * t69;
t30 = t110 * t49;
t50 = (t72 * t76 + t73 * t75) * t69;
t31 = t110 * t50;
t17 = t31 * mrSges(5,1) + t30 * mrSges(5,2);
t71 = sin(qJ(5));
t74 = cos(qJ(5));
t25 = t49 * t74 - t50 * t71;
t12 = t25 * qJD(5) + t30 * t74 - t31 * t71;
t26 = t49 * t71 + t50 * t74;
t13 = -t26 * qJD(5) - t30 * t71 - t31 * t74;
t4 = -t13 * mrSges(6,1) + t12 * mrSges(6,2);
t79 = -t17 - t4;
t112 = t51 - t79;
t111 = Ifges(5,5) * t30 - Ifges(5,6) * t31;
t70 = cos(pkin(5));
t100 = pkin(3) * t70;
t101 = pkin(2) * t70;
t66 = t76 * t101;
t87 = t69 * (-pkin(7) - pkin(8));
t84 = t73 * t87;
t42 = t66 + t84 + t100;
t65 = t73 * t101;
t97 = t69 * t76;
t57 = pkin(7) * t97 + t65;
t48 = pkin(8) * t97 + t57;
t21 = t72 * t42 + t75 * t48;
t53 = t57 * qJD(3);
t109 = 2 * m(6);
t108 = -2 * mrSges(4,3);
t107 = 0.2e1 * t25;
t106 = 0.2e1 * t26;
t105 = 0.2e1 * t49;
t104 = 0.2e1 * t50;
t61 = (-pkin(3) * t76 - pkin(2)) * t69;
t103 = 0.2e1 * t61;
t102 = m(5) * pkin(3);
t67 = pkin(3) * t75 + pkin(4);
t89 = qJD(5) * t74;
t90 = qJD(5) * t71;
t96 = t71 * t72;
t33 = t67 * t89 + (-t72 * t90 + (t74 * t75 - t96) * qJD(4)) * pkin(3);
t99 = t33 * mrSges(6,2);
t98 = t69 * t73;
t95 = t72 * t74;
t94 = Ifges(6,5) * t12 + Ifges(6,6) * t13;
t93 = qJD(3) * t73;
t92 = qJD(4) * t72;
t91 = qJD(4) * t75;
t86 = t73 * t88;
t34 = -t67 * t90 + (-t72 * t89 + (-t71 * t75 - t95) * qJD(4)) * pkin(3);
t32 = t34 * mrSges(6,1);
t85 = t32 - t99;
t20 = t75 * t42 - t48 * t72;
t18 = pkin(4) * t70 - pkin(9) * t50 + t20;
t19 = pkin(9) * t49 + t21;
t5 = t18 * t74 - t19 * t71;
t6 = t18 * t71 + t19 * t74;
t63 = qJD(3) * t66;
t43 = qJD(3) * t84 + t63;
t44 = (t76 * t87 - t65) * qJD(3);
t15 = t42 * t91 + t75 * t43 + t72 * t44 - t48 * t92;
t7 = -pkin(9) * t31 + t15;
t16 = -t21 * qJD(4) - t43 * t72 + t75 * t44;
t8 = -pkin(9) * t30 + t16;
t2 = t5 * qJD(5) + t7 * t74 + t71 * t8;
t3 = -t6 * qJD(5) - t7 * t71 + t74 * t8;
t80 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + t94;
t78 = (-mrSges(5,1) * t72 - mrSges(5,2) * t75) * qJD(4) * pkin(3);
t77 = t16 * mrSges(5,1) - t15 * mrSges(5,2) + t111 + t80;
t62 = Ifges(4,5) * t76 * t88;
t60 = -t70 * mrSges(4,2) + mrSges(4,3) * t97;
t59 = mrSges(4,1) * t70 - mrSges(4,3) * t98;
t58 = (-t71 * mrSges(6,1) - t74 * mrSges(6,2)) * qJD(5) * pkin(4);
t56 = -pkin(7) * t98 + t66;
t55 = pkin(3) * t95 + t67 * t71;
t54 = -pkin(3) * t96 + t67 * t74;
t52 = -pkin(7) * t86 + t63;
t46 = mrSges(5,1) * t70 - t50 * mrSges(5,3);
t45 = -mrSges(5,2) * t70 + t49 * mrSges(5,3);
t35 = -t49 * pkin(4) + t61;
t24 = mrSges(6,1) * t70 - t26 * mrSges(6,3);
t23 = -mrSges(6,2) * t70 + t25 * mrSges(6,3);
t22 = pkin(3) * t86 + pkin(4) * t31;
t1 = [0.2e1 * m(6) * (t12 * t26 + t13 * t25) + 0.2e1 * m(5) * (t30 * t50 - t31 * t49); t12 * t23 + t13 * t24 + t30 * t45 - t31 * t46 + (-t12 * t25 + t13 * t26) * mrSges(6,3) + (-t30 * t49 - t31 * t50) * mrSges(5,3) + t112 * t70 + (-t59 * t73 + t60 * t76 + (-t73 ^ 2 - t76 ^ 2) * mrSges(4,3) * t69) * t88 + m(6) * (t12 * t6 + t13 * t5 + t2 * t26 + t22 * t70 + t25 * t3) + m(4) * (t52 * t73 - t56 * t93) * t69 + (t93 * t100 * t69 + t15 * t50 + t16 * t49 - t20 * t31 + t21 * t30) * m(5); -0.2e1 * t53 * t59 + 0.2e1 * t52 * t60 + t17 * t103 + 0.2e1 * t15 * t45 + 0.2e1 * t16 * t46 + 0.2e1 * t35 * t4 + 0.2e1 * t2 * t23 + 0.2e1 * t3 * t24 + 0.2e1 * t22 * (-mrSges(6,1) * t25 + mrSges(6,2) * t26) + (t62 + t94 + t111) * t70 - (0.2e1 * t21 * mrSges(5,3) + Ifges(5,4) * t104 + Ifges(5,2) * t105 + Ifges(5,6) * t70) * t31 + (-0.2e1 * t20 * mrSges(5,3) + Ifges(5,1) * t104 + Ifges(5,4) * t105 + Ifges(5,5) * t70) * t30 + (0.2e1 * t6 * mrSges(6,3) + Ifges(6,4) * t106 + Ifges(6,2) * t107 + Ifges(6,6) * t70) * t13 + (-0.2e1 * t5 * mrSges(6,3) + Ifges(6,1) * t106 + Ifges(6,4) * t107 + Ifges(6,5) * t70) * t12 + 0.2e1 * m(5) * (t15 * t21 + t16 * t20) + 0.2e1 * m(4) * (t52 * t57 - t53 * t56) + (t2 * t6 + t22 * t35 + t3 * t5) * t109 + (-0.2e1 * pkin(2) * t51 + ((0.2e1 * Ifges(4,4) * t97 + Ifges(4,5) * t70 + t56 * t108) * t76 + (-0.2e1 * Ifges(4,6) * t70 + t102 * t103 + 0.2e1 * pkin(3) * (-mrSges(5,1) * t49 + mrSges(5,2) * t50) + t57 * t108 + 0.2e1 * (-Ifges(4,4) * t73 + (Ifges(4,1) - Ifges(4,2)) * t76) * t69) * t73) * qJD(3)) * t69; m(6) * (t12 * t55 + t13 * t54 + t25 * t34 + t26 * t33) + (t30 * t72 - t31 * t75 + (-t49 * t72 + t50 * t75) * qJD(4)) * t102 - t112; (-t54 * t12 + t55 * t13) * mrSges(6,3) + m(6) * (t2 * t55 + t3 * t54 + t33 * t6 + t34 * t5) + t77 + t62 - t52 * mrSges(4,2) - t53 * mrSges(4,1) + t33 * t23 + t34 * t24 - Ifges(4,6) * t86 + (t45 * t91 - t46 * t92 + m(5) * (t15 * t72 + t16 * t75 - t20 * t92 + t21 * t91) + (-t75 * t30 - t72 * t31) * mrSges(5,3)) * pkin(3); -0.2e1 * t99 + 0.2e1 * t32 + (t33 * t55 + t34 * t54) * t109 + 0.2e1 * t78; m(6) * (t12 * t71 + t13 * t74 + (-t25 * t71 + t26 * t74) * qJD(5)) * pkin(4) + t79; (t23 * t89 - t24 * t90 + m(6) * (t2 * t71 + t3 * t74 - t5 * t90 + t6 * t89) + (-t12 * t74 + t13 * t71) * mrSges(6,3)) * pkin(4) + t77; t78 + (m(6) * (t33 * t71 + t34 * t74 - t54 * t90 + t55 * t89) - mrSges(6,2) * t89 - mrSges(6,1) * t90) * pkin(4) + t85; 0.2e1 * t58; -t4; t80; t85; t58; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
