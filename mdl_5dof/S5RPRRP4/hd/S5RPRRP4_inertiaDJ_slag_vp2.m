% Calculate time derivative of joint inertia matrix for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:51
% EndTime: 2022-01-23 09:31:55
% DurationCPUTime: 1.16s
% Computational Cost: add. (1311->163), mult. (3178->257), div. (0->0), fcn. (2721->6), ass. (0->85)
t70 = sin(qJ(4));
t71 = sin(qJ(3));
t72 = cos(qJ(4));
t73 = cos(qJ(3));
t56 = t70 * t73 + t72 * t71;
t68 = sin(pkin(8));
t45 = t56 * t68;
t99 = t71 * t68;
t54 = pkin(3) * t99 + t68 * qJ(2);
t30 = t45 * pkin(4) + t54;
t111 = -0.2e1 * t30;
t110 = 0.2e1 * t68;
t109 = mrSges(5,1) + mrSges(6,1);
t108 = -Ifges(5,5) - Ifges(6,5);
t107 = -Ifges(5,6) - Ifges(6,6);
t106 = pkin(3) * t72;
t76 = t70 * t71 - t72 * t73;
t46 = t76 * t68;
t69 = cos(pkin(8));
t57 = -pkin(2) * t69 - pkin(6) * t68 - pkin(1);
t51 = t73 * t57;
t98 = t73 * t68;
t86 = pkin(7) * t98;
t92 = qJ(2) * t71;
t29 = -t86 + t51 + (-pkin(3) - t92) * t69;
t60 = t73 * t69 * qJ(2);
t44 = t71 * t57 + t60;
t38 = -pkin(7) * t99 + t44;
t11 = t70 * t29 + t72 * t38;
t104 = pkin(3) * qJD(4);
t103 = qJD(3) + qJD(4);
t102 = 0.2e1 * t69;
t101 = m(6) * pkin(4);
t22 = t103 * t46;
t100 = t22 * t70;
t97 = -mrSges(5,2) - mrSges(6,2);
t36 = t103 * t76;
t88 = qJD(4) * t72;
t96 = (-t36 * t70 + t56 * t88) * pkin(3);
t39 = mrSges(6,2) * t69 - t45 * mrSges(6,3);
t40 = mrSges(5,2) * t69 - t45 * mrSges(5,3);
t95 = t39 + t40;
t41 = -mrSges(6,1) * t69 + t46 * mrSges(6,3);
t42 = -mrSges(5,1) * t69 + t46 * mrSges(5,3);
t94 = t41 + t42;
t90 = qJD(3) * t73;
t91 = qJD(2) * t69;
t93 = t57 * t90 + t73 * t91;
t64 = t68 * qJD(2);
t49 = t68 * pkin(3) * t90 + t64;
t89 = qJD(4) * t70;
t87 = qJ(2) * qJD(2);
t85 = t69 * t92;
t84 = t76 * t89;
t83 = t71 * t91;
t82 = t97 * t72;
t10 = t72 * t29 - t38 * t70;
t37 = t103 * t56;
t21 = t37 * t68;
t81 = t107 * t22 - t108 * t21;
t80 = 2 * Ifges(5,4) + 2 * Ifges(6,4);
t23 = (-t85 - t86) * qJD(3) + t93;
t24 = -t83 + (-t60 + (pkin(7) * t68 - t57) * t71) * qJD(3);
t6 = -t11 * qJD(4) - t23 * t70 + t72 * t24;
t3 = t21 * qJ(5) + t46 * qJD(5) + t6;
t79 = m(6) * t3 + t21 * mrSges(6,3);
t78 = mrSges(4,1) * t71 + mrSges(4,2) * t73;
t75 = -t109 * t37 - t36 * t97;
t5 = t72 * t23 + t70 * t24 + t29 * t88 - t38 * t89;
t2 = t22 * qJ(5) - t45 * qJD(5) + t5;
t74 = t6 * mrSges(5,1) + t3 * mrSges(6,1) - t5 * mrSges(5,2) - t2 * mrSges(6,2) - t81;
t67 = t69 ^ 2;
t66 = t68 ^ 2;
t63 = pkin(4) + t106;
t62 = t66 * t87;
t53 = -t69 * mrSges(4,1) - mrSges(4,3) * t98;
t52 = mrSges(4,2) * t69 - mrSges(4,3) * t99;
t43 = t51 - t85;
t32 = -t44 * qJD(3) - t83;
t31 = -qJD(3) * t85 + t93;
t16 = t21 * mrSges(6,2);
t12 = -t22 * pkin(4) + t49;
t8 = -qJ(5) * t45 + t11;
t7 = -pkin(4) * t69 + t46 * qJ(5) + t10;
t1 = [0.2e1 * t12 * (t45 * mrSges(6,1) - t46 * mrSges(6,2)) + 0.2e1 * t49 * (t45 * mrSges(5,1) - t46 * mrSges(5,2)) + 0.2e1 * t31 * t52 + 0.2e1 * t32 * t53 + 0.2e1 * t2 * t39 + 0.2e1 * t5 * t40 + 0.2e1 * t3 * t41 + 0.2e1 * t6 * t42 + t16 * t111 + t81 * t69 + 0.2e1 * (t66 + t67) * qJD(2) * mrSges(3,3) + (-0.2e1 * t54 * mrSges(5,1) + mrSges(6,1) * t111 + 0.2e1 * t11 * mrSges(5,3) + 0.2e1 * t8 * mrSges(6,3) + t107 * t69 - t46 * t80 - 0.2e1 * (Ifges(5,2) + Ifges(6,2)) * t45) * t22 - (0.2e1 * t54 * mrSges(5,2) - 0.2e1 * t10 * mrSges(5,3) - 0.2e1 * t7 * mrSges(6,3) + t108 * t69 - 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t46 - t45 * t80) * t21 + 0.2e1 * m(4) * (t44 * t31 + t43 * t32 + t62) + 0.2e1 * m(3) * (t67 * t87 + t62) + 0.2e1 * m(5) * (t10 * t6 + t11 * t5 + t49 * t54) + 0.2e1 * m(6) * (t12 * t30 + t2 * t8 + t3 * t7) + (0.2e1 * t78 * t64 + ((-0.2e1 * t44 * mrSges(4,3) + Ifges(4,6) * t102 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t73) * t110) * t73 + (0.2e1 * t43 * mrSges(4,3) + Ifges(4,5) * t102 + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t71 + (-Ifges(4,1) + Ifges(4,2)) * t73) * t110) * t71) * qJD(3)) * t68; -t94 * t37 - t95 * t36 + (t73 * t52 - t71 * t53) * qJD(3) + m(6) * (t2 * t56 - t3 * t76 - t36 * t8 - t37 * t7) + m(5) * (-t10 * t37 - t11 * t36 + t5 * t56 - t6 * t76) + m(4) * (t71 * t31 + t32 * t73 + (-t43 * t71 + t44 * t73) * qJD(3)) + (mrSges(6,3) + mrSges(5,3)) * (-t21 * t76 + t22 * t56); 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (-t36 * t56 + t37 * t76); t32 * mrSges(4,1) - t31 * mrSges(4,2) + t79 * t63 + (-Ifges(4,5) * t71 - Ifges(4,6) * t73) * t68 * qJD(3) + (mrSges(6,3) * t100 + (t21 * t72 + t100) * mrSges(5,3) + (-t70 * t94 + t72 * t95) * qJD(4) + m(6) * (t2 * t70 - t7 * t89 + t8 * t88) + m(5) * (-t10 * t89 + t11 * t88 + t5 * t70 + t6 * t72)) * pkin(3) + t74; -t78 * qJD(3) + m(5) * ((-t37 * t72 + t84) * pkin(3) + t96) + m(6) * (pkin(3) * t84 - t37 * t63 + t96) + t75; 0.2e1 * (t82 + ((-t63 + t106) * m(6) - t109) * t70) * t104; pkin(4) * t79 + t74; -t101 * t37 + t75; (t82 + (-t101 - t109) * t70) * t104; 0; m(6) * t12 - t22 * mrSges(6,1) - t16; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
