% Calculate time derivative of joint inertia matrix for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR7_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR7_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:55
% EndTime: 2019-12-31 18:18:57
% DurationCPUTime: 0.85s
% Computational Cost: add. (1128->176), mult. (2442->281), div. (0->0), fcn. (2110->8), ass. (0->88)
t54 = cos(qJ(3));
t81 = cos(pkin(9));
t49 = sin(pkin(9));
t52 = sin(qJ(3));
t90 = t49 * t52;
t56 = t54 * t81 - t90;
t31 = t56 * qJD(3);
t51 = sin(qJ(5));
t65 = t81 * t52;
t34 = t49 * t54 + t65;
t53 = cos(qJ(5));
t76 = qJD(5) * t53;
t72 = t34 * t76;
t58 = -t51 * t31 - t72;
t77 = qJD(5) * t51;
t73 = t34 * t77;
t85 = t53 * t31;
t57 = t73 - t85;
t70 = (t51 ^ 2 + t53 ^ 2) * t31;
t44 = sin(pkin(8)) * pkin(1) + pkin(6);
t82 = qJ(4) + t44;
t64 = qJD(3) * t82;
t22 = qJD(4) * t54 - t52 * t64;
t55 = -t52 * qJD(4) - t54 * t64;
t10 = t22 * t49 - t55 * t81;
t102 = 0.2e1 * t10;
t71 = -cos(pkin(8)) * pkin(1) - pkin(2);
t38 = -pkin(3) * t54 + t71;
t101 = 0.2e1 * t38;
t100 = m(5) * pkin(3);
t99 = t53 / 0.2e1;
t98 = pkin(3) * t49;
t97 = Ifges(6,4) * t51;
t96 = Ifges(6,4) * t53;
t95 = Ifges(6,6) * t51;
t32 = t82 * t54;
t16 = t32 * t49 + t65 * t82;
t94 = t10 * t16;
t30 = t34 * qJD(3);
t93 = t30 * mrSges(5,3);
t92 = t30 * t56;
t89 = t51 * mrSges(6,3);
t40 = Ifges(6,2) * t53 + t97;
t87 = t51 * t40;
t86 = t53 * mrSges(6,3);
t41 = Ifges(6,1) * t51 + t96;
t84 = t53 * t41;
t83 = Ifges(6,5) * t85 + Ifges(6,3) * t30;
t80 = qJD(3) * t52;
t79 = qJD(3) * t54;
t78 = qJD(5) * t34;
t75 = 0.2e1 * t54;
t74 = pkin(3) * t80;
t69 = t81 * pkin(3);
t68 = -t77 / 0.2e1;
t67 = t30 * mrSges(5,1) + t31 * mrSges(5,2);
t66 = -(2 * Ifges(5,4)) - t95;
t39 = -mrSges(6,1) * t53 + mrSges(6,2) * t51;
t63 = mrSges(6,1) * t51 + mrSges(6,2) * t53;
t62 = Ifges(6,1) * t53 - t97;
t61 = -Ifges(6,2) * t51 + t96;
t60 = Ifges(6,5) * t51 + Ifges(6,6) * t53;
t59 = -t10 * t56 + t16 * t30;
t15 = -pkin(4) * t56 - t34 * pkin(7) + t38;
t17 = t32 * t81 - t82 * t90;
t5 = t15 * t53 - t17 * t51;
t6 = t15 * t51 + t17 * t53;
t46 = Ifges(6,5) * t76;
t45 = -t69 - pkin(4);
t43 = pkin(7) + t98;
t37 = t62 * qJD(5);
t36 = t61 * qJD(5);
t35 = t63 * qJD(5);
t20 = -mrSges(6,1) * t56 - t34 * t86;
t19 = mrSges(6,2) * t56 - t34 * t89;
t18 = t63 * t34;
t14 = pkin(4) * t30 - pkin(7) * t31 + t74;
t13 = -Ifges(6,5) * t56 + t34 * t62;
t12 = -Ifges(6,6) * t56 + t34 * t61;
t11 = t22 * t81 + t49 * t55;
t9 = -mrSges(6,2) * t30 + mrSges(6,3) * t58;
t8 = mrSges(6,1) * t30 + mrSges(6,3) * t57;
t7 = t58 * mrSges(6,1) + t57 * mrSges(6,2);
t4 = -Ifges(6,1) * t57 + Ifges(6,4) * t58 + t30 * Ifges(6,5);
t3 = -Ifges(6,4) * t57 + Ifges(6,2) * t58 + t30 * Ifges(6,6);
t2 = -t6 * qJD(5) - t11 * t51 + t14 * t53;
t1 = t5 * qJD(5) + t11 * t53 + t14 * t51;
t21 = [t67 * t101 - 0.2e1 * t17 * t93 - 0.2e1 * t16 * t7 + t18 * t102 + 0.2e1 * t1 * t19 + 0.2e1 * t2 * t20 + 0.2e1 * t5 * t8 + 0.2e1 * t6 * t9 + 0.2e1 * m(5) * (t11 * t17 + t94) + 0.2e1 * m(6) * (t1 * t6 + t2 * t5 + t94) + (0.2e1 * mrSges(5,3) * t16 - t51 * t12 + t53 * t13) * t31 - (-0.2e1 * t11 * mrSges(5,3) + t66 * t31 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t30 + t83) * t56 + (mrSges(5,3) * t102 + 0.2e1 * Ifges(5,1) * t31 - t51 * t3 + t53 * t4 + (Ifges(6,5) * t53 + t66) * t30 + (-t53 * t12 - t51 * t13 + t56 * t60) * qJD(5)) * t34 + ((t71 * mrSges(4,2) + Ifges(4,4) * t54) * t75 + (0.2e1 * t71 * mrSges(4,1) + t100 * t101 + 0.2e1 * pkin(3) * (-mrSges(5,1) * t56 + t34 * mrSges(5,2)) - 0.2e1 * Ifges(4,4) * t52 + (Ifges(4,1) - Ifges(4,2)) * t75) * t52) * qJD(3); t30 * t18 + t56 * t7 + m(6) * t59 + m(5) * (t11 * t34 + t17 * t31 + t59) + (t31 * t19 + t34 * t9 - t20 * t78 + m(6) * (t1 * t34 + t31 * t6 - t5 * t78)) * t53 + (-t19 * t78 - t31 * t20 - t34 * t8 + m(6) * (-t2 * t34 - t31 * t5 - t6 * t78)) * t51; 0.2e1 * m(6) * (t34 * t70 - t92) + 0.2e1 * m(5) * (t31 * t34 - t92); -t93 * t98 + t12 * t68 + t51 * t4 / 0.2e1 + t16 * t35 - t45 * t7 - t2 * t89 - Ifges(4,6) * t80 + t13 * t76 / 0.2e1 - t56 * (-Ifges(6,6) * t77 + t46) / 0.2e1 - t40 * t72 / 0.2e1 + t3 * t99 + Ifges(4,5) * t79 + t1 * t86 + (-mrSges(4,1) * t79 + mrSges(4,2) * t80) * t44 + (t60 / 0.2e1 - Ifges(5,6)) * t30 + (t49 * t100 - mrSges(5,2)) * t11 + (-t5 * t76 - t6 * t77) * mrSges(6,3) + (t37 * t99 - t51 * t36 / 0.2e1 + t41 * t68) * t34 + (m(6) * t45 - t81 * t100 - mrSges(5,1) + t39) * t10 + (-t20 * t76 - t19 * t77 + m(6) * (t1 * t53 - t2 * t51 + (-t5 * t53 - t51 * t6) * qJD(5)) + t53 * t9 - t51 * t8) * t43 + (-mrSges(5,3) * t69 + Ifges(5,5) + t84 / 0.2e1 - t87 / 0.2e1) * t31; m(6) * (t45 * t30 + t43 * t70) + t30 * t39 - t56 * t35 - mrSges(4,2) * t79 - mrSges(4,1) * t80 + (-t30 * t81 + t31 * t49) * t100 - t67 + mrSges(6,3) * t70; 0.2e1 * t35 * t45 + t36 * t53 + t37 * t51 + (t84 - t87) * qJD(5); -t20 * t77 + t53 * t8 + m(6) * (t1 * t51 + t2 * t53 + (-t5 * t51 + t53 * t6) * qJD(5)) + t19 * t76 + t51 * t9 + m(5) * t74 + t67; 0; 0; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 - Ifges(6,5) * t73 + Ifges(6,6) * t58 + t83; t7; t46 + (t39 * t43 - t95) * qJD(5); -t35; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t21(1), t21(2), t21(4), t21(7), t21(11); t21(2), t21(3), t21(5), t21(8), t21(12); t21(4), t21(5), t21(6), t21(9), t21(13); t21(7), t21(8), t21(9), t21(10), t21(14); t21(11), t21(12), t21(13), t21(14), t21(15);];
Mq = res;
