% Calculate time derivative of joint inertia matrix for
% S5RPPRR12
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR12_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR12_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:57
% EndTime: 2019-12-31 18:07:00
% DurationCPUTime: 0.71s
% Computational Cost: add. (1062->169), mult. (2188->261), div. (0->0), fcn. (1938->6), ass. (0->79)
t95 = 2 * qJD(2);
t49 = sin(pkin(8));
t50 = cos(pkin(8));
t87 = sin(qJ(4));
t88 = cos(qJ(4));
t27 = t87 * t49 - t88 * t50;
t61 = (t49 ^ 2 + t50 ^ 2) * qJD(3);
t52 = sin(qJ(5));
t53 = cos(qJ(5));
t34 = -mrSges(6,1) * t53 + t52 * mrSges(6,2);
t94 = m(6) * pkin(4) + mrSges(5,1) - t34;
t93 = -2 * mrSges(5,3);
t51 = -pkin(1) - qJ(3);
t89 = -pkin(6) + t51;
t32 = t89 * t49;
t33 = t89 * t50;
t19 = t87 * t32 - t88 * t33;
t92 = 0.2e1 * t19;
t91 = t27 / 0.2e1;
t84 = Ifges(6,4) * t53;
t36 = Ifges(6,1) * t52 + t84;
t90 = t36 / 0.2e1;
t86 = mrSges(6,3) * t27;
t85 = Ifges(6,4) * t52;
t62 = qJD(4) * t87;
t63 = qJD(4) * t88;
t25 = -t49 * t62 + t50 * t63;
t83 = Ifges(6,5) * t25;
t82 = Ifges(6,6) * t25;
t28 = t88 * t49 + t87 * t50;
t81 = Ifges(6,6) * t28;
t80 = Ifges(6,6) * t52;
t20 = t88 * t32 + t87 * t33;
t13 = -t27 * qJD(3) + t20 * qJD(4);
t79 = t13 * t19;
t24 = -t49 * t63 - t50 * t62;
t78 = t24 * t27;
t77 = t53 * t24;
t76 = Ifges(6,5) * t77 + Ifges(6,3) * t25;
t41 = t49 * pkin(3) + qJ(2);
t15 = pkin(4) * t28 + pkin(7) * t27 + t41;
t6 = t15 * t53 - t52 * t20;
t74 = qJD(5) * t6;
t7 = t52 * t15 + t20 * t53;
t73 = qJD(5) * t7;
t72 = qJD(5) * t52;
t71 = qJD(5) * t53;
t70 = t25 * t93;
t69 = t27 * t72;
t66 = (t52 ^ 2 + t53 ^ 2) * t25;
t65 = t25 * mrSges(5,1) + t24 * mrSges(5,2);
t64 = -(2 * Ifges(5,4)) - t80;
t60 = -t52 * t6 + t53 * t7;
t59 = mrSges(6,1) * t52 + mrSges(6,2) * t53;
t58 = Ifges(6,1) * t53 - t85;
t57 = -Ifges(6,2) * t52 + t84;
t56 = t27 * t13 - t24 * t19;
t55 = -t52 * t24 + t27 * t71;
t54 = t69 + t77;
t43 = Ifges(6,5) * t71;
t35 = Ifges(6,2) * t53 + t85;
t31 = t58 * qJD(5);
t30 = t57 * qJD(5);
t29 = t59 * qJD(5);
t18 = t28 * mrSges(6,1) + t53 * t86;
t17 = -mrSges(6,2) * t28 + t52 * t86;
t16 = t59 * t27;
t14 = pkin(4) * t25 - pkin(7) * t24 + qJD(2);
t12 = -t28 * qJD(3) - t19 * qJD(4);
t11 = Ifges(6,5) * t28 - t58 * t27;
t10 = -t57 * t27 + t81;
t9 = -t25 * mrSges(6,2) + t55 * mrSges(6,3);
t8 = t25 * mrSges(6,1) - t54 * mrSges(6,3);
t5 = -t55 * mrSges(6,1) + t54 * mrSges(6,2);
t4 = t54 * Ifges(6,1) + t55 * Ifges(6,4) + t83;
t3 = t54 * Ifges(6,4) + t55 * Ifges(6,2) + t82;
t2 = -t52 * t12 + t14 * t53 - t73;
t1 = t12 * t53 + t52 * t14 + t74;
t21 = [t20 * t70 + 0.2e1 * t41 * t65 - 0.2e1 * t13 * t16 + 0.2e1 * t1 * t17 + 0.2e1 * t2 * t18 + t5 * t92 + 0.2e1 * t6 * t8 + 0.2e1 * t7 * t9 + (mrSges(5,3) * t92 - t52 * t10 + t53 * t11) * t24 + 0.2e1 * m(5) * (qJD(2) * t41 + t12 * t20 + t79) + 0.2e1 * m(4) * (qJ(2) * qJD(2) - t51 * t61) + 0.2e1 * m(6) * (t1 * t7 + t2 * t6 + t79) + (mrSges(5,1) * t95 + t12 * t93 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t25 + t64 * t24 + t76) * t28 + (-0.2e1 * qJD(2) * mrSges(5,2) + t13 * t93 - 0.2e1 * Ifges(5,1) * t24 + t52 * t3 - t53 * t4 + (-Ifges(6,5) * t53 - t64) * t25 + (t52 * t11 + t53 * t10 + t28 * (Ifges(6,5) * t52 + Ifges(6,6) * t53)) * qJD(5)) * t27 + (m(3) * qJ(2) + t49 * mrSges(4,1) + t50 * mrSges(4,2) + mrSges(3,3)) * t95 + 0.2e1 * mrSges(4,3) * t61; t27 * t5 + (t53 * t17 - t52 * t18) * t25 + (0.2e1 * t27 * mrSges(5,3) + t16) * t24 + (t70 - t52 * t8 + t53 * t9 + (-t17 * t52 - t18 * t53) * qJD(5)) * t28 + m(6) * (t60 * t25 + (t1 * t53 - t2 * t52 + (-t52 * t7 - t53 * t6) * qJD(5)) * t28 + t56) + m(5) * (t12 * t28 + t20 * t25 + t56) - m(4) * t61; 0.2e1 * m(5) * (t25 * t28 - t78) + 0.2e1 * m(6) * (t28 * t66 - t78); m(6) * (t60 * qJD(5) + t52 * t1 + t2 * t53) + t17 * t71 + t52 * t9 - t18 * t72 + t53 * t8 + (m(5) + m(4)) * qJD(2) + t65; 0; 0; t28 * t43 / 0.2e1 + t19 * t29 + Ifges(5,5) * t24 - Ifges(5,6) * t25 - pkin(4) * t5 - t12 * mrSges(5,2) - t94 * t13 + (t30 * t91 - t24 * t35 / 0.2e1 - t2 * mrSges(6,3) + t4 / 0.2e1 + t83 / 0.2e1 + (t27 * t90 - t7 * mrSges(6,3) - t10 / 0.2e1 - t81 / 0.2e1) * qJD(5) + (-qJD(5) * t17 - t8 + m(6) * (-t2 - t73)) * pkin(7)) * t52 + (-t27 * t31 / 0.2e1 + t24 * t90 + t1 * mrSges(6,3) + t82 / 0.2e1 + t3 / 0.2e1 + (t35 * t91 - t6 * mrSges(6,3) + t11 / 0.2e1) * qJD(5) + (m(6) * (t1 - t74) + t9 - qJD(5) * t18) * pkin(7)) * t53; -t25 * mrSges(5,2) + t27 * t29 + (m(6) * pkin(7) + mrSges(6,3)) * t66 + t94 * t24; 0; -0.2e1 * pkin(4) * t29 + t30 * t53 + t52 * t31 + (-t52 * t35 + t53 * t36) * qJD(5); t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t69 + t55 * Ifges(6,6) + t76; (-t53 * t25 + t28 * t72) * mrSges(6,2) + (-t25 * t52 - t28 * t71) * mrSges(6,1); -t29; t43 + (t34 * pkin(7) - t80) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t21(1), t21(2), t21(4), t21(7), t21(11); t21(2), t21(3), t21(5), t21(8), t21(12); t21(4), t21(5), t21(6), t21(9), t21(13); t21(7), t21(8), t21(9), t21(10), t21(14); t21(11), t21(12), t21(13), t21(14), t21(15);];
Mq = res;
