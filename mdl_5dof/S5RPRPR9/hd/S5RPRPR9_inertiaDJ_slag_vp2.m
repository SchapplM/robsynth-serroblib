% Calculate time derivative of joint inertia matrix for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR9_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR9_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:43
% EndTime: 2019-12-31 18:23:45
% DurationCPUTime: 0.74s
% Computational Cost: add. (580->158), mult. (1260->237), div. (0->0), fcn. (836->6), ass. (0->83)
t39 = sin(qJ(5));
t40 = sin(qJ(3));
t42 = cos(qJ(3));
t63 = qJD(5) * t42;
t41 = cos(qJ(5));
t67 = qJD(3) * t41;
t45 = t39 * t63 + t40 * t67;
t57 = -cos(pkin(8)) * pkin(1) - pkin(2);
t93 = 0.2e1 * t57;
t92 = -mrSges(4,1) + mrSges(5,2);
t26 = t39 * mrSges(6,1) + t41 * mrSges(6,2);
t91 = mrSges(5,3) + t26;
t68 = qJD(3) * t40;
t59 = t39 * t68;
t60 = t41 * t63;
t90 = -t59 + t60;
t89 = 2 * m(6);
t70 = t40 * qJ(4);
t46 = t57 - t70;
t84 = t42 * pkin(3);
t17 = t46 - t84;
t88 = -0.2e1 * t17;
t87 = -t39 / 0.2e1;
t86 = -t41 / 0.2e1;
t85 = t41 / 0.2e1;
t43 = -pkin(3) - pkin(7);
t33 = sin(pkin(8)) * pkin(1) + pkin(6);
t83 = pkin(4) + t33;
t82 = mrSges(6,3) * t42;
t81 = Ifges(6,4) * t39;
t80 = Ifges(6,4) * t41;
t79 = Ifges(6,5) * t40;
t49 = Ifges(6,1) * t39 + t80;
t13 = -t49 * t42 + t79;
t78 = t39 * t13;
t28 = Ifges(6,1) * t41 - t81;
t77 = t39 * t28;
t76 = t39 * t43;
t75 = t40 * mrSges(5,3);
t48 = Ifges(6,2) * t41 + t81;
t12 = Ifges(6,6) * t40 - t48 * t42;
t74 = t41 * t12;
t27 = -Ifges(6,2) * t39 + t80;
t73 = t41 * t27;
t72 = t41 * t43;
t71 = t39 ^ 2 + t41 ^ 2;
t20 = t83 * t42;
t15 = qJD(3) * t20;
t69 = qJD(3) * t39;
t66 = qJD(3) * t42;
t65 = qJD(5) * t39;
t64 = qJD(5) * t41;
t62 = t40 * qJD(4);
t56 = m(6) * t71;
t55 = m(5) * t33 + mrSges(5,1);
t54 = Ifges(6,5) * t59 + t45 * Ifges(6,6) + Ifges(6,3) * t66;
t53 = pkin(3) * t68 - t62;
t10 = (pkin(7) * t40 - qJ(4) * t42) * qJD(3) + t53;
t11 = t43 * t42 + t46;
t19 = t83 * t40;
t3 = -t39 * t11 + t41 * t19;
t1 = t3 * qJD(5) + t41 * t10 + t39 * t15;
t4 = t41 * t11 + t39 * t19;
t2 = -t4 * qJD(5) - t39 * t10 + t41 * t15;
t52 = t39 * t1 + t41 * t2;
t51 = t3 * t39 - t4 * t41;
t50 = mrSges(6,1) * t41 - mrSges(6,2) * t39;
t47 = -Ifges(6,5) * t39 - Ifges(6,6) * t41;
t24 = t40 * mrSges(6,1) + t39 * t82;
t25 = -t40 * mrSges(6,2) - t41 * t82;
t8 = -mrSges(6,2) * t66 + t45 * mrSges(6,3);
t9 = mrSges(6,1) * t66 + mrSges(6,3) * t90;
t44 = -t24 * t65 + t25 * t64 + t39 * t8 + t41 * t9;
t23 = t49 * qJD(5);
t22 = t48 * qJD(5);
t21 = t50 * qJD(5);
t18 = t50 * t42;
t16 = qJ(4) * t66 - t53;
t14 = t83 * t68;
t7 = t45 * mrSges(6,1) + mrSges(6,2) * t90;
t6 = -t28 * t63 + (Ifges(6,5) * t42 + t49 * t40) * qJD(3);
t5 = -t27 * t63 + (Ifges(6,6) * t42 + t48 * t40) * qJD(3);
t29 = [t40 * t54 + 0.2e1 * t2 * t24 + 0.2e1 * t1 * t25 + 0.2e1 * t4 * t8 + 0.2e1 * t3 * t9 - 0.2e1 * t14 * t18 - 0.2e1 * t20 * t7 + (t4 * t1 - t20 * t14 + t3 * t2) * t89 + 0.2e1 * (-m(5) * t17 + t75) * t16 + (-0.2e1 * t16 * mrSges(5,2) - t39 * t6 - t41 * t5 + (t12 * t39 + (-t13 - t79) * t41) * qJD(5)) * t42 + ((mrSges(4,1) * t93 + mrSges(5,2) * t88 + t78 + t74 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t40) * t40 + (mrSges(4,2) * t93 + mrSges(5,3) * t88 + (0.2e1 * Ifges(4,4) + 0.2e1 * Ifges(5,6) + t47) * t42 + ((2 * Ifges(4,1)) - (2 * Ifges(4,2)) + (2 * Ifges(5,2)) - (2 * Ifges(5,3)) + Ifges(6,3)) * t40) * t42) * qJD(3); (t25 * t69 + t24 * t67 - t7 + m(6) * (t3 * t67 + t4 * t69 - t14)) * t40 + (qJD(3) * t18 + m(6) * (t3 * t65 - t4 * t64 + t15 - t52) - t44) * t42; (0.1e1 - t71) * t40 * t66 * t89; t6 * t85 + t5 * t87 + t20 * t21 - t14 * t26 + qJD(4) * t18 - qJ(4) * t7 + m(6) * (-qJ(4) * t14 + qJD(4) * t20 + t1 * t76 + t2 * t72) + t9 * t72 + t8 * t76 - t52 * mrSges(6,3) + (t55 * qJD(4) - t22 * t86 - t23 * t87) * t42 + (t40 * t47 / 0.2e1 - t74 / 0.2e1 - t78 / 0.2e1 + (t28 * t86 + t39 * t27 / 0.2e1) * t42 + t51 * mrSges(6,3) + (-m(6) * t51 - t39 * t24 + t41 * t25) * t43) * qJD(5) + ((t73 / 0.2e1 + t77 / 0.2e1 - qJ(4) * mrSges(5,1) + Ifges(5,5) - Ifges(4,6)) * t40 + (t40 * mrSges(4,2) - t75 + m(5) * (-t70 - t84)) * t33 + (-pkin(3) * mrSges(5,1) + Ifges(6,5) * t85 + Ifges(6,6) * t87 + t33 * t92 - Ifges(5,4) + Ifges(4,5)) * t42) * qJD(3); t40 * t21 + m(6) * t62 + m(5) * t16 + ((m(6) * qJ(4) - mrSges(4,2) + t91) * t42 + (-t71 * mrSges(6,3) + t43 * t56 + t92) * t40) * qJD(3); 0.2e1 * qJ(4) * t21 + t39 * t22 - t41 * t23 + (-t73 - t77) * qJD(5) + 0.2e1 * ((m(5) + m(6)) * qJ(4) + t91) * qJD(4); m(6) * (-t51 * qJD(5) + t52) + t55 * t66 + t44; (m(5) + t56) * t68; 0; 0; t2 * mrSges(6,1) - t1 * mrSges(6,2) - Ifges(6,5) * t60 + t54; t7; ((-mrSges(6,2) * t43 - Ifges(6,6)) * t41 + (-mrSges(6,1) * t43 - Ifges(6,5)) * t39) * qJD(5); -t26 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t29(1), t29(2), t29(4), t29(7), t29(11); t29(2), t29(3), t29(5), t29(8), t29(12); t29(4), t29(5), t29(6), t29(9), t29(13); t29(7), t29(8), t29(9), t29(10), t29(14); t29(11), t29(12), t29(13), t29(14), t29(15);];
Mq = res;
