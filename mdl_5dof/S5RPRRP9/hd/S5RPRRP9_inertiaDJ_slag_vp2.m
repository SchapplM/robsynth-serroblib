% Calculate time derivative of joint inertia matrix for
% S5RPRRP9
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP9_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP9_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:37
% EndTime: 2019-12-31 18:48:39
% DurationCPUTime: 0.70s
% Computational Cost: add. (1402->124), mult. (3132->188), div. (0->0), fcn. (2955->6), ass. (0->59)
t81 = -mrSges(5,1) - mrSges(6,1);
t80 = (mrSges(6,2) + mrSges(5,3));
t78 = 2 * t80;
t50 = sin(pkin(8));
t69 = pkin(6) + qJ(2);
t41 = t69 * t50;
t51 = cos(pkin(8));
t42 = t69 * t51;
t53 = sin(qJ(3));
t54 = cos(qJ(3));
t28 = -t53 * t41 + t54 * t42;
t71 = t51 * t54;
t39 = -t53 * t50 + t71;
t65 = -pkin(2) * t51 - pkin(1);
t29 = -t39 * pkin(3) + t65;
t77 = 0.2e1 * t29;
t76 = m(5) * pkin(3);
t52 = sin(qJ(4));
t75 = pkin(3) * t52;
t40 = t50 * t54 + t53 * t51;
t34 = t40 * qJD(3);
t74 = t34 * pkin(3);
t73 = cos(qJ(4));
t33 = t39 * qJD(3);
t59 = t73 * t39 - t52 * t40;
t16 = t59 * qJD(4) + t73 * t33 - t52 * t34;
t72 = t16 * mrSges(6,2);
t37 = t54 * t41;
t68 = qJD(4) * t75;
t20 = -qJD(3) * t37 + qJD(2) * t71 + (-qJD(2) * t50 - qJD(3) * t42) * t53;
t19 = -t34 * pkin(7) + t20;
t21 = -t40 * qJD(2) - t28 * qJD(3);
t56 = -t33 * pkin(7) + t21;
t27 = -t42 * t53 - t37;
t23 = -t40 * pkin(7) + t27;
t24 = pkin(7) * t39 + t28;
t60 = t73 * t23 - t52 * t24;
t4 = t60 * qJD(4) + t73 * t19 + t52 * t56;
t9 = t52 * t23 + t73 * t24;
t5 = t9 * qJD(4) + t52 * t19 - t73 * t56;
t67 = t9 * t4 - t5 * t60;
t66 = t73 * pkin(3);
t26 = t52 * t39 + t73 * t40;
t63 = t26 * t68;
t62 = -2 * Ifges(5,4) + 2 * Ifges(6,5);
t61 = qJD(4) * t66;
t17 = t26 * qJD(4) + t52 * t33 + t73 * t34;
t57 = t81 * t5 + (-mrSges(5,2) + mrSges(6,3)) * t4 + (-Ifges(5,6) + Ifges(6,6)) * t17 + (Ifges(6,4) + Ifges(5,5)) * t16;
t43 = t61 + qJD(5);
t55 = -mrSges(5,2) * t61 + t43 * mrSges(6,3) + t81 * t68;
t49 = qJD(5) * mrSges(6,3);
t46 = -t66 - pkin(4);
t45 = qJ(5) + t75;
t31 = t33 * mrSges(4,2);
t11 = t17 * mrSges(6,1);
t10 = t16 * mrSges(5,2);
t7 = -pkin(4) * t59 - t26 * qJ(5) + t29;
t6 = pkin(4) * t17 - qJ(5) * t16 - qJD(5) * t26 + t74;
t1 = [0.2e1 * t7 * t11 + 0.2e1 * t65 * (t34 * mrSges(4,1) + t31) - 0.2e1 * t39 * Ifges(4,2) * t34 + 0.2e1 * t40 * t33 * Ifges(4,1) + 0.2e1 * t6 * (-mrSges(6,1) * t59 - t26 * mrSges(6,3)) + t10 * t77 + 0.2e1 * (-mrSges(5,1) * t59 + mrSges(5,2) * t26) * t74 + 0.2e1 * (t39 * t33 - t40 * t34) * Ifges(4,4) + 0.2e1 * (t20 * t39 - t21 * t40 - t27 * t33 - t28 * t34) * mrSges(4,3) + (mrSges(5,1) * t77 + t26 * t62 - t9 * t78 - 0.2e1 * (Ifges(5,2) + Ifges(6,3)) * t59) * t17 + (-0.2e1 * mrSges(6,3) * t7 - t59 * t62 - t60 * t78 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t26) * t16 + 0.2e1 * m(4) * (t20 * t28 + t21 * t27) + 0.2e1 * m(5) * (t29 * t74 + t67) + 0.2e1 * m(6) * (t6 * t7 + t67) + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t50 ^ 2 + t51 ^ 2) * qJD(2) + 0.2e1 * t80 * (t5 * t26 + t4 * t59); t10 + t17 * mrSges(5,1) + m(6) * t6 - t16 * mrSges(6,3) + t11 + t31 - (-mrSges(4,1) - t76) * t34; 0; m(6) * (t45 * t4 + t43 * t9 + t46 * t5 - t60 * t68) - Ifges(4,6) * t34 + Ifges(4,5) * t33 - t20 * mrSges(4,2) + t21 * mrSges(4,1) + t46 * t72 + t57 + (-t73 * t5 + t4 * t52 + (-t52 * t60 + t73 * t9) * qJD(4)) * t76 + (-t45 * t17 + t43 * t59 + t63) * mrSges(6,2) + (-t16 * t66 - t17 * t75 + t59 * t61 + t63) * mrSges(5,3); 0; 0.2e1 * m(6) * (t45 * t43 + t46 * t68) + 0.2e1 * t55; m(6) * (-pkin(4) * t5 + qJ(5) * t4 + qJD(5) * t9) + (-pkin(4) * t16 - qJ(5) * t17 + qJD(5) * t59) * mrSges(6,2) + t57; 0; t49 + m(6) * (-pkin(4) * t68 + qJ(5) * t43 + qJD(5) * t45) + t55; 0.2e1 * m(6) * qJ(5) * qJD(5) + 0.2e1 * t49; m(6) * t5 + t72; 0; m(6) * t68; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
