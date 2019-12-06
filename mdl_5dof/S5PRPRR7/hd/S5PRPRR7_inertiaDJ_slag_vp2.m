% Calculate time derivative of joint inertia matrix for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR7_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR7_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:45
% EndTime: 2019-12-05 15:59:46
% DurationCPUTime: 0.40s
% Computational Cost: add. (545->101), mult. (1292->169), div. (0->0), fcn. (1026->6), ass. (0->56)
t39 = sin(qJ(4));
t42 = cos(qJ(4));
t26 = t39 * mrSges(5,1) + t42 * mrSges(5,2);
t71 = mrSges(4,3) + t26;
t70 = qJD(4) + qJD(5);
t38 = sin(qJ(5));
t57 = qJD(4) * t39;
t41 = cos(qJ(5));
t62 = t41 * t42;
t63 = t38 * t39;
t11 = -qJD(5) * t63 - t38 * t57 + t70 * t62;
t48 = t38 * t42 + t41 * t39;
t12 = t70 * t48;
t47 = -t62 + t63;
t69 = (-t38 * t47 - t41 * t48) * qJD(5) - t11 * t38 + t12 * t41;
t43 = cos(qJ(2));
t68 = t70 * t43;
t37 = t42 ^ 2;
t67 = m(6) * pkin(4);
t44 = -pkin(2) - pkin(6);
t66 = pkin(7) - t44;
t65 = t48 * t11;
t64 = t47 * t12;
t40 = sin(qJ(2));
t58 = qJD(2) * t43;
t61 = qJ(3) * t58 + t40 * qJD(3);
t60 = t39 ^ 2 + t37;
t59 = qJD(2) * t40;
t56 = qJD(4) * t43;
t55 = 2 * mrSges(6,3);
t29 = t40 * t58;
t6 = t47 * t68 + t48 * t59;
t7 = -t47 * t59 + t48 * t68;
t54 = t7 * mrSges(6,1) - t6 * mrSges(6,2);
t25 = t66 * t42;
t53 = -t12 * mrSges(6,1) - t11 * mrSges(6,2);
t24 = t66 * t39;
t14 = t38 * t24 - t41 * t25;
t19 = t66 * t57;
t20 = qJD(4) * t25;
t2 = qJD(5) * t14 + t38 * t19 - t41 * t20;
t15 = -t41 * t24 - t38 * t25;
t3 = -qJD(5) * t15 + t41 * t19 + t38 * t20;
t52 = t3 * mrSges(6,1) - t2 * mrSges(6,2) - Ifges(6,5) * t12 - Ifges(6,6) * t11;
t51 = -t64 - t65;
t46 = t11 * t15 - t12 * t14 + t2 * t48 - t3 * t47;
t16 = t47 * t43;
t17 = t48 * t43;
t45 = -t11 * t17 - t12 * t16 - t47 * t7 + t48 * t6;
t31 = t39 * pkin(4) + qJ(3);
t27 = qJD(4) * t42 * pkin(4) + qJD(3);
t23 = (mrSges(5,1) * t42 - mrSges(5,2) * t39) * qJD(4);
t18 = (-mrSges(6,1) * t38 - mrSges(6,2) * t41) * qJD(5) * pkin(4);
t13 = mrSges(6,1) * t48 - mrSges(6,2) * t47;
t4 = t11 * mrSges(6,1) - t12 * mrSges(6,2);
t1 = [0.2e1 * m(5) * (-t60 + 0.1e1) * t29 + 0.2e1 * m(6) * (t16 * t7 - t17 * t6 + t29); (t23 + t4) * t40 + ((-mrSges(3,2) + t13 + t71) * t43 + (-mrSges(5,3) * t60 - mrSges(3,1) + mrSges(4,2)) * t40) * qJD(2) + m(4) * (-pkin(2) * t59 + t61) + m(5) * (t44 * t59 * t60 + t61) + m(6) * (t14 * t7 + t15 * t6 + t3 * t16 - t2 * t17 + t27 * t40 + t31 * t58) - t45 * mrSges(6,3); 0.2e1 * Ifges(6,1) * t64 + 0.2e1 * Ifges(6,2) * t65 + 0.2e1 * m(6) * (t14 * t3 + t15 * t2 + t31 * t27) + 0.2e1 * t27 * t13 + 0.2e1 * t31 * t4 - t46 * t55 + 0.2e1 * (t11 * t47 + t12 * t48) * Ifges(6,4) + 0.2e1 * (-Ifges(5,1) + Ifges(5,2)) * t42 * t57 + 0.2e1 * (t23 + (m(4) + m(5)) * qJD(3)) * qJ(3) + 0.2e1 * t71 * qJD(3) + 0.2e1 * (-t37 * qJD(4) + t39 * t57) * Ifges(5,4); m(6) * t45 + (m(5) * t60 + m(4)) * t59; m(6) * t46 + t51 * t55; -0.2e1 * m(6) * t51; (-t39 * t59 + t42 * t56) * mrSges(5,2) + (t39 * t56 + t42 * t59) * mrSges(5,1) + (t38 * t6 + t41 * t7 + (-t16 * t38 - t17 * t41) * qJD(5)) * t67 + t54; ((-mrSges(5,2) * t44 - Ifges(5,6)) * t42 + (-mrSges(5,1) * t44 - Ifges(5,5)) * t39) * qJD(4) + (m(6) * (t2 * t38 + t3 * t41 + (-t14 * t38 + t15 * t41) * qJD(5)) + t69 * mrSges(6,3)) * pkin(4) + t52; -t26 * qJD(4) - t69 * t67 + t53; 0.2e1 * t18; t54; t52; t53; t18; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
