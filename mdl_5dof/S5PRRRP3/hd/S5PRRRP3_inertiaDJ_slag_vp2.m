% Calculate time derivative of joint inertia matrix for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:32
% EndTime: 2019-12-05 16:43:34
% DurationCPUTime: 0.58s
% Computational Cost: add. (666->110), mult. (1615->172), div. (0->0), fcn. (1274->4), ass. (0->50)
t43 = cos(qJ(3));
t39 = -t43 * pkin(3) - pkin(2);
t59 = 0.2e1 * t39;
t67 = 2 * mrSges(5,3);
t66 = 2 * mrSges(6,3);
t42 = cos(qJ(4));
t65 = t42 * pkin(3);
t40 = sin(qJ(4));
t41 = sin(qJ(3));
t30 = t40 * t43 + t42 * t41;
t60 = qJD(3) + qJD(4);
t20 = t60 * t30;
t29 = -t40 * t41 + t42 * t43;
t64 = t29 * t20;
t19 = t60 * t29;
t63 = t30 * t19;
t62 = -mrSges(5,1) - mrSges(6,1);
t58 = -pkin(7) - pkin(6);
t36 = t58 * t41;
t37 = t58 * t43;
t22 = t40 * t36 - t42 * t37;
t61 = pkin(3) * qJD(4);
t57 = t20 * pkin(4);
t53 = qJD(4) * t42;
t56 = (t19 * t40 + t30 * t53) * pkin(3);
t55 = -t20 * mrSges(6,1) - t19 * mrSges(6,2);
t54 = qJD(4) * t40;
t52 = 0.2e1 * t43;
t51 = t29 * t54;
t50 = qJD(3) * t58;
t49 = (-mrSges(5,2) - mrSges(6,2)) * t42;
t21 = t42 * t36 + t40 * t37;
t48 = 2 * Ifges(5,4) + 2 * Ifges(6,4);
t34 = t41 * t50;
t35 = t43 * t50;
t7 = -t22 * qJD(4) - t40 * t34 + t42 * t35;
t3 = -t19 * qJ(5) - t30 * qJD(5) + t7;
t47 = m(6) * t3 - t19 * mrSges(6,3);
t14 = t20 * mrSges(5,1);
t46 = -t19 * mrSges(5,2) - t14 + t55;
t6 = t42 * t34 + t40 * t35 + t36 * t53 + t37 * t54;
t45 = -t40 * t20 + (t29 * t42 + t30 * t40) * qJD(4);
t2 = -t20 * qJ(5) + t29 * qJD(5) + t6;
t44 = t7 * mrSges(5,1) + t3 * mrSges(6,1) - t6 * mrSges(5,2) - t2 * mrSges(6,2) - (Ifges(5,6) + Ifges(6,6)) * t20 + (Ifges(5,5) + Ifges(6,5)) * t19;
t38 = pkin(4) + t65;
t23 = -t29 * pkin(4) + t39;
t10 = qJD(3) * t41 * pkin(3) + t57;
t9 = t29 * qJ(5) + t22;
t8 = -t30 * qJ(5) + t21;
t1 = [0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t63 - t64); m(5) * (t22 * t19 - t21 * t20 + t7 * t29 + t6 * t30) + m(6) * (t9 * t19 + t2 * t30 - t8 * t20 + t3 * t29); t14 * t59 + 0.2e1 * t10 * (-t29 * mrSges(6,1) + t30 * mrSges(6,2)) - 0.2e1 * t23 * t55 + 0.2e1 * m(5) * (t21 * t7 + t22 * t6) + 0.2e1 * m(6) * (t23 * t10 + t9 * t2 + t8 * t3) + (mrSges(5,2) * t59 - 0.2e1 * t21 * mrSges(5,3) - 0.2e1 * t8 * mrSges(6,3) + t29 * t48) * t19 - (t22 * t67 + t30 * t48 + t9 * t66) * t20 + ((-pkin(2) * mrSges(4,2) + Ifges(4,4) * t43) * t52 + (m(5) * pkin(3) * t59 + 0.2e1 * pkin(3) * (-t29 * mrSges(5,1) + t30 * mrSges(5,2)) - 0.2e1 * pkin(2) * mrSges(4,1) - 0.2e1 * Ifges(4,4) * t41 + (Ifges(4,1) - Ifges(4,2)) * t52) * t41) * qJD(3) + (t2 * t29 - t3 * t30) * t66 + (t6 * t29 - t7 * t30) * t67 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t63 - 0.2e1 * (Ifges(5,2) + Ifges(6,2)) * t64; (-t41 * mrSges(4,1) - t43 * mrSges(4,2)) * qJD(3) + m(5) * ((-t20 * t42 - t51) * pkin(3) + t56) + m(6) * (-pkin(3) * t51 - t38 * t20 + t56) + t46; t47 * t38 + (Ifges(4,5) * t43 - Ifges(4,6) * t41 + (-mrSges(4,1) * t43 + mrSges(4,2) * t41) * pkin(6)) * qJD(3) + (m(6) * (t2 * t40 + t9 * t53 - t8 * t54) + m(5) * (-t21 * t54 + t22 * t53 + t40 * t6 + t42 * t7) + t45 * mrSges(6,3) + (-t42 * t19 + t45) * mrSges(5,3)) * pkin(3) + t44; 0.2e1 * (t49 + ((-t38 + t65) * m(6) + t62) * t40) * t61; -m(6) * t57 + t46; t47 * pkin(4) + t44; (t49 + (-m(6) * pkin(4) + t62) * t40) * t61; 0; 0; m(6) * t10 - t55; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
