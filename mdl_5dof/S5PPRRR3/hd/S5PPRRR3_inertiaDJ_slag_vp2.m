% Calculate time derivative of joint inertia matrix for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:23
% EndTime: 2019-12-05 15:16:26
% DurationCPUTime: 0.72s
% Computational Cost: add. (631->129), mult. (1853->225), div. (0->0), fcn. (1623->8), ass. (0->67)
t44 = sin(qJ(4));
t47 = cos(qJ(4));
t29 = (mrSges(5,1) * t44 + mrSges(5,2) * t47) * qJD(4);
t43 = sin(qJ(5));
t46 = cos(qJ(5));
t52 = t43 * t44 - t46 * t47;
t74 = qJD(4) + qJD(5);
t15 = t74 * t52;
t28 = t43 * t47 + t46 * t44;
t16 = t74 * t28;
t4 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t76 = t29 + t4;
t75 = m(5) * pkin(3);
t45 = sin(qJ(3));
t23 = t52 * t45;
t38 = t44 ^ 2;
t68 = t47 ^ 2 + t38;
t73 = 2 * m(6);
t72 = m(6) * pkin(4);
t71 = -pkin(7) - pkin(6);
t41 = sin(pkin(9));
t48 = cos(qJ(3));
t69 = t41 * t48;
t67 = qJD(3) * t41;
t66 = qJD(3) * t48;
t65 = qJD(4) * t44;
t64 = qJD(4) * t47;
t17 = mrSges(6,1) * t52 + t28 * mrSges(6,2);
t33 = -t47 * mrSges(5,1) + t44 * mrSges(5,2);
t63 = -mrSges(4,1) + t17 + t33;
t62 = pkin(4) * t65;
t61 = t45 * t67;
t60 = t44 * t66;
t59 = t45 * t66;
t58 = t47 * t66;
t57 = t45 * t65;
t42 = cos(pkin(9));
t24 = -t42 * t47 - t44 * t69;
t25 = -t42 * t44 + t47 * t69;
t11 = t46 * t24 - t43 * t25;
t18 = -t25 * qJD(4) + t44 * t61;
t19 = t24 * qJD(4) - t47 * t61;
t2 = t11 * qJD(5) + t43 * t18 + t46 * t19;
t12 = t43 * t24 + t46 * t25;
t3 = -t12 * qJD(5) + t46 * t18 - t43 * t19;
t56 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t6 = -t16 * t45 - t52 * t66;
t7 = t74 * t23 - t28 * t66;
t55 = t7 * mrSges(6,1) - t6 * mrSges(6,2);
t54 = qJD(4) * t71;
t53 = (t45 ^ 2 - t48 ^ 2) * t67;
t34 = t71 * t44;
t35 = t71 * t47;
t20 = t46 * t34 + t43 * t35;
t21 = t43 * t34 - t46 * t35;
t30 = t44 * t54;
t31 = t47 * t54;
t10 = -t21 * qJD(5) - t43 * t30 + t46 * t31;
t9 = t20 * qJD(5) + t46 * t30 + t43 * t31;
t51 = t10 * mrSges(6,1) - t9 * mrSges(6,2) - Ifges(6,5) * t15 - Ifges(6,6) * t16;
t50 = -t18 * t44 + t19 * t47 + (-t24 * t47 - t25 * t44) * qJD(4);
t49 = m(5) * t50;
t37 = -t47 * pkin(4) - pkin(3);
t32 = t41 ^ 2 * t59;
t26 = (-mrSges(6,1) * t43 - mrSges(6,2) * t46) * qJD(5) * pkin(4);
t22 = t28 * t45;
t1 = [0.2e1 * m(6) * (t11 * t3 + t12 * t2 + t32) + 0.2e1 * m(5) * (t24 * t18 + t25 * t19 + t32); m(6) * (t7 * t11 + t6 * t12 - t23 * t2 - t22 * t3 + t53) + m(5) * (-t24 * t60 + t25 * t58 + t53) + t45 * t49; (-t22 * t7 - t23 * t6) * t73 + 0.4e1 * (m(5) * (-0.1e1 + t68) / 0.2e1 - m(6) / 0.2e1) * t59; m(6) * (t10 * t11 + t9 * t12 + t21 * t2 + t20 * t3) + pkin(6) * t49 + (t11 * t15 - t12 * t16 - t2 * t52 - t3 * t28) * mrSges(6,3) + t50 * mrSges(5,3) + (t76 * t45 + (mrSges(4,2) * t45 + t63 * t48) * qJD(3) + m(6) * (pkin(4) * t57 + t37 * t66) - t66 * t75) * t41; m(6) * (-t10 * t22 + t20 * t7 + t21 * t6 - t9 * t23) + (-t22 * t15 + t23 * t16 - t7 * t28 - t52 * t6) * mrSges(6,3) + (m(6) * t37 + t63 - t75) * t45 * qJD(3) + (-m(6) * t62 + (-mrSges(4,2) + (m(5) * pkin(6) + mrSges(5,3)) * t68) * qJD(3) - t76) * t48; (t20 * t10 + t21 * t9 + t37 * t62) * t73 - 0.2e1 * t15 * t28 * Ifges(6,1) + 0.2e1 * t16 * Ifges(6,2) * t52 + 0.2e1 * t17 * t62 + 0.2e1 * t37 * t4 - 0.2e1 * pkin(3) * t29 + 0.2e1 * (-t38 * qJD(4) + t47 * t64) * Ifges(5,4) + 0.2e1 * (t15 * t52 - t28 * t16) * Ifges(6,4) + 0.2e1 * (-t10 * t28 + t20 * t15 - t21 * t16 - t52 * t9) * mrSges(6,3) + 0.2e1 * (Ifges(5,1) - Ifges(5,2)) * t44 * t64; t18 * mrSges(5,1) - t19 * mrSges(5,2) + (t2 * t43 + t3 * t46 + (-t11 * t43 + t12 * t46) * qJD(5)) * t72 + t56; (t57 - t58) * mrSges(5,2) + (-t45 * t64 - t60) * mrSges(5,1) + (t43 * t6 + t46 * t7 + (t22 * t43 - t23 * t46) * qJD(5)) * t72 + t55; (Ifges(5,5) * t47 - Ifges(5,6) * t44 + t33 * pkin(6)) * qJD(4) + (m(6) * (t10 * t46 + t43 * t9 + (-t20 * t43 + t21 * t46) * qJD(5)) + (t46 * t15 - t43 * t16 + (t28 * t43 - t46 * t52) * qJD(5)) * mrSges(6,3)) * pkin(4) + t51; 0.2e1 * t26; t56; t55; t51; t26; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
