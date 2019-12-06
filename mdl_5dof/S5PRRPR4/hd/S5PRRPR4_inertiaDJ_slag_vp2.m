% Calculate time derivative of joint inertia matrix for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:31
% EndTime: 2019-12-05 16:21:34
% DurationCPUTime: 0.73s
% Computational Cost: add. (1111->154), mult. (2795->252), div. (0->0), fcn. (2545->8), ass. (0->69)
t60 = sin(pkin(9));
t61 = cos(pkin(9));
t63 = sin(qJ(3));
t66 = cos(qJ(3));
t49 = t60 * t66 + t61 * t63;
t46 = t49 * qJD(3);
t69 = t60 * t63 - t61 * t66;
t47 = t69 * qJD(3);
t26 = t46 * mrSges(5,1) - t47 * mrSges(5,2);
t62 = sin(qJ(5));
t65 = cos(qJ(5));
t27 = -t49 * t62 - t65 * t69;
t13 = t27 * qJD(5) - t46 * t62 - t47 * t65;
t28 = t49 * t65 - t62 * t69;
t14 = -t28 * qJD(5) - t46 * t65 + t47 * t62;
t4 = -t14 * mrSges(6,1) + t13 * mrSges(6,2);
t84 = -t26 - t4;
t80 = t63 ^ 2 + t66 ^ 2;
t83 = pkin(3) * t60;
t81 = -qJ(4) - pkin(6);
t72 = qJD(3) * t81;
t44 = qJD(4) * t66 + t63 * t72;
t45 = -qJD(4) * t63 + t66 * t72;
t23 = t61 * t44 + t60 * t45;
t54 = t81 * t63;
t55 = t81 * t66;
t31 = t60 * t54 - t61 * t55;
t67 = cos(qJ(2));
t79 = qJD(2) * t67;
t64 = sin(qJ(2));
t78 = qJD(3) * t64;
t77 = qJD(3) * t66;
t76 = 0.2e1 * t63;
t75 = pkin(3) * qJD(3) * t63;
t74 = t64 * t79;
t39 = t49 * t64;
t40 = t69 * t64;
t18 = -t39 * t65 + t40 * t62;
t20 = -t49 * t79 + t69 * t78;
t21 = -t64 * t46 - t69 * t79;
t6 = t18 * qJD(5) + t20 * t62 + t21 * t65;
t19 = -t39 * t62 - t40 * t65;
t7 = -t19 * qJD(5) + t20 * t65 - t21 * t62;
t73 = t7 * mrSges(6,1) - t6 * mrSges(6,2);
t57 = -pkin(3) * t66 - pkin(2);
t22 = -t44 * t60 + t61 * t45;
t30 = t61 * t54 + t55 * t60;
t71 = -t66 * mrSges(4,1) + t63 * mrSges(4,2);
t56 = pkin(3) * t61 + pkin(4);
t42 = t56 * t65 - t62 * t83;
t37 = t42 * qJD(5);
t43 = t56 * t62 + t65 * t83;
t38 = t43 * qJD(5);
t70 = -t38 * mrSges(6,1) - t37 * mrSges(6,2);
t24 = -pkin(7) * t49 + t30;
t25 = -pkin(7) * t69 + t31;
t8 = t24 * t65 - t25 * t62;
t9 = t24 * t62 + t25 * t65;
t16 = pkin(7) * t47 + t22;
t17 = -pkin(7) * t46 + t23;
t2 = t8 * qJD(5) + t16 * t62 + t17 * t65;
t3 = -t9 * qJD(5) + t16 * t65 - t17 * t62;
t68 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t13 + Ifges(6,6) * t14;
t53 = (mrSges(4,1) * t63 + mrSges(4,2) * t66) * qJD(3);
t33 = pkin(4) * t69 + t57;
t32 = pkin(4) * t46 + t75;
t29 = mrSges(5,1) * t69 + mrSges(5,2) * t49;
t15 = -mrSges(6,1) * t27 + mrSges(6,2) * t28;
t1 = [0.2e1 * m(6) * (t18 * t7 + t19 * t6 - t74) + 0.2e1 * m(5) * (-t39 * t20 - t40 * t21 - t74) + 0.2e1 * m(4) * (-0.1e1 + t80) * t74; (-t53 + t84) * t67 + m(6) * (t3 * t18 + t2 * t19 - t32 * t67 + t9 * t6 + t8 * t7) + m(5) * (t30 * t20 + t31 * t21 - t22 * t39 - t23 * t40 - t67 * t75) + (-t13 * t18 + t14 * t19 + t27 * t6 - t28 * t7) * mrSges(6,3) + (-t20 * t49 - t21 * t69 - t39 * t47 + t40 * t46) * mrSges(5,3) + ((-m(4) * pkin(2) + m(5) * t57 + m(6) * t33 - mrSges(3,1) + t15 + t29 + t71) * t64 + (-mrSges(3,2) + (m(4) * pkin(6) + mrSges(4,3)) * t80) * t67) * qJD(2); -0.2e1 * t49 * t47 * Ifges(5,1) + 0.2e1 * t13 * t28 * Ifges(6,1) + 0.2e1 * t69 * Ifges(5,2) * t46 + 0.2e1 * t27 * Ifges(6,2) * t14 - 0.2e1 * pkin(2) * t53 + 0.2e1 * t32 * t15 + 0.2e1 * t57 * t26 + 0.2e1 * t33 * t4 + (-Ifges(4,4) * t63 + pkin(3) * t29) * qJD(3) * t76 + 0.2e1 * m(5) * (t22 * t30 + t23 * t31 + t57 * t75) + 0.2e1 * m(6) * (t2 * t9 + t3 * t8 + t32 * t33) + (0.2e1 * Ifges(4,4) * t66 + (Ifges(4,1) - Ifges(4,2)) * t76) * t77 + 0.2e1 * (t27 * t13 + t14 * t28) * Ifges(6,4) + 0.2e1 * (-t49 * t46 + t47 * t69) * Ifges(5,4) + 0.2e1 * (-t8 * t13 + t9 * t14 + t2 * t27 - t3 * t28) * mrSges(6,3) + 0.2e1 * (-t22 * t49 - t23 * t69 + t30 * t47 - t31 * t46) * mrSges(5,3); t20 * mrSges(5,1) - t21 * mrSges(5,2) + (t63 * t78 - t66 * t79) * mrSges(4,2) + (-t63 * t79 - t64 * t77) * mrSges(4,1) + m(6) * (-t18 * t38 + t19 * t37 + t42 * t7 + t43 * t6) + m(5) * (t20 * t61 + t21 * t60) * pkin(3) + t73; m(6) * (t2 * t43 + t3 * t42 + t37 * t9 - t38 * t8) + t22 * mrSges(5,1) - t23 * mrSges(5,2) - Ifges(5,5) * t47 - Ifges(5,6) * t46 + (Ifges(4,5) * t66 - Ifges(4,6) * t63 + t71 * pkin(6)) * qJD(3) + (m(5) * (t22 * t61 + t23 * t60) + (-t46 * t60 + t47 * t61) * mrSges(5,3)) * pkin(3) + (-t13 * t42 + t14 * t43 + t27 * t37 + t28 * t38) * mrSges(6,3) + t68; 0.2e1 * m(6) * (t37 * t43 - t38 * t42) + 0.2e1 * t70; (m(5) + m(6)) * t64 * qJD(2); m(5) * t75 + m(6) * t32 - t84; 0; 0; t73; t68; t70; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
