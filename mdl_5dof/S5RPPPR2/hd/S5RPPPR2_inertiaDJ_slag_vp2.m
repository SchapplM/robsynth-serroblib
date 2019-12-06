% Calculate time derivative of joint inertia matrix for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:30:59
% EndTime: 2019-12-05 17:31:02
% DurationCPUTime: 0.54s
% Computational Cost: add. (774->120), mult. (1833->223), div. (0->0), fcn. (1773->8), ass. (0->61)
t70 = 2 * m(6);
t44 = sin(pkin(7));
t40 = t44 ^ 2;
t42 = sin(pkin(9));
t45 = cos(pkin(9));
t47 = cos(pkin(7));
t46 = cos(pkin(8));
t63 = t44 * t46;
t28 = -t42 * t47 + t45 * t63;
t48 = sin(qJ(5));
t43 = sin(pkin(8));
t49 = cos(qJ(5));
t64 = t43 * t49;
t18 = -t28 * t48 + t44 * t64;
t69 = 0.2e1 * t18;
t65 = t43 * t48;
t19 = t28 * t49 + t44 * t65;
t68 = 0.2e1 * t19;
t67 = t42 * t43;
t66 = t43 * t44;
t55 = qJD(3) * t44;
t56 = qJD(2) * t47;
t30 = -t43 * t55 + t46 * t56;
t22 = -qJD(4) * t47 + t30;
t50 = (-qJD(4) * t46 + qJD(2)) * t44;
t11 = t42 * t22 - t45 * t50;
t62 = t45 * t11;
t29 = t43 * t56 + t46 * t55;
t61 = t46 * t29;
t15 = t18 * qJD(5);
t16 = t19 * qJD(5);
t60 = Ifges(6,5) * t15 - Ifges(6,6) * t16;
t34 = -pkin(2) * t47 - qJ(3) * t44 - pkin(1);
t57 = qJ(2) * t47;
t58 = t43 * t34 + t46 * t57;
t20 = -qJ(4) * t47 + t58;
t23 = (pkin(3) * t43 - qJ(4) * t46 + qJ(2)) * t44;
t59 = t45 * t20 + t42 * t23;
t54 = qJ(2) * qJD(2);
t53 = t34 * t46 - t43 * t57;
t52 = t47 * pkin(3) - t53;
t27 = t42 * t63 + t45 * t47;
t6 = pkin(4) * t27 - pkin(6) * t28 + t52;
t8 = pkin(6) * t66 + t59;
t3 = -t48 * t8 + t49 * t6;
t4 = t48 * t6 + t49 * t8;
t51 = -t20 * t42 + t23 * t45;
t32 = t45 * t64 - t46 * t48;
t31 = -t45 * t65 - t46 * t49;
t41 = t47 ^ 2;
t38 = t40 * t54;
t25 = t32 * qJD(5);
t24 = t31 * qJD(5);
t12 = t45 * t22 + t42 * t50;
t10 = mrSges(6,1) * t27 - mrSges(6,3) * t19;
t9 = -mrSges(6,2) * t27 + mrSges(6,3) * t18;
t7 = -pkin(4) * t66 - t51;
t5 = mrSges(6,1) * t16 + mrSges(6,2) * t15;
t2 = -qJD(5) * t4 - t48 * t12 + t49 * t29;
t1 = qJD(5) * t3 + t49 * t12 + t48 * t29;
t13 = [0.2e1 * t11 * (-t18 * mrSges(6,1) + t19 * mrSges(6,2)) + 0.2e1 * t7 * t5 + 0.2e1 * t1 * t9 + 0.2e1 * t2 * t10 + t27 * t60 + 0.2e1 * t12 * (-mrSges(5,2) * t66 - t27 * mrSges(5,3)) - 0.2e1 * t11 * (mrSges(5,1) * t66 - t28 * mrSges(5,3)) - 0.2e1 * t29 * (-t47 * mrSges(4,1) - mrSges(4,3) * t63) + 0.2e1 * t29 * (t27 * mrSges(5,1) + t28 * mrSges(5,2)) + 0.2e1 * t30 * (t47 * mrSges(4,2) - mrSges(4,3) * t66) - (0.2e1 * mrSges(6,3) * t4 + Ifges(6,4) * t68 + Ifges(6,2) * t69 + Ifges(6,6) * t27) * t16 + (-0.2e1 * mrSges(6,3) * t3 + Ifges(6,1) * t68 + Ifges(6,4) * t69 + Ifges(6,5) * t27) * t15 + 0.2e1 * (t40 * (mrSges(4,1) * t43 + mrSges(4,2) * t46) + (t41 + t40) * mrSges(3,3)) * qJD(2) + (t1 * t4 + t11 * t7 + t2 * t3) * t70 + 0.2e1 * m(5) * (-t11 * t51 + t12 * t59 + t29 * t52) + 0.2e1 * m(4) * (-t29 * t53 + t30 * t58 + t38) + 0.2e1 * m(3) * (t41 * t54 + t38); t5 * t67 - t25 * t10 + t24 * t9 + (-t15 * t31 - t16 * t32) * mrSges(6,3) + m(6) * (t1 * t32 + t11 * t67 + t2 * t31 + t24 * t4 - t25 * t3) + m(5) * (-t61 + (t11 * t42 + t12 * t45) * t43) + m(4) * (t30 * t43 - t61); (t24 * t32 - t25 * t31) * t70; m(4) * t44 * qJD(2) - t45 * t5 + ((-t10 * t49 - t48 * t9) * qJD(5) + (t15 * t48 - t16 * t49) * mrSges(6,3)) * t42 + m(6) * (-t62 + (t1 * t49 - t2 * t48 + (-t3 * t49 - t4 * t48) * qJD(5)) * t42) + m(5) * (t12 * t42 - t62); m(6) * (t24 * t49 + t25 * t48 + (-t31 * t49 - t32 * t48) * qJD(5)) * t42; 0; (-t48 * t10 + t49 * t9) * qJD(5) + (-t15 * t49 - t16 * t48) * mrSges(6,3) + m(6) * (t48 * t1 + t49 * t2 + (-t3 * t48 + t4 * t49) * qJD(5)) + m(5) * t29; m(6) * (t48 * t24 - t49 * t25 + (-t31 * t48 + t32 * t49) * qJD(5)); 0; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 + t60; -mrSges(6,1) * t25 - mrSges(6,2) * t24; (-mrSges(6,1) * t49 + mrSges(6,2) * t48) * t42 * qJD(5); (-mrSges(6,1) * t48 - mrSges(6,2) * t49) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11); t13(2), t13(3), t13(5), t13(8), t13(12); t13(4), t13(5), t13(6), t13(9), t13(13); t13(7), t13(8), t13(9), t13(10), t13(14); t13(11), t13(12), t13(13), t13(14), t13(15);];
Mq = res;
