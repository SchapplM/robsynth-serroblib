% Calculate time derivative of joint inertia matrix for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:47
% EndTime: 2019-12-05 15:44:48
% DurationCPUTime: 0.27s
% Computational Cost: add. (413->70), mult. (1053->116), div. (0->0), fcn. (924->8), ass. (0->43)
t34 = sin(pkin(9));
t35 = cos(pkin(9));
t38 = sin(qJ(2));
t41 = cos(qJ(2));
t24 = t34 * t41 + t35 * t38;
t18 = t24 * qJD(2);
t23 = -t34 * t38 + t35 * t41;
t19 = t23 * qJD(2);
t37 = sin(qJ(4));
t40 = cos(qJ(4));
t44 = t23 * t40 - t24 * t37;
t5 = qJD(4) * t44 - t37 * t18 + t40 * t19;
t36 = sin(qJ(5));
t39 = cos(qJ(5));
t51 = t36 ^ 2 + t39 ^ 2;
t48 = t51 * t5;
t26 = -mrSges(6,1) * t39 + mrSges(6,2) * t36;
t59 = -mrSges(5,1) + t26;
t61 = Ifges(6,1) - Ifges(6,2);
t30 = pkin(2) * t35 + pkin(3);
t56 = pkin(2) * t34;
t15 = t30 * t40 - t37 * t56;
t10 = t15 * qJD(4);
t60 = t10 * mrSges(5,2);
t47 = t51 * t10;
t16 = t37 * t30 + t40 * t56;
t9 = t23 * t37 + t24 * t40;
t6 = qJD(4) * t9 + t40 * t18 + t37 * t19;
t58 = t44 * t6;
t11 = t16 * qJD(4);
t55 = t11 * t44;
t52 = Ifges(6,6) * t36;
t50 = qJD(5) * t36;
t49 = qJD(5) * t39;
t14 = pkin(7) + t16;
t46 = t51 * t14;
t45 = mrSges(6,1) * t36 + mrSges(6,2) * t39;
t25 = t45 * qJD(5);
t43 = -t5 * mrSges(5,2) + mrSges(6,3) * t48 - t44 * t25 + t59 * t6;
t42 = (-0.2e1 * Ifges(6,4) * t36 + t39 * t61) * t50 + (0.2e1 * Ifges(6,4) * t39 + t36 * t61) * t49;
t31 = Ifges(6,5) * t49;
t13 = -pkin(4) - t15;
t1 = [0.2e1 * m(6) * (t48 * t9 - t58) + 0.2e1 * m(5) * (t5 * t9 - t58) + 0.2e1 * m(4) * (-t18 * t23 + t19 * t24); -t18 * mrSges(4,1) - t19 * mrSges(4,2) + (-mrSges(3,1) * t38 - mrSges(3,2) * t41) * qJD(2) + m(6) * (t13 * t6 + t46 * t5 + t47 * t9 - t55) + m(5) * (t10 * t9 - t15 * t6 + t16 * t5 - t55) + m(4) * (-t18 * t35 + t19 * t34) * pkin(2) + t43; 0.2e1 * t13 * t25 - 0.2e1 * t60 + 0.2e1 * m(6) * (t10 * t46 + t13 * t11) + 0.2e1 * m(5) * (t10 * t16 - t11 * t15) + t42 + 0.2e1 * t59 * t11 + 0.2e1 * mrSges(6,3) * t47; 0; 0; 0; m(6) * (-pkin(4) * t6 + pkin(7) * t48) + t43; -t60 + (-pkin(4) + t13) * t25 + t42 + (m(6) * pkin(7) + mrSges(6,3)) * t47 + (-m(6) * pkin(4) + t59) * t11; 0; -0.2e1 * pkin(4) * t25 + t42; (-t39 * t5 + t50 * t9) * mrSges(6,2) + (-t36 * t5 - t49 * t9) * mrSges(6,1); t31 - t45 * t10 + (t14 * t26 - t52) * qJD(5); -t25; t31 + (pkin(7) * t26 - t52) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
