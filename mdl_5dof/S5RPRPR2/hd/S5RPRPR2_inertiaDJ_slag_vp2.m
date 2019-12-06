% Calculate time derivative of joint inertia matrix for
% S5RPRPR2
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:21
% EndTime: 2019-12-05 17:49:22
% DurationCPUTime: 0.35s
% Computational Cost: add. (665->93), mult. (1450->138), div. (0->0), fcn. (1138->8), ass. (0->53)
t40 = sin(pkin(9));
t42 = cos(pkin(9));
t43 = sin(qJ(5));
t45 = cos(qJ(5));
t27 = t45 * t40 + t43 * t42;
t25 = t27 * qJD(5);
t26 = -t43 * t40 + t45 * t42;
t71 = t26 * t25;
t24 = t26 * qJD(5);
t70 = t27 * t24;
t35 = cos(pkin(8)) * pkin(1) + pkin(2);
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t60 = pkin(1) * sin(pkin(8));
t50 = t46 * t35 - t44 * t60;
t69 = 2 * mrSges(5,3);
t68 = 2 * mrSges(6,3);
t19 = t50 * qJD(3);
t18 = qJD(4) + t19;
t54 = t40 ^ 2 + t42 ^ 2;
t51 = t54 * t18;
t67 = -t42 * mrSges(5,1) - t26 * mrSges(6,1) + t40 * mrSges(5,2) + t27 * mrSges(6,2);
t55 = t44 * t35 + t46 * t60;
t48 = t54 * qJD(4);
t66 = 2 * m(5);
t65 = 2 * m(6);
t11 = t25 * mrSges(6,1) + t24 * mrSges(6,2);
t62 = 0.2e1 * t11;
t20 = t55 * qJD(3);
t61 = 0.2e1 * t20;
t59 = t42 * pkin(4);
t58 = t19 * mrSges(4,2);
t56 = Ifges(6,5) * t24 - Ifges(6,6) * t25;
t52 = 0.2e1 * Ifges(6,1) * t70 - 0.2e1 * Ifges(6,2) * t71 + 0.2e1 * (t24 * t26 - t25 * t27) * Ifges(6,4);
t49 = t54 * qJ(4);
t47 = -pkin(3) - t50;
t21 = qJ(4) + t55;
t15 = (-pkin(7) - t21) * t40;
t37 = t42 * pkin(7);
t16 = t42 * t21 + t37;
t3 = t45 * t15 - t43 * t16;
t4 = t43 * t15 + t45 * t16;
t29 = (-pkin(7) - qJ(4)) * t40;
t31 = t42 * qJ(4) + t37;
t13 = t45 * t29 - t43 * t31;
t14 = t43 * t29 + t45 * t31;
t36 = -pkin(3) - t59;
t17 = t47 - t59;
t10 = -t27 * qJD(4) - t14 * qJD(5);
t9 = t26 * qJD(4) + t13 * qJD(5);
t2 = -t4 * qJD(5) - t27 * t18;
t1 = t3 * qJD(5) + t26 * t18;
t5 = [(t4 * t1 + t17 * t20 + t3 * t2) * t65 + t17 * t62 + (t47 * t20 + t21 * t51) * t66 - 0.2e1 * t58 - 0.2e1 * t20 * mrSges(4,1) + 0.2e1 * m(4) * (t55 * t19 - t50 * t20) + t52 - 0.2e1 * (t2 * t27 + t3 * t24) * mrSges(6,3) + (t1 * t26 - t4 * t25) * t68 + t67 * t61 + t51 * t69; m(6) * (t1 * t27 + t2 * t26 + t4 * t24 - t3 * t25); (t70 - t71) * t65; -t58 + (t17 + t36) * t11 + (-mrSges(4,1) + t67) * t20 + m(6) * (t14 * t1 + t10 * t3 + t13 * t2 + t36 * t20 + t9 * t4) + m(5) * (-pkin(3) * t20 + t18 * t49 + t21 * t48) + (t51 + t48) * mrSges(5,3) + ((-t10 - t2) * t27 + (t1 + t9) * t26 - (t14 + t4) * t25 + (-t13 - t3) * t24) * mrSges(6,3) + t52; m(6) * (t10 * t26 - t13 * t25 + t14 * t24 + t9 * t27); t36 * t62 + (t13 * t10 + t14 * t9) * t65 + t49 * t66 * qJD(4) + t52 + t48 * t69 + (-t10 * t27 - t13 * t24 - t14 * t25 + t9 * t26) * t68; (m(6) / 0.2e1 + m(5) / 0.2e1) * t61 + t11; 0; t11; 0; t2 * mrSges(6,1) - t1 * mrSges(6,2) + t56; -t11; t10 * mrSges(6,1) - t9 * mrSges(6,2) + t56; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
