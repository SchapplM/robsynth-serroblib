% Calculate time derivative of joint inertia matrix for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:00
% EndTime: 2019-12-05 17:38:01
% DurationCPUTime: 0.38s
% Computational Cost: add. (411->80), mult. (823->130), div. (0->0), fcn. (601->4), ass. (0->46)
t28 = sin(qJ(4));
t29 = cos(qJ(5));
t30 = cos(qJ(4));
t52 = sin(qJ(5));
t15 = -t52 * t28 + t29 * t30;
t47 = qJD(4) * t30;
t48 = qJD(4) * t28;
t60 = -mrSges(5,1) * t48 - mrSges(5,2) * t47;
t59 = qJD(4) + qJD(5);
t24 = t28 ^ 2;
t25 = t30 ^ 2;
t57 = m(6) * pkin(4);
t56 = mrSges(6,3) * pkin(4);
t14 = -t29 * t28 - t30 * t52;
t7 = t15 * t59;
t55 = t14 * t7;
t8 = t59 * t14;
t54 = t15 * t8;
t26 = qJ(2) - pkin(6);
t53 = pkin(7) - t26;
t27 = pkin(1) + qJ(3);
t49 = t24 + t25;
t46 = t27 * qJD(3);
t45 = t28 * qJD(2);
t44 = t30 * qJD(2);
t43 = 2 * mrSges(6,3);
t42 = t29 * t56;
t39 = t8 * mrSges(6,1) - t7 * mrSges(6,2);
t17 = t53 * t30;
t37 = t52 * t56;
t36 = t7 * mrSges(6,1) + t8 * mrSges(6,2);
t35 = -t54 + t55;
t11 = t48 * t53 + t44;
t12 = -qJD(4) * t17 + t45;
t16 = t53 * t28;
t9 = t16 * t52 - t29 * t17;
t2 = qJD(5) * t9 + t11 * t52 + t29 * t12;
t32 = t29 * t16 + t17 * t52;
t3 = qJD(5) * t32 + t29 * t11 - t12 * t52;
t34 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t8 - Ifges(6,6) * t7;
t33 = t30 * mrSges(5,1) - t28 * mrSges(5,2);
t31 = -t14 * t2 + t15 * t3 - t32 * t7 + t8 * t9;
t21 = t28 * pkin(4) + t27;
t18 = pkin(4) * t47 + qJD(3);
t13 = (-mrSges(6,1) * t52 - mrSges(6,2) * t29) * qJD(5) * pkin(4);
t1 = [0.2e1 * Ifges(6,1) * t54 - 0.2e1 * Ifges(6,2) * t55 + 0.2e1 * t18 * (-t14 * mrSges(6,1) + t15 * mrSges(6,2)) + 0.2e1 * t21 * t36 + 0.2e1 * m(6) * (t21 * t18 - t2 * t32 + t9 * t3) + 0.2e1 * m(5) * (qJD(2) * t26 * t49 + t46) + 0.2e1 * m(4) * (qJ(2) * qJD(2) + t46) + 0.2e1 * (t8 * t14 - t15 * t7) * Ifges(6,4) - t31 * t43 + 0.2e1 * (t27 * t33 + (-t25 + t24) * Ifges(5,4)) * qJD(4) + 0.2e1 * (m(3) * qJ(2) - mrSges(5,3) * t49 + mrSges(4,2) + mrSges(3,3)) * qJD(2) + 0.2e1 * (-Ifges(5,1) + Ifges(5,2)) * t28 * t47 + 0.2e1 * (t28 * mrSges(5,1) + t30 * mrSges(5,2) + mrSges(4,3)) * qJD(3); -m(6) * t18 - t33 * qJD(4) + (-m(5) - m(4)) * qJD(3) - t36; 0; m(6) * t31 + t35 * t43 + (m(5) * t49 + m(4)) * qJD(2); 0; -0.2e1 * m(6) * t35; -t7 * t37 - t8 * t42 - mrSges(5,2) * t45 + mrSges(5,1) * t44 + (t52 * t2 + t29 * t3) * t57 - Ifges(5,5) * t48 - Ifges(5,6) * t47 + t34 + t60 * t26 + (t14 * t42 + t15 * t37 + (-t29 * t32 - t52 * t9) * t57) * qJD(5); 0; (t52 * t7 + t29 * t8 + (-t14 * t29 - t15 * t52) * qJD(5)) * t57 + t39 + t60; 0.2e1 * t13; t34; 0; t39; t13; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
