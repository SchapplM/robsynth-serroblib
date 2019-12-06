% Calculate time derivative of joint inertia matrix for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:35
% EndTime: 2019-12-05 15:37:37
% DurationCPUTime: 0.41s
% Computational Cost: add. (417->85), mult. (1160->132), div. (0->0), fcn. (959->6), ass. (0->46)
t54 = (mrSges(6,2) + mrSges(5,3));
t63 = -2 * t54;
t62 = m(5) + m(6);
t32 = sin(pkin(8));
t33 = cos(pkin(8));
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t61 = -t34 * t32 + t36 * t33;
t60 = m(6) * qJD(5);
t59 = m(6) * qJ(5) + mrSges(6,3);
t58 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t57 = -mrSges(5,2) + t59;
t53 = pkin(6) + qJ(3);
t20 = t61 * qJD(4);
t24 = t36 * t32 + t34 * t33;
t21 = t24 * qJD(4);
t10 = t21 * mrSges(6,1) - t20 * mrSges(6,3);
t11 = t21 * mrSges(5,1) + t20 * mrSges(5,2);
t52 = -t10 - t11;
t51 = t32 ^ 2 + t33 ^ 2;
t35 = sin(qJ(2));
t50 = qJD(2) * t35;
t37 = cos(qJ(2));
t49 = qJD(2) * t37;
t48 = qJD(4) * t35;
t47 = t35 * t49;
t29 = -t33 * pkin(3) - pkin(2);
t26 = t53 * t32;
t27 = t53 * t33;
t15 = -t34 * t26 + t36 * t27;
t42 = -t36 * t26 - t34 * t27;
t6 = qJD(3) * t61 + qJD(4) * t42;
t7 = qJD(3) * t24 + qJD(4) * t15;
t46 = t15 * t6 - t42 * t7;
t45 = t51 * mrSges(4,3);
t44 = t51 * qJ(3);
t43 = -2 * Ifges(5,4) + 2 * Ifges(6,5);
t16 = t24 * t35;
t17 = t61 * t35;
t8 = -t24 * t48 + t49 * t61;
t9 = t24 * t49 + t61 * t48;
t41 = t15 * t8 + t7 * t16 + t6 * t17 - t42 * t9;
t13 = -mrSges(6,1) * t61 - t24 * mrSges(6,3);
t12 = -pkin(4) * t61 - t24 * qJ(5) + t29;
t5 = t21 * pkin(4) - t20 * qJ(5) - t24 * qJD(5);
t1 = [0.2e1 * m(4) * (-0.1e1 + t51) * t47 + 0.2e1 * t62 * (t16 * t9 + t17 * t8 - t47); t52 * t37 + ((-mrSges(3,2) + t45) * t37 + (-t33 * mrSges(4,1) - mrSges(5,1) * t61 + t32 * mrSges(4,2) + t24 * mrSges(5,2) - mrSges(3,1) + t13) * t35) * qJD(2) + m(4) * (t51 * t35 * qJD(3) + (-pkin(2) * t35 + t37 * t44) * qJD(2)) + m(5) * (t29 * t50 + t41) + m(6) * (t12 * t50 - t5 * t37 + t41) + t54 * (t16 * t20 - t17 * t21 + t9 * t24 + t61 * t8); 0.2e1 * t12 * t10 + 0.2e1 * t29 * t11 + 0.2e1 * t5 * t13 + 0.2e1 * m(5) * t46 + 0.2e1 * m(6) * (t12 * t5 + t46) + 0.2e1 * (m(4) * t44 + t45) * qJD(3) + (t24 * t43 - 0.2e1 * (Ifges(5,2) + Ifges(6,3)) * t61 + t15 * t63) * t21 + (-t61 * t43 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t24 + t42 * t63) * t20 + 0.2e1 * t54 * (t7 * t24 + t6 * t61); (m(4) + t62) * t50; m(6) * t5 - t52; 0; t17 * t60 + t57 * t8 + t58 * t9; t15 * t60 + (-Ifges(5,6) + Ifges(6,6)) * t21 + (Ifges(6,4) + Ifges(5,5)) * t20 + (-pkin(4) * t20 - qJ(5) * t21 + qJD(5) * t61) * mrSges(6,2) + t58 * t7 + t57 * t6; 0; 0.2e1 * t59 * qJD(5); m(6) * t9; m(6) * t7 + t20 * mrSges(6,2); 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
