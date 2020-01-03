% Calculate time derivative of joint inertia matrix for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:38
% EndTime: 2019-12-31 17:47:39
% DurationCPUTime: 0.39s
% Computational Cost: add. (447->89), mult. (998->159), div. (0->0), fcn. (849->6), ass. (0->47)
t59 = pkin(3) + qJ(2);
t34 = cos(pkin(7));
t30 = t34 ^ 2;
t32 = sin(pkin(7));
t58 = t32 ^ 2 + t30;
t35 = sin(qJ(5));
t36 = cos(qJ(5));
t31 = sin(pkin(8));
t54 = t31 * t34;
t18 = t32 * t36 + t35 * t54;
t57 = 0.2e1 * t18;
t40 = -t32 * t35 + t36 * t54;
t56 = -0.2e1 * t40;
t45 = qJD(3) * t32;
t22 = -qJD(4) * t34 - t45;
t33 = cos(pkin(8));
t47 = qJD(2) * t32;
t11 = t22 * t31 - t33 * t47;
t55 = t31 * t11;
t53 = t33 * t11;
t52 = t33 * t34;
t16 = t18 * qJD(5);
t17 = t40 * qJD(5);
t51 = Ifges(6,5) * t16 + Ifges(6,6) * t17;
t43 = -qJ(3) * t32 - pkin(1);
t20 = (-pkin(2) - qJ(4)) * t34 + t43;
t23 = t59 * t32;
t50 = t33 * t20 + t31 * t23;
t49 = t58 * qJ(2) * qJD(2);
t48 = t59 * t34;
t46 = qJD(2) * t34;
t10 = (pkin(4) * t33 + pkin(6) * t31) * t34 + t48;
t7 = pkin(6) * t32 + t50;
t3 = t10 * t36 - t35 * t7;
t4 = t35 * t10 + t36 * t7;
t41 = -t20 * t31 + t23 * t33;
t39 = qJD(5) * (-mrSges(6,1) * t36 + mrSges(6,2) * t35);
t12 = t22 * t33 + t31 * t47;
t1 = qJD(5) * t3 + t12 * t36 + t35 * t46;
t2 = -qJD(5) * t4 - t35 * t12 + t36 * t46;
t38 = t1 * t36 - t2 * t35 + (-t3 * t36 - t35 * t4) * qJD(5);
t8 = -mrSges(6,2) * t52 + mrSges(6,3) * t18;
t9 = mrSges(6,1) * t52 + mrSges(6,3) * t40;
t37 = (-t35 * t8 - t36 * t9) * qJD(5) + (t16 * t35 + t17 * t36) * mrSges(6,3);
t6 = -pkin(4) * t32 - t41;
t5 = -mrSges(6,1) * t17 + mrSges(6,2) * t16;
t13 = [0.2e1 * t11 * (-t18 * mrSges(6,1) - mrSges(6,2) * t40) + 0.2e1 * t6 * t5 + 0.2e1 * t1 * t8 + 0.2e1 * t2 * t9 - 0.2e1 * t11 * (t32 * mrSges(5,1) + mrSges(5,3) * t54) + 0.2e1 * t12 * (-t32 * mrSges(5,2) - mrSges(5,3) * t52) + t51 * t52 - 0.2e1 * (mrSges(4,2) * t34 - mrSges(4,3) * t32) * t45 + (0.2e1 * mrSges(6,3) * t4 + Ifges(6,4) * t56 + Ifges(6,2) * t57 + Ifges(6,6) * t52) * t17 + (-0.2e1 * mrSges(6,3) * t3 + Ifges(6,1) * t56 + Ifges(6,4) * t57 + Ifges(6,5) * t52) * t16 + 0.2e1 * (t30 * (mrSges(5,1) * t33 - mrSges(5,2) * t31) + (mrSges(3,3) + mrSges(4,1)) * t58) * qJD(2) + 0.2e1 * m(6) * (t1 * t4 + t11 * t6 + t2 * t3) + 0.2e1 * m(5) * (-t41 * t11 + t50 * t12 + t48 * t46) + 0.2e1 * m(4) * (-(-pkin(2) * t34 + t43) * t45 + t49) + 0.2e1 * m(3) * t49; -m(4) * t45 + t31 * t5 + t37 * t33 + m(6) * (t38 * t33 + t55) + m(5) * (t12 * t33 + t55); 0; m(4) * t47 - t33 * t5 + t37 * t31 + m(6) * (t38 * t31 - t53) + m(5) * (t12 * t31 - t53); 0; 0; m(6) * (t35 * t1 + t2 * t36) + m(5) * t46 + (-t16 * t36 + t17 * t35) * mrSges(6,3) + (m(6) * (-t3 * t35 + t36 * t4) + t36 * t8 - t35 * t9) * qJD(5); 0; 0; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 + t51; t33 * t39; t31 * t39; (-mrSges(6,1) * t35 - mrSges(6,2) * t36) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11); t13(2), t13(3), t13(5), t13(8), t13(12); t13(4), t13(5), t13(6), t13(9), t13(13); t13(7), t13(8), t13(9), t13(10), t13(14); t13(11), t13(12), t13(13), t13(14), t13(15);];
Mq = res;
