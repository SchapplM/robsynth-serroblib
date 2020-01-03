% Calculate time derivative of joint inertia matrix for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:41
% EndTime: 2019-12-31 19:27:42
% DurationCPUTime: 0.34s
% Computational Cost: add. (277->83), mult. (609->121), div. (0->0), fcn. (337->6), ass. (0->47)
t62 = Ifges(6,1) - Ifges(6,2);
t29 = cos(pkin(8));
t30 = sin(qJ(5));
t32 = cos(qJ(5));
t46 = t30 ^ 2 + t32 ^ 2;
t61 = t29 * t46;
t60 = t46 * mrSges(6,3);
t43 = qJD(5) * t30;
t31 = sin(qJ(2));
t33 = cos(qJ(2));
t45 = pkin(1) * qJD(2);
t59 = (-mrSges(3,2) * t33 + (-mrSges(3,1) - mrSges(4,1)) * t31) * t45;
t58 = 2 * m(4);
t57 = 2 * m(5);
t56 = 2 * m(6);
t37 = -mrSges(6,1) * t30 - mrSges(6,2) * t32;
t10 = t37 * qJD(5);
t55 = 0.2e1 * t10;
t20 = t33 * t45 + qJD(3);
t28 = sin(pkin(8));
t42 = t31 * t45;
t5 = t28 * t20 - t29 * t42;
t54 = t29 * t5;
t49 = t20 * mrSges(4,3);
t48 = t30 * mrSges(6,2);
t17 = t32 * mrSges(6,1) - t48;
t47 = mrSges(5,1) + t17;
t41 = -t33 * pkin(1) - pkin(2);
t21 = -pkin(3) + t41;
t22 = t31 * pkin(1) + qJ(3);
t4 = t28 * t21 + t29 * t22;
t34 = -pkin(2) - pkin(3);
t14 = t29 * qJ(3) + t28 * t34;
t44 = qJD(3) * t29;
t6 = t29 * t20 + t28 * t42;
t40 = t46 * t6;
t38 = 0.2e1 * t47;
t3 = t29 * t21 - t28 * t22;
t13 = -t28 * qJ(3) + t29 * t34;
t36 = (2 * mrSges(5,2)) - 0.2e1 * t60;
t35 = (0.2e1 * Ifges(6,4) * t32 + t62 * t30) * qJD(5) * t32 + (-0.2e1 * Ifges(6,4) * t30 + t62 * t32) * t43;
t24 = Ifges(6,6) * t43;
t9 = -pkin(7) + t14;
t8 = pkin(4) - t13;
t2 = -pkin(7) + t4;
t1 = pkin(4) - t3;
t7 = [0.2e1 * t49 + t1 * t55 + t5 * t38 + t36 * t6 + 0.2e1 * t59 + (t1 * t5 + t2 * t40) * t56 + (-t3 * t5 + t4 * t6) * t57 + (t22 * t20 + t41 * t42) * t58 + t35; t6 * mrSges(5,2) + t49 + t47 * t5 + (t8 + t1) * t10 + t59 + (t29 * mrSges(5,2) + t47 * t28 + mrSges(4,3)) * qJD(3) + m(6) * (t8 * t5 + t9 * t40 + (t1 * t28 + t2 * t61) * qJD(3)) + m(5) * (-t13 * t5 + t14 * t6 + (-t28 * t3 + t29 * t4) * qJD(3)) + m(4) * (-pkin(2) * t42 + qJ(3) * t20 + qJD(3) * t22) + (-t44 - t6) * t60 + t35; t8 * t55 + (qJ(3) * t58 + 0.2e1 * mrSges(4,3) + t28 * t38 + t36 * t29 + (t28 * t8 + t9 * t61) * t56 + (-t13 * t28 + t14 * t29) * t57) * qJD(3) + t35; m(4) * t42 - t29 * t10 + m(6) * (t28 * t40 - t54) + m(5) * (t28 * t6 - t54); (-t10 + m(6) * (-0.1e1 + t46) * t28 * qJD(3)) * t29; 0; 0; 0; 0; 0; t24 + t37 * t6 + (-Ifges(6,5) * t32 - t17 * t2) * qJD(5); t24 + t37 * t44 + (t9 * t48 + (-mrSges(6,1) * t9 - Ifges(6,5)) * t32) * qJD(5); -t17 * t28 * qJD(5); t10; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
