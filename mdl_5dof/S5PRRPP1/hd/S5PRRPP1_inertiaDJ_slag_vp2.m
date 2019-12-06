% Calculate time derivative of joint inertia matrix for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:13
% EndTime: 2019-12-05 16:06:14
% DurationCPUTime: 0.42s
% Computational Cost: add. (356->86), mult. (895->131), div. (0->0), fcn. (661->4), ass. (0->44)
t46 = (mrSges(6,2) + mrSges(5,3));
t54 = 2 * t46;
t53 = m(5) + m(6);
t25 = sin(pkin(8));
t49 = t25 * pkin(3);
t21 = qJ(5) + t49;
t52 = m(6) * t21 + mrSges(6,3);
t27 = cos(qJ(3));
t24 = -t27 * pkin(3) - pkin(2);
t51 = 0.2e1 * t24;
t50 = m(5) * pkin(3);
t44 = cos(pkin(8));
t26 = sin(qJ(3));
t47 = t25 * t26;
t31 = t27 * t44 - t47;
t16 = t31 * qJD(3);
t48 = t16 * mrSges(6,2);
t45 = -qJ(4) - pkin(6);
t43 = qJD(3) * t26;
t42 = qJD(3) * t27;
t35 = t44 * t26;
t18 = t25 * t27 + t35;
t41 = t18 * qJD(5);
t40 = 0.2e1 * t27;
t39 = pkin(3) * t43;
t20 = t45 * t27;
t10 = -t20 * t44 + t45 * t47;
t34 = qJD(3) * t45;
t14 = t27 * qJD(4) + t26 * t34;
t28 = -t26 * qJD(4) + t27 * t34;
t5 = t25 * t14 - t28 * t44;
t6 = t14 * t44 + t25 * t28;
t9 = -t25 * t20 - t35 * t45;
t38 = t10 * t6 + t9 * t5;
t36 = t44 * pkin(3);
t33 = -2 * Ifges(5,4) + 2 * Ifges(6,5);
t15 = t18 * qJD(3);
t12 = t15 * mrSges(6,1);
t13 = t16 * mrSges(5,2);
t30 = -t15 * mrSges(5,1) + t16 * mrSges(6,3) - t12 - t13;
t23 = -t36 - pkin(4);
t7 = -pkin(4) * t31 - t18 * qJ(5) + t24;
t2 = t15 * pkin(4) - t16 * qJ(5) + t39 - t41;
t1 = [0.2e1 * t53 * (-t15 * t31 + t18 * t16); t53 * (t10 * t16 + t9 * t15 + t6 * t18 - t31 * t5); 0.2e1 * t2 * (-mrSges(6,1) * t31 - t18 * mrSges(6,3)) + t13 * t51 + 0.2e1 * t7 * t12 + 0.2e1 * m(5) * t38 + 0.2e1 * m(6) * (t7 * t2 + t38) + (-0.2e1 * t7 * mrSges(6,3) - t31 * t33 + t9 * t54 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t18) * t16 + (mrSges(5,1) * t51 + t18 * t33 - 0.2e1 * (Ifges(5,2) + Ifges(6,3)) * t31 - 0.2e1 * t46 * t10) * t15 + ((-pkin(2) * mrSges(4,2) + Ifges(4,4) * t27) * t40 + (-0.2e1 * pkin(2) * mrSges(4,1) + 0.2e1 * pkin(3) * (-mrSges(5,1) * t31 + t18 * mrSges(5,2)) + t50 * t51 - 0.2e1 * Ifges(4,4) * t26 + (Ifges(4,1) - Ifges(4,2)) * t40) * t26) * qJD(3) + (t5 * t18 + t31 * t6) * t54; -mrSges(4,2) * t42 - mrSges(4,1) * t43 + (-t15 * t44 + t16 * t25) * t50 + m(6) * (t23 * t15 + t21 * t16 + t41) + t30; Ifges(4,5) * t42 - Ifges(4,6) * t43 + t23 * t48 + (m(6) * t10 + mrSges(6,2) * t31) * qJD(5) + (t25 * t50 - mrSges(5,2) + t52) * t6 + (m(6) * t23 - t44 * t50 - mrSges(5,1) - mrSges(6,1)) * t5 + (-mrSges(4,1) * t42 + mrSges(4,2) * t43) * pkin(6) + (-mrSges(5,3) * t36 + Ifges(6,4) + Ifges(5,5)) * t16 + (-mrSges(6,2) * t21 - mrSges(5,3) * t49 - Ifges(5,6) + Ifges(6,6)) * t15; 0.2e1 * t52 * qJD(5); 0; m(5) * t39 + m(6) * t2 - t30; 0; 0; m(6) * t15; m(6) * t5 + t48; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
