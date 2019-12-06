% Calculate time derivative of joint inertia matrix for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:49
% EndTime: 2019-12-05 16:39:51
% DurationCPUTime: 0.47s
% Computational Cost: add. (281->101), mult. (706->138), div. (0->0), fcn. (376->4), ass. (0->52)
t34 = sin(qJ(4));
t69 = -0.2e1 * t34;
t80 = Ifges(6,4) + Ifges(5,4);
t79 = Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,2);
t36 = cos(qJ(4));
t77 = t34 ^ 2 + t36 ^ 2;
t68 = 0.2e1 * t36;
t37 = cos(qJ(3));
t74 = t77 * t37;
t35 = sin(qJ(3));
t73 = (t77 * mrSges(5,3) - mrSges(4,2)) * t37 + (-t36 * mrSges(5,1) + t34 * mrSges(5,2) - mrSges(4,1)) * t35;
t48 = qJD(4) * t36;
t49 = qJD(4) * t34;
t38 = (t79 * t34 + t80 * t68) * t48 + (t79 * t36 + t80 * t69) * t49;
t72 = 2 * m(6);
t13 = mrSges(6,1) * t49 + mrSges(6,2) * t48;
t71 = 0.2e1 * t13;
t18 = -t36 * mrSges(6,1) + t34 * mrSges(6,2);
t70 = 0.2e1 * t18;
t67 = m(6) * pkin(4);
t66 = mrSges(6,3) * pkin(4);
t65 = t37 * pkin(2);
t64 = mrSges(5,2) * t36;
t54 = -Ifges(5,6) - Ifges(6,6);
t53 = -qJ(5) - pkin(7);
t52 = (Ifges(5,5) + Ifges(6,5)) * t48;
t51 = pkin(2) * qJD(3);
t23 = t35 * pkin(2) + pkin(7);
t50 = -qJ(5) - t23;
t47 = 0.2e1 * qJD(4);
t46 = t37 * t51;
t45 = pkin(4) * t49;
t44 = mrSges(6,1) + t67;
t25 = -t36 * pkin(4) - pkin(3);
t43 = qJD(4) * t53;
t42 = qJD(4) * t50;
t39 = mrSges(5,1) * t34 + t64;
t31 = t36 * qJ(5);
t30 = t36 * qJD(5);
t24 = -pkin(3) - t65;
t20 = t36 * pkin(7) + t31;
t17 = t53 * t34;
t16 = t25 - t65;
t15 = t35 * t51 + t45;
t14 = t39 * qJD(4);
t10 = t36 * t23 + t31;
t9 = t50 * t34;
t4 = -t34 * qJD(5) + t36 * t43;
t3 = t34 * t43 + t30;
t2 = (-qJD(5) - t46) * t34 + t36 * t42;
t1 = t34 * t42 + t36 * t46 + t30;
t5 = [0; m(6) * (t1 * t34 + t2 * t36 + (t10 * t36 - t34 * t9) * qJD(4)); 0.2e1 * t24 * t14 + t15 * t70 + t16 * t71 + (t10 * t1 + t16 * t15 + t9 * t2) * t72 + (t1 * t68 + t2 * t69 + (-t10 * t34 - t36 * t9) * t47) * mrSges(6,3) + 0.2e1 * (m(5) * (t74 * t23 + t24 * t35) + t73) * t51 + t38; m(6) * (t3 * t34 + t4 * t36 + (-t17 * t34 + t20 * t36) * qJD(4)); m(6) * (t20 * t1 + t3 * t10 + t25 * t15 + t16 * t45 + t17 * t2 + t4 * t9) + (t15 + t45) * t18 + (-pkin(3) + t24) * t14 + (t16 + t25) * t13 + (m(5) * (-pkin(3) * t35 + t74 * pkin(7)) + t73) * t51 + ((t1 + t3) * t36 + (-t2 - t4) * t34 + ((-t17 - t9) * t36 + (-t10 - t20) * t34) * qJD(4)) * mrSges(6,3) + t38; -0.2e1 * pkin(3) * t14 + t45 * t70 + t25 * t71 + (t17 * t4 + t20 * t3 + t25 * t45) * t72 + (t3 * t68 + t4 * t69 + (-t17 * t36 - t20 * t34) * t47) * mrSges(6,3) + t38; (-t64 + (-mrSges(5,1) - t67) * t34) * qJD(4) - t13; -t1 * mrSges(6,2) + t44 * t2 - t39 * t46 + ((-mrSges(5,1) * t23 - t66) * t36 + (mrSges(5,2) * t23 + t54) * t34) * qJD(4) + t52; -t3 * mrSges(6,2) + t44 * t4 + ((-mrSges(5,1) * pkin(7) - t66) * t36 + (mrSges(5,2) * pkin(7) + t54) * t34) * qJD(4) + t52; 0; 0; m(6) * t15 + t13; m(6) * t45 + t13; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
