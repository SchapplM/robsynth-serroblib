% Calculate time derivative of joint inertia matrix for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP4_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP4_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:10
% EndTime: 2019-12-31 17:15:11
% DurationCPUTime: 0.42s
% Computational Cost: add. (507->90), mult. (1219->143), div. (0->0), fcn. (944->4), ass. (0->43)
t55 = 2 * mrSges(4,3);
t54 = 2 * mrSges(5,3);
t37 = cos(qJ(3));
t53 = t37 * pkin(2);
t52 = -mrSges(4,1) - mrSges(5,1);
t36 = sin(qJ(2));
t48 = -pkin(6) - pkin(5);
t31 = t48 * t36;
t38 = cos(qJ(2));
t32 = t48 * t38;
t35 = sin(qJ(3));
t18 = t35 * t31 - t37 * t32;
t51 = pkin(2) * qJD(3);
t50 = qJD(2) + qJD(3);
t34 = -t38 * pkin(2) - pkin(1);
t49 = 0.2e1 * t34;
t47 = qJD(3) * t35;
t46 = qJD(3) * t37;
t45 = 0.2e1 * t38;
t44 = qJD(2) * t48;
t43 = (-mrSges(4,2) - mrSges(5,2)) * t37;
t17 = t37 * t31 + t35 * t32;
t42 = 2 * Ifges(4,4) + 2 * Ifges(5,4);
t24 = -t35 * t36 + t37 * t38;
t15 = t50 * t24;
t25 = t35 * t38 + t37 * t36;
t29 = t36 * t44;
t30 = t38 * t44;
t6 = -t18 * qJD(3) - t35 * t29 + t37 * t30;
t3 = -t15 * qJ(4) - t25 * qJD(4) + t6;
t41 = m(5) * t3 - t15 * mrSges(5,3);
t5 = t37 * t29 + t35 * t30 + t31 * t46 + t32 * t47;
t16 = t50 * t25;
t40 = -t35 * t16 + (t24 * t37 + t25 * t35) * qJD(3);
t2 = -t16 * qJ(4) + t24 * qJD(4) + t5;
t39 = t6 * mrSges(4,1) + t3 * mrSges(5,1) - t5 * mrSges(4,2) - t2 * mrSges(5,2) - (Ifges(4,6) + Ifges(5,6)) * t16 + (Ifges(4,5) + Ifges(5,5)) * t15;
t33 = pkin(3) + t53;
t19 = -t24 * pkin(3) + t34;
t10 = t15 * mrSges(5,2);
t9 = qJD(2) * t36 * pkin(2) + t16 * pkin(3);
t8 = t24 * qJ(4) + t18;
t7 = -t25 * qJ(4) + t17;
t1 = [0.2e1 * t9 * (-t24 * mrSges(5,1) + t25 * mrSges(5,2)) + 0.2e1 * t19 * t10 + 0.2e1 * m(4) * (t17 * t6 + t18 * t5) + 0.2e1 * m(5) * (t19 * t9 + t8 * t2 + t7 * t3) + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t38) * t45 + (m(4) * pkin(2) * t49 - 0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * pkin(2) * (-t24 * mrSges(4,1) + t25 * mrSges(4,2)) - 0.2e1 * Ifges(3,4) * t36 + (Ifges(3,1) - Ifges(3,2)) * t45) * t36) * qJD(2) + (t2 * t24 - t3 * t25) * t54 + (t5 * t24 - t6 * t25) * t55 - (-0.2e1 * t34 * mrSges(4,1) - 0.2e1 * t19 * mrSges(5,1) + t18 * t55 + t25 * t42 + t8 * t54 + 0.2e1 * (Ifges(4,2) + Ifges(5,2)) * t24) * t16 + (mrSges(4,2) * t49 - 0.2e1 * t17 * mrSges(4,3) - 0.2e1 * t7 * mrSges(5,3) + t24 * t42 + 0.2e1 * (Ifges(4,1) + Ifges(5,1)) * t25) * t15; t41 * t33 + (Ifges(3,5) * t38 - Ifges(3,6) * t36 + (-mrSges(3,1) * t38 + mrSges(3,2) * t36) * pkin(5)) * qJD(2) + (m(5) * (t2 * t35 + t46 * t8 - t47 * t7) + m(4) * (-t17 * t47 + t18 * t46 + t35 * t5 + t37 * t6) + t40 * mrSges(5,3) + (-t37 * t15 + t40) * mrSges(4,3)) * pkin(2) + t39; 0.2e1 * (t43 + ((-t33 + t53) * m(5) + t52) * t35) * t51; pkin(3) * t41 + t39; (t43 + (-m(5) * pkin(3) + t52) * t35) * t51; 0; m(5) * t9 + t16 * mrSges(5,1) + t10; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
