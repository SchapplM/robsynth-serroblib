% Calculate time derivative of joint inertia matrix for
% S4RRRP3
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
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP3_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP3_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:01
% EndTime: 2019-12-31 17:14:02
% DurationCPUTime: 0.29s
% Computational Cost: add. (190->67), mult. (518->93), div. (0->0), fcn. (249->4), ass. (0->40)
t31 = sin(qJ(3));
t33 = cos(qJ(3));
t66 = t31 ^ 2 + t33 ^ 2;
t65 = 2 * Ifges(4,4) - 2 * Ifges(5,5);
t34 = cos(qJ(2));
t52 = pkin(1) * qJD(2);
t47 = t34 * t52;
t64 = t66 * t47;
t63 = Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,3);
t50 = qJD(3) * t33;
t62 = (-mrSges(3,2) + (mrSges(5,2) + mrSges(4,3)) * t66) * t47;
t16 = -t33 * mrSges(5,1) - t31 * mrSges(5,3);
t61 = 0.2e1 * t16;
t60 = t34 * pkin(1);
t32 = sin(qJ(2));
t21 = pkin(1) * t32 + pkin(6);
t55 = t64 * t21;
t54 = t64 * pkin(6);
t51 = qJD(3) * t31;
t49 = qJD(4) * t33;
t48 = m(5) * t49;
t46 = t32 * t52;
t45 = mrSges(5,2) * t49 + Ifges(5,6) * t51 + (Ifges(5,4) + Ifges(4,5)) * t50;
t17 = -t33 * mrSges(4,1) + t31 * mrSges(4,2);
t41 = t31 * mrSges(4,1) + t33 * mrSges(4,2);
t40 = t31 * mrSges(5,1) - t33 * mrSges(5,3);
t39 = -pkin(3) * t33 - qJ(4) * t31;
t38 = (-mrSges(3,1) + t17) * t46;
t15 = -pkin(2) + t39;
t2 = pkin(3) * t51 - qJ(4) * t50 - qJD(4) * t31;
t37 = mrSges(5,2) * t39 - Ifges(4,6) * t31;
t36 = (-t31 * t65 + t63 * t33) * t51 + (t63 * t31 + t33 * t65) * t50;
t35 = m(5) * t39 + t16 + t17;
t24 = mrSges(5,2) * t50;
t22 = -pkin(2) - t60;
t13 = t41 * qJD(3);
t12 = t40 * qJD(3);
t6 = t15 - t60;
t1 = t2 + t46;
t3 = [t1 * t61 + 0.2e1 * t6 * t12 + 0.2e1 * t22 * t13 + 0.2e1 * t38 + 0.2e1 * m(4) * (t22 * t46 + t55) + 0.2e1 * m(5) * (t1 * t6 + t55) + 0.2e1 * t62 + t36; (t1 + t2) * t16 + (t22 - pkin(2)) * t13 + (t6 + t15) * t12 + t38 + m(4) * (-pkin(2) * t46 + t54) + m(5) * (t1 * t15 + t2 * t6 + t54) + t62 + t36; -0.2e1 * pkin(2) * t13 + t2 * t61 + 0.2e1 * (m(5) * t2 + t12) * t15 + t36; t21 * t48 + (m(5) * (-pkin(3) * t31 + qJ(4) * t33) - t40 - t41) * t47 + (t21 * t35 + t37) * qJD(3) + t45; t37 * qJD(3) + (qJD(3) * t35 + t48) * pkin(6) + t45; 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4); t24 + (t21 * t50 + t31 * t47) * m(5); m(5) * pkin(6) * t50 + t24; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;
