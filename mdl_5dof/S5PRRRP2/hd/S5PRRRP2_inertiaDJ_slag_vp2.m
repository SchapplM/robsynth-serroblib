% Calculate time derivative of joint inertia matrix for
% S5PRRRP2
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:38
% EndTime: 2019-12-05 16:41:40
% DurationCPUTime: 0.33s
% Computational Cost: add. (196->71), mult. (538->98), div. (0->0), fcn. (257->4), ass. (0->40)
t31 = sin(qJ(4));
t33 = cos(qJ(4));
t66 = t31 ^ 2 + t33 ^ 2;
t65 = 2 * Ifges(5,4) - 2 * Ifges(6,5);
t34 = cos(qJ(3));
t52 = pkin(2) * qJD(3);
t47 = t34 * t52;
t64 = t66 * t47;
t63 = Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3);
t50 = qJD(4) * t33;
t62 = (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t66) * t47;
t51 = qJD(4) * t31;
t2 = -pkin(4) * t51 + qJ(5) * t50 + t31 * qJD(5);
t61 = m(6) * t2;
t60 = t34 * pkin(2);
t32 = sin(qJ(3));
t21 = t32 * pkin(2) + pkin(7);
t55 = t64 * t21;
t54 = t64 * pkin(7);
t49 = qJD(5) * t33;
t48 = m(6) * t49;
t46 = t32 * t52;
t45 = mrSges(6,2) * t49 + Ifges(6,6) * t51 + (Ifges(6,4) + Ifges(5,5)) * t50;
t17 = -t33 * mrSges(5,1) + t31 * mrSges(5,2);
t41 = t31 * mrSges(5,1) + t33 * mrSges(5,2);
t16 = -t33 * mrSges(6,1) - t31 * mrSges(6,3);
t40 = t31 * mrSges(6,1) - t33 * mrSges(6,3);
t39 = -t33 * pkin(4) - t31 * qJ(5);
t38 = (-mrSges(4,1) + t17) * t46;
t15 = -pkin(3) + t39;
t37 = t39 * mrSges(6,2) - Ifges(5,6) * t31;
t36 = (-t65 * t31 + t63 * t33) * t51 + (t63 * t31 + t65 * t33) * t50;
t35 = m(6) * t39 + t16 + t17;
t24 = mrSges(6,2) * t50;
t22 = -pkin(3) - t60;
t13 = t41 * qJD(4);
t12 = t40 * qJD(4);
t6 = t15 - t60;
t1 = -t2 + t46;
t3 = [0; 0; 0.2e1 * t1 * t16 + 0.2e1 * t6 * t12 + 0.2e1 * t22 * t13 + 0.2e1 * t38 + 0.2e1 * m(6) * (t6 * t1 + t55) + 0.2e1 * m(5) * (t22 * t46 + t55) + 0.2e1 * t62 + t36; 0; (t1 - t2) * t16 + (t22 - pkin(3)) * t13 + (t6 + t15) * t12 + t38 + m(6) * (t15 * t1 - t2 * t6 + t54) + m(5) * (-pkin(3) * t46 + t54) + t62 + t36; -0.2e1 * pkin(3) * t13 - 0.2e1 * t2 * t16 + 0.2e1 * (t12 - t61) * t15 + t36; t61 + ((-mrSges(5,2) + mrSges(6,3)) * t33 + (-mrSges(5,1) - mrSges(6,1)) * t31) * qJD(4); t21 * t48 + (m(6) * (-pkin(4) * t31 + qJ(5) * t33) - t40 - t41) * t47 + (t35 * t21 + t37) * qJD(4) + t45; t37 * qJD(4) + (t35 * qJD(4) + t48) * pkin(7) + t45; 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); m(6) * t51; t24 + (t21 * t50 + t31 * t47) * m(6); m(6) * pkin(7) * t50 + t24; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
