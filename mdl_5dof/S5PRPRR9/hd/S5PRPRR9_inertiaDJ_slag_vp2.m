% Calculate time derivative of joint inertia matrix for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR9_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR9_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:38
% EndTime: 2019-12-31 17:39:38
% DurationCPUTime: 0.26s
% Computational Cost: add. (213->59), mult. (466->97), div. (0->0), fcn. (242->4), ass. (0->33)
t23 = cos(qJ(4));
t20 = sin(qJ(5));
t22 = cos(qJ(5));
t35 = t20 ^ 2 + t22 ^ 2;
t47 = mrSges(6,3) * t35;
t49 = mrSges(5,2) - t47;
t50 = t49 * t23;
t48 = m(6) * pkin(7);
t10 = -t22 * mrSges(6,1) + t20 * mrSges(6,2);
t36 = mrSges(5,1) - t10;
t21 = sin(qJ(4));
t24 = -pkin(2) - pkin(3);
t9 = t23 * qJ(3) + t21 * t24;
t45 = t35 * t23;
t32 = qJD(5) * t22;
t33 = qJD(5) * t20;
t25 = -(t20 * (Ifges(6,4) * t20 + Ifges(6,2) * t22) - t22 * (Ifges(6,1) * t20 + Ifges(6,4) * t22)) * qJD(5) + t20 * (Ifges(6,1) * t32 - Ifges(6,4) * t33) + t22 * (Ifges(6,4) * t32 - Ifges(6,2) * t33);
t44 = -m(6) * pkin(4) - t36;
t43 = 0.2e1 * m(6);
t27 = mrSges(6,1) * t20 + mrSges(6,2) * t22;
t3 = t27 * qJD(5);
t42 = -0.2e1 * t3;
t8 = -t21 * qJ(3) + t23 * t24;
t1 = t23 * qJD(3) + qJD(4) * t8;
t41 = t21 * t1;
t2 = t21 * qJD(3) + t9 * qJD(4);
t40 = t23 * t2;
t39 = t23 * t3;
t34 = qJD(4) * t23;
t31 = t35 * t1;
t7 = -pkin(7) + t9;
t6 = pkin(4) - t8;
t4 = [0; 0; t6 * t42 + 0.2e1 * t1 * mrSges(5,2) + (t6 * t2 + t7 * t31) * t43 + 0.2e1 * m(5) * (t9 * t1 - t8 * t2) + t25 + 0.2e1 * t36 * t2 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * qJD(3) - 0.2e1 * t47 * t1; 0; t39 + m(6) * (t35 * t41 - t40) + m(5) * (-t40 + t41) + (t36 * t21 + t50 + m(6) * (t21 * t6 + t45 * t7) + m(5) * (-t21 * t8 + t23 * t9)) * qJD(4); (-0.1e1 + t35) * t21 * t34 * t43; 0; t31 * t48 + (pkin(4) + t6) * t3 + t44 * t2 - t49 * t1 - t25; -t39 + (t44 * t21 + t45 * t48 - t50) * qJD(4); pkin(4) * t42 + t25; t3; -t27 * t1 + ((-mrSges(6,1) * t7 - Ifges(6,5)) * t22 + (mrSges(6,2) * t7 + Ifges(6,6)) * t20) * qJD(5); (t21 * t33 - t22 * t34) * mrSges(6,2) + (-t20 * t34 - t21 * t32) * mrSges(6,1); (Ifges(6,5) * t22 - Ifges(6,6) * t20 + t10 * pkin(7)) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
