% Calculate time derivative of joint inertia matrix for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:52
% EndTime: 2019-12-31 17:43:53
% DurationCPUTime: 0.24s
% Computational Cost: add. (219->57), mult. (505->98), div. (0->0), fcn. (398->6), ass. (0->26)
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t34 = qJD(3) * (t19 ^ 2 + t20 ^ 2);
t27 = t19 * qJD(4);
t33 = 2 * m(6);
t22 = sin(qJ(5));
t23 = cos(qJ(5));
t25 = t19 * t22 + t20 * t23;
t13 = t19 * t23 - t20 * t22;
t9 = qJD(5) * t13;
t32 = t25 * t9;
t8 = t25 * qJD(5);
t31 = t13 * t8;
t16 = sin(pkin(7)) * pkin(1) + qJ(3);
t30 = -pkin(6) + t16;
t29 = t16 * t34;
t5 = t9 * mrSges(6,1) - t8 * mrSges(6,2);
t10 = t30 * t19;
t11 = t30 * t20;
t3 = t23 * t10 - t22 * t11;
t4 = t22 * t10 + t23 * t11;
t24 = cos(pkin(7)) * pkin(1) + t19 * qJ(4) + pkin(2);
t6 = (pkin(3) + pkin(4)) * t20 + t24;
t2 = t13 * qJD(3) - qJD(5) * t4;
t1 = t25 * qJD(3) + qJD(5) * t3;
t7 = [-0.2e1 * Ifges(6,1) * t31 + 0.2e1 * Ifges(6,2) * t32 + 0.2e1 * t6 * t5 + 0.2e1 * (t20 * mrSges(5,1) + mrSges(6,1) * t25 + t13 * mrSges(6,2) + t19 * mrSges(5,3)) * t27 + (t4 * t1 + t3 * t2 + t6 * t27) * t33 + 0.2e1 * m(5) * (-(-t20 * pkin(3) - t24) * t27 + t29) + 0.2e1 * m(4) * t29 + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * t34 + 0.2e1 * (-t13 * t9 + t25 * t8) * Ifges(6,4) + 0.2e1 * (-t1 * t25 - t2 * t13 + t3 * t8 - t4 * t9) * mrSges(6,3); m(6) * (t1 * t13 - t2 * t25 - t3 * t9 - t4 * t8); (-t31 + t32) * t33; (-m(5) - m(6)) * t27 - t5; 0; 0; m(6) * (t22 * t1 + t23 * t2 + (-t22 * t3 + t23 * t4) * qJD(5)) + m(5) * t19 * qJD(3) + (-t22 * t9 + t23 * t8 + (t13 * t22 - t23 * t25) * qJD(5)) * mrSges(6,3); m(6) * (-t22 * t8 - t23 * t9 + (t13 * t23 + t22 * t25) * qJD(5)); 0; 0; t2 * mrSges(6,1) - t1 * mrSges(6,2) - Ifges(6,5) * t8 - Ifges(6,6) * t9; -t5; 0; (-mrSges(6,1) * t22 - mrSges(6,2) * t23) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
