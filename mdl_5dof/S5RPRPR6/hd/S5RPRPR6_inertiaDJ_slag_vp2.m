% Calculate time derivative of joint inertia matrix for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:42
% EndTime: 2019-12-31 18:17:43
% DurationCPUTime: 0.22s
% Computational Cost: add. (194->53), mult. (453->75), div. (0->0), fcn. (252->6), ass. (0->31)
t14 = cos(pkin(8)) * pkin(1) + pkin(2);
t19 = sin(qJ(3));
t21 = cos(qJ(3));
t37 = pkin(1) * sin(pkin(8));
t24 = t19 * t14 + t21 * t37;
t3 = t24 * qJD(3);
t43 = 0.2e1 * t3;
t20 = cos(qJ(5));
t42 = t20 ^ 2;
t28 = t21 * t14 - t19 * t37;
t41 = -mrSges(4,1) + mrSges(5,2);
t18 = sin(qJ(5));
t10 = t18 * mrSges(6,1) + t20 * mrSges(6,2);
t40 = mrSges(5,3) + t10;
t26 = mrSges(6,1) * t20 - mrSges(6,2) * t18;
t7 = t26 * qJD(5);
t39 = 0.2e1 * t7;
t2 = t28 * qJD(3);
t1 = -qJD(4) - t2;
t5 = qJ(4) + t24;
t38 = t5 * t1;
t36 = t2 * mrSges(4,2);
t32 = t18 ^ 2 + t42;
t30 = t32 * t3;
t29 = t32 * mrSges(6,3);
t27 = -pkin(3) - t28;
t25 = -qJ(4) * t1 + qJD(4) * t5;
t23 = (-0.2e1 * t42 * Ifges(6,4) + (0.2e1 * Ifges(6,4) * t18 + 0.2e1 * (-Ifges(6,1) + Ifges(6,2)) * t20) * t18) * qJD(5);
t22 = -pkin(3) - pkin(7);
t4 = -pkin(7) + t27;
t6 = [-0.2e1 * t36 + t5 * t39 + t41 * t43 + 0.2e1 * m(5) * (t27 * t3 - t38) + 0.2e1 * m(6) * (t4 * t30 - t38) + 0.2e1 * m(4) * (t24 * t2 - t28 * t3) + t23 - 0.2e1 * t40 * t1 - 0.2e1 * t29 * t3; 0; 0; -t36 + (qJ(4) + t5) * t7 + (-t29 + t41) * t3 + m(5) * (-pkin(3) * t3 + t25) + m(6) * (t22 * t30 + t25) + t23 + t40 * (qJD(4) - t1); 0; qJ(4) * t39 + 0.2e1 * ((m(5) + m(6)) * qJ(4) + t40) * qJD(4) + t23; (m(5) / 0.2e1 + m(6) * t32 / 0.2e1) * t43; 0; 0; 0; t26 * t3 + ((-mrSges(6,2) * t4 - Ifges(6,6)) * t20 + (-mrSges(6,1) * t4 - Ifges(6,5)) * t18) * qJD(5); -t7; ((-mrSges(6,2) * t22 - Ifges(6,6)) * t20 + (-mrSges(6,1) * t22 - Ifges(6,5)) * t18) * qJD(5); -t10 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t6(1), t6(2), t6(4), t6(7), t6(11); t6(2), t6(3), t6(5), t6(8), t6(12); t6(4), t6(5), t6(6), t6(9), t6(13); t6(7), t6(8), t6(9), t6(10), t6(14); t6(11), t6(12), t6(13), t6(14), t6(15);];
Mq = res;
