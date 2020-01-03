% Calculate time derivative of joint inertia matrix for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:27
% EndTime: 2019-12-31 17:34:28
% DurationCPUTime: 0.24s
% Computational Cost: add. (119->59), mult. (355->95), div. (0->0), fcn. (190->4), ass. (0->30)
t42 = m(6) * pkin(4);
t17 = sin(qJ(4));
t29 = qJD(4) * t17;
t19 = cos(qJ(4));
t27 = qJD(4) * t19;
t4 = mrSges(6,1) * t29 + mrSges(6,2) * t27;
t25 = t29 * t42 + t4;
t5 = mrSges(5,1) * t29 + mrSges(5,2) * t27;
t41 = t25 + t5;
t31 = t17 ^ 2 + t19 ^ 2;
t32 = -qJ(5) - pkin(6);
t6 = t32 * t17;
t8 = t32 * t19;
t40 = t17 * t6 + t19 * t8;
t20 = cos(qJ(3));
t30 = qJD(3) * t20;
t10 = -t19 * pkin(4) - pkin(3);
t39 = m(6) * t10 - t19 * mrSges(6,1) + t17 * mrSges(6,2);
t38 = -2 * mrSges(6,3);
t34 = mrSges(5,2) + mrSges(6,2);
t33 = Ifges(5,4) + Ifges(6,4);
t18 = sin(qJ(3));
t28 = qJD(4) * t18;
t26 = 0.2e1 * t17;
t23 = mrSges(6,1) + t42;
t22 = qJD(4) * t32;
t21 = -mrSges(5,1) - t23;
t3 = -t17 * qJD(5) + t19 * t22;
t2 = t19 * qJD(5) + t17 * t22;
t1 = [0; 0; 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (-0.1e1 + t31) * t18 * t30; m(6) * (t40 * qJD(4) - t2 * t17 - t3 * t19); m(6) * t17 * t28 * t8 + m(5) * t31 * pkin(6) * t30 + (m(6) * (-t17 * t3 + t19 * t2 - t27 * t6) + (-m(5) * pkin(3) - t19 * mrSges(5,1) + t17 * mrSges(5,2) - mrSges(4,1) + t39) * qJD(3)) * t18 + ((-m(6) * t40 - mrSges(4,2) + (mrSges(5,3) + mrSges(6,3)) * t31) * qJD(3) - t41) * t20; 0.2e1 * t10 * t4 + 0.2e1 * m(6) * (-t8 * t2 + t6 * t3) - 0.2e1 * pkin(3) * t5 + (t3 * t38 + (0.2e1 * t39 * pkin(4) - t33 * t26 - t8 * t38) * qJD(4)) * t17 + (0.2e1 * t2 * mrSges(6,3) + (t6 * t38 + 0.2e1 * t33 * t19 + (Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,2)) * t26) * qJD(4)) * t19; t41; (t17 * t34 + t19 * t21) * t28 + (t17 * t21 - t19 * t34) * t30; -t2 * mrSges(6,2) + t23 * t3 + ((mrSges(5,2) * pkin(6) - Ifges(5,6) - Ifges(6,6)) * t17 + (-mrSges(5,1) * pkin(6) - mrSges(6,3) * pkin(4) + Ifges(5,5) + Ifges(6,5)) * t19) * qJD(4); 0; 0; m(6) * qJD(3) * t18; t25; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
