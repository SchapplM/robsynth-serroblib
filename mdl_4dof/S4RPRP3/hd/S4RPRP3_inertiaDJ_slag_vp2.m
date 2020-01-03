% Calculate time derivative of joint inertia matrix for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP3_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP3_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:32
% EndTime: 2019-12-31 16:42:33
% DurationCPUTime: 0.13s
% Computational Cost: add. (92->41), mult. (208->61), div. (0->0), fcn. (116->4), ass. (0->19)
t13 = -cos(pkin(6)) * pkin(1) - pkin(2);
t22 = 0.2e1 * t13;
t11 = cos(qJ(3));
t21 = -0.2e1 * t11 * pkin(3) + 0.2e1 * t13;
t20 = -2 * mrSges(5,3);
t19 = m(5) * pkin(3);
t10 = sin(qJ(3));
t15 = qJD(3) * t10;
t18 = qJD(3) * t11 * mrSges(5,2) + mrSges(5,1) * t15;
t17 = Ifges(5,4) + Ifges(4,4);
t6 = sin(pkin(6)) * pkin(1) + pkin(5);
t16 = qJ(4) + t6;
t14 = 0.2e1 * t11;
t12 = qJD(3) * t16;
t4 = t16 * t11;
t3 = t16 * t10;
t2 = -t10 * qJD(4) - t11 * t12;
t1 = t11 * qJD(4) - t10 * t12;
t5 = [t18 * t21 + 0.2e1 * m(5) * (t4 * t1 - t3 * t2) + 0.2e1 * (t1 * t11 - t2 * t10) * mrSges(5,3) + ((mrSges(4,2) * t22 + t17 * t14 - t3 * t20) * t11 + (t4 * t20 + mrSges(4,1) * t22 + t19 * t21 + 0.2e1 * (pkin(3) * mrSges(5,2) - t17) * t10 + (-pkin(3) * mrSges(5,1) + Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,2)) * t14) * t10) * qJD(3); m(5) * (t1 * t10 + t2 * t11 + (t10 * t3 + t11 * t4) * qJD(3)); 0; -t1 * mrSges(5,2) + (mrSges(5,1) + t19) * t2 + ((mrSges(4,2) * t6 - Ifges(4,6) - Ifges(5,6)) * t10 + (-mrSges(4,1) * t6 - mrSges(5,3) * pkin(3) + Ifges(4,5) + Ifges(5,5)) * t11) * qJD(3); (-mrSges(4,2) * t11 + (-mrSges(4,1) - t19) * t10) * qJD(3) - t18; 0; t15 * t19 + t18; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t5(1), t5(2), t5(4), t5(7); t5(2), t5(3), t5(5), t5(8); t5(4), t5(5), t5(6), t5(9); t5(7), t5(8), t5(9), t5(10);];
Mq = res;
