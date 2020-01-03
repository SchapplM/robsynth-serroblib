% Calculate time derivative of joint inertia matrix for
% S4RPRP4
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP4_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP4_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:38
% EndTime: 2019-12-31 16:43:38
% DurationCPUTime: 0.11s
% Computational Cost: add. (62->35), mult. (151->50), div. (0->0), fcn. (75->4), ass. (0->14)
t11 = -cos(pkin(6)) * pkin(1) - pkin(2);
t15 = 0.2e1 * t11;
t14 = Ifges(5,5) - Ifges(4,4);
t6 = sin(qJ(3));
t13 = qJD(3) * t6;
t7 = cos(qJ(3));
t12 = 0.2e1 * t7;
t3 = sin(pkin(6)) * pkin(1) + pkin(5);
t10 = (m(5) * t3 + mrSges(5,2)) * t7;
t9 = -t7 * mrSges(5,1) - t6 * mrSges(5,3);
t8 = -t7 * pkin(3) - t6 * qJ(4);
t2 = t11 + t8;
t1 = qJD(3) * t7 * qJ(4) - pkin(3) * t13 + t6 * qJD(4);
t4 = [0.2e1 * (-m(5) * t2 - t9) * t1 + ((mrSges(4,2) * t15 - 0.2e1 * t2 * mrSges(5,3) - t14 * t12) * t7 + (0.2e1 * t2 * mrSges(5,1) + mrSges(4,1) * t15 + 0.2e1 * t14 * t6 + (Ifges(5,1) - Ifges(5,3) + Ifges(4,1) - Ifges(4,2)) * t12) * t6) * qJD(3); 0; 0; qJD(4) * t10 + ((-pkin(3) * mrSges(5,2) + Ifges(5,4) + Ifges(4,5)) * t7 + (-qJ(4) * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * t6 + (m(5) * t8 - t7 * mrSges(4,1) + t6 * mrSges(4,2) + t9) * t3) * qJD(3); m(5) * t1 + ((-mrSges(4,2) + mrSges(5,3)) * t7 + (-mrSges(4,1) - mrSges(5,1)) * t6) * qJD(3); 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4); qJD(3) * t10; m(5) * t13; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t4(1), t4(2), t4(4), t4(7); t4(2), t4(3), t4(5), t4(8); t4(4), t4(5), t4(6), t4(9); t4(7), t4(8), t4(9), t4(10);];
Mq = res;
