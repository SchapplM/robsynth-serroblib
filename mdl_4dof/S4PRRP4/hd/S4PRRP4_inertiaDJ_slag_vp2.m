% Calculate time derivative of joint inertia matrix for
% S4PRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP4_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP4_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:46
% EndTime: 2019-12-31 16:27:46
% DurationCPUTime: 0.10s
% Computational Cost: add. (51->33), mult. (140->48), div. (0->0), fcn. (64->2), ass. (0->12)
t12 = -2 * pkin(2);
t11 = Ifges(5,5) - Ifges(4,4);
t4 = sin(qJ(3));
t10 = qJD(3) * t4;
t5 = cos(qJ(3));
t9 = 0.2e1 * t5;
t8 = (m(5) * pkin(5) + mrSges(5,2)) * t5;
t7 = -t5 * mrSges(5,1) - t4 * mrSges(5,3);
t6 = -t5 * pkin(3) - t4 * qJ(4);
t2 = -pkin(2) + t6;
t1 = qJD(3) * t5 * qJ(4) - pkin(3) * t10 + t4 * qJD(4);
t3 = [0; 0; 0.2e1 * (-m(5) * t2 - t7) * t1 + (((mrSges(4,2) * t12) - 0.2e1 * t2 * mrSges(5,3) - t11 * t9) * t5 + ((mrSges(4,1) * t12) + 0.2e1 * t2 * mrSges(5,1) + 0.2e1 * t11 * t4 + (Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,3)) * t9) * t4) * qJD(3); m(5) * t1 + ((-mrSges(4,2) + mrSges(5,3)) * t5 + (-mrSges(4,1) - mrSges(5,1)) * t4) * qJD(3); qJD(4) * t8 + ((-pkin(3) * mrSges(5,2) + Ifges(5,4) + Ifges(4,5)) * t5 + (-qJ(4) * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * t4 + (m(5) * t6 - t5 * mrSges(4,1) + t4 * mrSges(4,2) + t7) * pkin(5)) * qJD(3); 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4); m(5) * t10; qJD(3) * t8; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;
