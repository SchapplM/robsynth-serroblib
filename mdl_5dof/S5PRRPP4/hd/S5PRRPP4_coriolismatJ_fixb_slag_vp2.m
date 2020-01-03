% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:48
% EndTime: 2019-12-31 17:40:49
% DurationCPUTime: 0.25s
% Computational Cost: add. (354->83), mult. (698->92), div. (0->0), fcn. (429->2), ass. (0->39)
t50 = pkin(6) - qJ(5);
t49 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t48 = m(5) + m(6);
t47 = pkin(3) + pkin(4);
t29 = sin(qJ(3));
t30 = cos(qJ(3));
t38 = t30 * qJ(4);
t19 = -t29 * pkin(3) + t38;
t13 = -pkin(4) * t29 + t19;
t46 = m(6) * t13;
t15 = -t29 * t47 + t38;
t45 = m(6) * t15;
t44 = m(6) * t29;
t43 = -mrSges(5,2) + mrSges(6,3);
t40 = t29 * qJ(4);
t12 = t30 * t47 + pkin(2) + t40;
t32 = -pkin(3) * t30 - t40;
t16 = -pkin(2) + t32;
t24 = t30 * mrSges(6,2);
t18 = -t30 * mrSges(5,1) - t29 * mrSges(5,3);
t33 = m(5) * t16 + t18;
t1 = -t33 * t19 + (t24 + t46) * t12 + (t13 * mrSges(6,1) - pkin(2) * mrSges(4,2) - t16 * mrSges(5,3) + t49 * t30) * t30 + (-pkin(2) * mrSges(4,1) + t16 * mrSges(5,1) - t12 * mrSges(6,1) + t13 * mrSges(6,2) + (Ifges(4,1) + Ifges(5,1) + Ifges(6,1) - Ifges(4,2) - Ifges(6,2) - Ifges(5,3)) * t30 - t49 * t29) * t29;
t42 = t1 * qJD(2);
t2 = -t12 * t44 + (-t30 * mrSges(6,1) - t29 * mrSges(6,2) + t33) * t29;
t41 = t2 * qJD(2);
t3 = t29 * mrSges(6,1) - t24 + 0.2e1 * (-t15 / 0.4e1 - t13 / 0.4e1) * m(6);
t39 = t3 * qJD(2);
t17 = t50 * t29;
t20 = t50 * t30;
t4 = m(6) * (-t17 * t29 - t20 * t30) + (t29 ^ 2 + t30 ^ 2) * mrSges(6,3);
t37 = t4 * qJD(2);
t36 = qJD(3) * t29;
t22 = qJ(4) * t48 + mrSges(6,2) + mrSges(5,3);
t35 = t22 * qJD(3);
t34 = qJD(2) * t44;
t21 = t48 * t29;
t7 = m(6) * t20 + (m(5) * pkin(6) - t43) * t30;
t5 = -t45 / 0.2e1 + t46 / 0.2e1;
t6 = [0, 0, t21 * qJD(4) + (-mrSges(4,1) - mrSges(5,1) - mrSges(6,1)) * t36 + (m(5) * t19 + t24 + t45 + (-mrSges(4,2) + mrSges(5,3)) * t30) * qJD(3), t21 * qJD(3), 0; 0, qJD(3) * t1 - qJD(4) * t2 + qJD(5) * t4, t42 + t7 * qJD(4) + t5 * qJD(5) + (qJ(4) * t43 - Ifges(4,6) + Ifges(5,6) - Ifges(6,6)) * t36 + (-t20 * mrSges(6,1) - t17 * mrSges(6,2) + m(6) * (-qJ(4) * t17 - t20 * t47) + (-pkin(3) * mrSges(5,2) + mrSges(6,3) * t47 + Ifges(5,4) + Ifges(4,5) - Ifges(6,5)) * t30 + (m(5) * t32 - t30 * mrSges(4,1) + t29 * mrSges(4,2) + t18) * pkin(6)) * qJD(3), qJD(3) * t7 - t41, qJD(3) * t5 + t37; 0, t3 * qJD(5) - t42, t22 * qJD(4), t35, t39; 0, -qJD(5) * t44 + t41, -t35, 0, -t34; 0, -qJD(3) * t3 + qJD(4) * t44 - t37, -t39, t34, 0;];
Cq = t6;
