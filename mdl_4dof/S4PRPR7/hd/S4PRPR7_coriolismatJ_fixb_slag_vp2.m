% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRPR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:37
% EndTime: 2019-12-31 16:25:38
% DurationCPUTime: 0.13s
% Computational Cost: add. (110->39), mult. (295->56), div. (0->0), fcn. (174->4), ass. (0->24)
t18 = sin(qJ(2));
t29 = m(5) * t18;
t17 = sin(qJ(4));
t19 = cos(qJ(4));
t10 = t17 * mrSges(5,1) + t19 * mrSges(5,2);
t32 = mrSges(4,3) + t10;
t31 = m(5) / 0.2e1;
t20 = cos(qJ(2));
t15 = t17 ^ 2;
t16 = t19 ^ 2;
t26 = t15 + t16;
t30 = (0.1e1 - t26) * t20 * t29;
t1 = (qJ(3) * mrSges(5,1) - Ifges(5,4) * t19) * t19 + (-qJ(3) * mrSges(5,2) + Ifges(5,4) * t17 + (-Ifges(5,1) + Ifges(5,2)) * t19) * t17;
t25 = t1 * qJD(2);
t24 = t10 * qJD(4);
t23 = qJD(1) * t30;
t6 = (t15 / 0.2e1 + t16 / 0.2e1 - 0.1e1 / 0.2e1) * t29;
t8 = (m(4) + m(5)) * qJ(3) + t32;
t22 = t6 * qJD(1) - t8 * qJD(2);
t21 = -pkin(2) - pkin(5);
t14 = t20 * qJ(3);
t4 = t29 / 0.2e1 + (t26 * t31 + m(4)) * t18;
t3 = t18 * (t19 * mrSges(5,1) - t17 * mrSges(5,2));
t2 = [qJD(2) * t30, t4 * qJD(3) + t3 * qJD(4) + t23 + ((-t26 * mrSges(5,3) - mrSges(3,1) + mrSges(4,2)) * t18 + m(4) * (-t18 * pkin(2) + t14) + 0.2e1 * (t26 * t21 * t18 + t14) * t31 + (-mrSges(3,2) + t32) * t20) * qJD(2), t4 * qJD(2), t3 * qJD(2) + t20 * t24; -t6 * qJD(3) - t23, t8 * qJD(3) + t1 * qJD(4), -t22, t25 + ((-mrSges(5,2) * t21 - Ifges(5,6)) * t19 + (-mrSges(5,1) * t21 - Ifges(5,5)) * t17) * qJD(4); t6 * qJD(2), t22, 0, -t24; 0, -t25, 0, 0;];
Cq = t2;
