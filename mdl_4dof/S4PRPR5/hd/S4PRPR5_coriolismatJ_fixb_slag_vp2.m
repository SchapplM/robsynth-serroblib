% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:59
% EndTime: 2019-12-31 16:22:59
% DurationCPUTime: 0.13s
% Computational Cost: add. (149->33), mult. (360->56), div. (0->0), fcn. (338->6), ass. (0->26)
t20 = cos(qJ(4));
t18 = t20 ^ 2;
t19 = sin(qJ(4));
t27 = -t19 ^ 2 - t18;
t24 = sin(pkin(7));
t25 = cos(pkin(7));
t30 = sin(qJ(2));
t31 = cos(qJ(2));
t8 = t24 * t30 - t25 * t31;
t34 = t27 * t8;
t23 = t20 * mrSges(5,1) - t19 * mrSges(5,2);
t9 = -t24 * t31 - t25 * t30;
t33 = t9 * t23;
t29 = t19 * mrSges(5,1);
t28 = t20 * mrSges(5,2);
t1 = m(5) * (-0.1e1 - t27) * t9 * t8;
t26 = t1 * qJD(1);
t10 = t28 + t29;
t15 = -t25 * pkin(2) - pkin(3);
t2 = t15 * t10 + t18 * Ifges(5,4) + (-Ifges(5,4) * t19 + (Ifges(5,1) - Ifges(5,2)) * t20) * t19;
t21 = t28 / 0.2e1 + t29 / 0.2e1;
t3 = (-t10 / 0.2e1 + t21) * t8;
t22 = t3 * qJD(1) - t2 * qJD(2);
t14 = t24 * pkin(2) + pkin(5);
t4 = (t10 / 0.2e1 + t21) * t8;
t5 = [t1 * qJD(2), t26 + (-t31 * mrSges(3,2) - t30 * mrSges(3,1) + t8 * mrSges(4,2) + t9 * mrSges(4,1) + m(4) * (-t24 * t8 + t25 * t9) * pkin(2) + t33 + m(5) * (t14 * t34 - t15 * t9) + mrSges(5,3) * t34) * qJD(2) + t4 * qJD(4), 0, t4 * qJD(2) + qJD(4) * t33; -t3 * qJD(4) - t26, t2 * qJD(4), 0, (Ifges(5,5) * t20 - Ifges(5,6) * t19 - t23 * t14) * qJD(4) - t22; 0, 0, 0, -t10 * qJD(4); t3 * qJD(2), t22, 0, 0;];
Cq = t5;
