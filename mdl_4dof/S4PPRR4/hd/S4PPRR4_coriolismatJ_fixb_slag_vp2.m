% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PPRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:28
% EndTime: 2019-12-31 16:18:28
% DurationCPUTime: 0.10s
% Computational Cost: add. (135->28), mult. (338->48), div. (0->0), fcn. (318->6), ass. (0->24)
t21 = cos(qJ(4));
t17 = t21 ^ 2;
t32 = cos(qJ(3));
t19 = sin(qJ(4));
t31 = t19 * mrSges(5,1);
t30 = t21 * mrSges(5,2);
t29 = -t19 ^ 2 - t17;
t18 = sin(pkin(7));
t20 = sin(qJ(3));
t27 = cos(pkin(7));
t12 = t20 * t18 - t32 * t27;
t13 = t32 * t18 + t20 * t27;
t1 = m(5) * (0.1e1 + t29) * t13 * t12;
t28 = t1 * qJD(1);
t26 = t29 * t12;
t25 = -t21 * mrSges(5,1) + t19 * mrSges(5,2);
t14 = t30 + t31;
t2 = -Ifges(5,4) * t17 + pkin(3) * t14 + (Ifges(5,4) * t19 + (-Ifges(5,1) + Ifges(5,2)) * t21) * t19;
t23 = t30 / 0.2e1 + t31 / 0.2e1;
t3 = (-t14 / 0.2e1 + t23) * t12;
t24 = t3 * qJD(1) + t2 * qJD(3);
t22 = t13 * t25;
t4 = (t14 / 0.2e1 + t23) * t12;
t5 = [t1 * qJD(3), 0, t28 + (t12 * mrSges(4,2) - t13 * mrSges(4,1) + t22 + m(5) * (-pkin(3) * t13 + pkin(5) * t26) + mrSges(5,3) * t26) * qJD(3) + t4 * qJD(4), t4 * qJD(3) + qJD(4) * t22; 0, 0, 0, -t14 * qJD(4); -t3 * qJD(4) - t28, 0, -t2 * qJD(4), (Ifges(5,5) * t21 - Ifges(5,6) * t19 + t25 * pkin(5)) * qJD(4) - t24; t3 * qJD(3), 0, t24, 0;];
Cq = t5;
