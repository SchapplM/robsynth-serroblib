% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:11
% EndTime: 2019-12-31 17:00:12
% DurationCPUTime: 0.22s
% Computational Cost: add. (291->64), mult. (518->67), div. (0->0), fcn. (340->2), ass. (0->29)
t20 = sin(qJ(2));
t21 = cos(qJ(2));
t19 = pkin(2) + qJ(4);
t29 = t20 * qJ(3);
t6 = -t19 * t21 - pkin(1) - t29;
t24 = m(5) * t6 - t20 * mrSges(5,2) - t21 * mrSges(5,3);
t22 = -t21 * pkin(2) - t29;
t8 = -pkin(1) + t22;
t9 = t21 * mrSges(4,2) - t20 * mrSges(4,3);
t36 = m(4) * t8 + t24 + t9;
t35 = Ifges(5,6) - Ifges(3,4) - Ifges(4,6);
t34 = pkin(3) + pkin(5);
t33 = -mrSges(4,1) - mrSges(5,1);
t32 = m(5) * qJD(2);
t1 = t36 * (t20 * pkin(2) - t21 * qJ(3)) + (-pkin(1) * mrSges(3,2) - t6 * mrSges(5,2) - t8 * mrSges(4,3) - t35 * t21) * t21 + (t24 * qJ(4) - pkin(1) * mrSges(3,1) - t8 * mrSges(4,2) + t6 * mrSges(5,3) + (Ifges(3,1) - Ifges(3,2) + Ifges(4,2) - Ifges(5,2) - Ifges(4,3) + Ifges(5,3)) * t21 + t35 * t20) * t20;
t31 = t1 * qJD(1);
t2 = t36 * t20;
t30 = t2 * qJD(1);
t3 = t24 * t21;
t28 = t3 * qJD(1);
t14 = m(5) * t19 + mrSges(5,3);
t27 = t14 * qJD(2);
t15 = mrSges(5,2) + mrSges(4,3) + (m(4) + m(5)) * qJ(3);
t26 = t15 * qJD(2);
t13 = t34 * t21;
t12 = t34 * t20;
t5 = -m(5) * t12 - t20 * mrSges(5,1);
t4 = m(5) * t13 + (m(4) * pkin(5) - t33) * t21;
t7 = [t1 * qJD(2) - t2 * qJD(3) - t3 * qJD(4), t4 * qJD(3) + t5 * qJD(4) + t31 + (m(5) * (-qJ(3) * t12 - t19 * t13) - t12 * mrSges(5,2) - t13 * mrSges(5,3) + (-pkin(2) * mrSges(4,1) - t19 * mrSges(5,1) - Ifges(4,4) + Ifges(3,5) + Ifges(5,5)) * t21 + (t33 * qJ(3) + Ifges(5,4) + Ifges(4,5) - Ifges(3,6)) * t20 + (m(4) * t22 - t21 * mrSges(3,1) + t20 * mrSges(3,2) + t9) * pkin(5)) * qJD(2), t4 * qJD(2) - t30, t5 * qJD(2) - t28; -t31, t15 * qJD(3) + t14 * qJD(4), t26, t27; t30, -m(5) * qJD(4) - t26, 0, -t32; t28, m(5) * qJD(3) - t27, t32, 0;];
Cq = t7;
