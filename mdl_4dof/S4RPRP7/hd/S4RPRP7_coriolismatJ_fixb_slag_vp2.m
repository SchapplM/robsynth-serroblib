% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRP7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP7_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP7_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:03
% EndTime: 2019-12-31 16:47:04
% DurationCPUTime: 0.12s
% Computational Cost: add. (158->41), mult. (276->45), div. (0->0), fcn. (174->2), ass. (0->21)
t16 = sin(qJ(3));
t17 = cos(qJ(3));
t31 = t16 * mrSges(5,1) - t17 * mrSges(5,3);
t30 = Ifges(5,5) - Ifges(4,4);
t29 = -t16 * mrSges(4,1) - t17 * mrSges(4,2);
t20 = -t16 * pkin(3) + t17 * qJ(4);
t28 = m(5) * t20 + t29 - t31;
t6 = qJ(2) - t20;
t21 = m(5) * t6 + t31;
t1 = (-qJ(2) * mrSges(4,2) + t6 * mrSges(5,3) + t21 * qJ(4) - t30 * t16) * t16 + (t21 * pkin(3) + qJ(2) * mrSges(4,1) + t6 * mrSges(5,1) + (-Ifges(4,1) - Ifges(5,1) + Ifges(4,2) + Ifges(5,3)) * t16 + t30 * t17) * t17;
t26 = t1 * qJD(1);
t2 = -mrSges(3,3) - t21 + (-m(3) - m(4)) * qJ(2) + t29;
t25 = t2 * qJD(1);
t3 = t21 * t17;
t24 = t3 * qJD(1);
t9 = m(5) * qJ(4) + mrSges(5,3);
t23 = t9 * qJD(3);
t22 = qJD(3) * t16;
t18 = -pkin(1) - pkin(5);
t5 = (m(5) * t18 - mrSges(5,2)) * t16;
t4 = [-t2 * qJD(2) + t1 * qJD(3) - t3 * qJD(4), -t25, t26 + t5 * qJD(4) + (pkin(3) * mrSges(5,2) - Ifges(5,4) - Ifges(4,5)) * t22 + ((-qJ(4) * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * t17 + t28 * t18) * qJD(3), t5 * qJD(3) - t24; t25, 0, m(5) * t16 * qJD(4) + t28 * qJD(3), m(5) * t22; -t26, 0, t9 * qJD(4), t23; t24, 0, -t23, 0;];
Cq = t4;
