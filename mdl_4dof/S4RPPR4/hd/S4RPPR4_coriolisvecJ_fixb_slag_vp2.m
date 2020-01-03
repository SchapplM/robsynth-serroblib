% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:49
% EndTime: 2019-12-31 16:38:49
% DurationCPUTime: 0.26s
% Computational Cost: add. (215->67), mult. (473->108), div. (0->0), fcn. (181->4), ass. (0->30)
t31 = m(4) + m(5);
t16 = sin(qJ(4));
t30 = t16 / 0.2e1;
t17 = cos(qJ(4));
t29 = -t17 / 0.2e1;
t28 = Ifges(5,4) * t16;
t27 = Ifges(5,4) * t17;
t26 = (mrSges(4,3) * qJD(1));
t25 = Ifges(5,5) * qJD(4);
t24 = Ifges(5,6) * qJD(4);
t23 = qJD(1) * mrSges(5,3);
t13 = sin(pkin(6)) * pkin(1) + qJ(3);
t9 = qJD(1) * t13;
t22 = qJD(1) * t17;
t21 = qJD(1) * qJD(3);
t12 = -cos(pkin(6)) * pkin(1) - pkin(2) - pkin(5);
t20 = m(5) * t12 - mrSges(5,3);
t19 = 0.2e1 * t9;
t7 = t12 * qJD(1) + qJD(3);
t3 = -t16 * qJD(2) + t17 * t7;
t4 = t17 * qJD(2) + t16 * t7;
t18 = -t16 * t3 + t17 * t4;
t11 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t22;
t10 = -qJD(4) * mrSges(5,2) - t16 * t23;
t8 = qJD(1) * (t16 * mrSges(5,1) + t17 * mrSges(5,2));
t6 = t25 + (Ifges(5,1) * t17 - t28) * qJD(1);
t5 = t24 + (-Ifges(5,2) * t16 + t27) * qJD(1);
t2 = qJD(4) * t4;
t1 = qJD(4) * t3;
t14 = [(t31 * t19 + (2 * t26) + t8) * qJD(3) + (mrSges(5,2) * t21 - t20 * t2 + (-t4 * mrSges(5,3) - t5 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,4) * t22 - t24 / 0.2e1 + (m(5) * t4 + t10) * t12 + t19 * mrSges(5,1)) * qJD(4)) * t17 + (mrSges(5,1) * t21 + t20 * t1 + (t3 * mrSges(5,3) - t6 / 0.2e1 - t9 * mrSges(5,2) - t25 / 0.2e1 + (-m(5) * t3 - t11) * t12 + (-t13 * mrSges(5,2) + 0.3e1 / 0.2e1 * t28 + (-0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,2)) * t17) * qJD(1)) * qJD(4)) * t16; m(5) * (t1 * t17 + t2 * t16) + (m(5) * (-t16 * t4 - t17 * t3) - t16 * t10 - t17 * t11 + (-t16 ^ 2 - t17 ^ 2) * t23) * qJD(4); m(5) * (t1 * t16 - t2 * t17) + (-t31 * t9 - t26 - t8) * qJD(1) + (m(5) * t18 + t17 * t10 - t16 * t11) * qJD(4); -t2 * mrSges(5,1) - t1 * mrSges(5,2) - t3 * t10 + t4 * t11 + (-t9 * (mrSges(5,1) * t17 - mrSges(5,2) * t16) + t6 * t30 + t17 * t5 / 0.2e1 + ((-Ifges(5,1) * t16 - t27) * t29 + (-Ifges(5,2) * t17 - t28) * t30) * qJD(1) + t18 * mrSges(5,3) + (-Ifges(5,5) * t16 / 0.2e1 + Ifges(5,6) * t29) * qJD(4)) * qJD(1);];
tauc = t14(:);
