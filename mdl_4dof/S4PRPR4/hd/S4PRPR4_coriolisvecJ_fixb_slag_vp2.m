% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRPR4
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:56
% EndTime: 2019-12-31 16:21:57
% DurationCPUTime: 0.19s
% Computational Cost: add. (172->53), mult. (393->86), div. (0->0), fcn. (140->2), ass. (0->27)
t33 = ((m(4) + m(5)) * qJ(3));
t11 = sin(qJ(4));
t12 = cos(qJ(4));
t26 = -pkin(2) - pkin(5);
t10 = t26 * qJD(2) + qJD(3);
t3 = -t11 * qJD(1) + t12 * t10;
t4 = t12 * qJD(1) + t11 * t10;
t20 = t11 * t3 - t12 * t4;
t24 = Ifges(5,4) * t12;
t25 = Ifges(5,4) * t11;
t28 = -t12 / 0.2e1;
t29 = -t11 / 0.2e1;
t32 = -t20 * mrSges(5,3) - (Ifges(5,6) * qJD(4) + (-Ifges(5,2) * t11 + t24) * qJD(2)) * t28 - (Ifges(5,5) * qJD(4) + (Ifges(5,1) * t12 - t25) * qJD(2)) * t29;
t31 = t11 ^ 2;
t30 = t12 ^ 2;
t22 = qJD(2) * mrSges(5,3);
t1 = qJD(4) * t3;
t2 = qJD(4) * t4;
t21 = -t1 * t11 + t2 * t12;
t19 = t11 * mrSges(5,1) + t12 * mrSges(5,2);
t17 = qJ(3) * (mrSges(5,1) * t12 - mrSges(5,2) * t11);
t16 = (Ifges(5,5) * t29 + Ifges(5,6) * t28) * qJD(4);
t8 = -qJD(4) * mrSges(5,2) - t11 * t22;
t9 = qJD(4) * mrSges(5,1) - t12 * t22;
t15 = -m(5) * t21 + (-m(5) * t20 - t11 * t9 + t12 * t8) * qJD(4);
t7 = qJD(2) * t19;
t5 = [m(5) * (t1 * t12 + t2 * t11) + (-t11 * t8 - t12 * t9 + m(5) * (-t11 * t4 - t12 * t3) + (-t30 - t31) * t22) * qJD(4); qJD(3) * t7 + t21 * mrSges(5,3) + (t16 - t32) * qJD(4) + t15 * t26 + (((2 * mrSges(4,3)) + t19 + (2 * t33)) * qJD(3) + (0.2e1 * t17 + (-0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,2)) * t11 * t12 + (0.3e1 / 0.2e1 * t31 - 0.3e1 / 0.2e1 * t30) * Ifges(5,4)) * qJD(4)) * qJD(2); t15 + (-t7 + (-mrSges(4,3) - t33) * qJD(2)) * qJD(2); -t2 * mrSges(5,1) - t1 * mrSges(5,2) - t3 * t8 + t4 * t9 + ((-t17 + (-Ifges(5,1) * t11 - t24) * t28 + t11 * (-Ifges(5,2) * t12 - t25) / 0.2e1) * qJD(2) + t16 + t32) * qJD(2);];
tauc = t5(:);
