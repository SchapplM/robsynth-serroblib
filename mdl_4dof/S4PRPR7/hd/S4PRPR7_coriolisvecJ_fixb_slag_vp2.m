% Calculate vector of centrifugal and Coriolis load on the joints for
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:37
% EndTime: 2019-12-31 16:25:38
% DurationCPUTime: 0.32s
% Computational Cost: add. (193->85), mult. (484->126), div. (0->0), fcn. (192->4), ass. (0->43)
t15 = sin(qJ(4));
t10 = -t15 * qJD(2) * mrSges(5,3) - qJD(4) * mrSges(5,2);
t17 = cos(qJ(4));
t32 = qJD(2) * t17;
t11 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t32;
t22 = t17 * t10 - t15 * t11;
t46 = t22 * qJD(4);
t19 = -pkin(2) - pkin(5);
t18 = cos(qJ(2));
t28 = t18 * qJD(1);
t25 = qJD(3) - t28;
t6 = t19 * qJD(2) + t25;
t45 = (t15 ^ 2 + t17 ^ 2) * t6;
t44 = t15 / 0.2e1;
t43 = -t17 / 0.2e1;
t16 = sin(qJ(2));
t27 = qJD(1) * qJD(2);
t26 = t16 * t27;
t1 = -qJD(4) * t15 * t6 + t17 * t26;
t42 = t1 * t17;
t40 = Ifges(5,4) * t15;
t39 = Ifges(5,4) * t17;
t29 = t16 * qJD(1);
t12 = qJD(2) * qJ(3) + t29;
t38 = t12 * t18;
t37 = t17 * t11;
t36 = qJD(2) * pkin(2);
t34 = Ifges(5,5) * qJD(4);
t33 = Ifges(5,6) * qJD(4);
t31 = qJD(4) * t17;
t30 = t12 * qJD(2);
t2 = t15 * t26 + t6 * t31;
t24 = t2 * t15 + t42;
t23 = mrSges(5,1) * t17 - mrSges(5,2) * t15;
t8 = (qJD(3) + t28) * qJD(2);
t21 = t8 * qJ(3) + t12 * qJD(3);
t20 = qJD(2) ^ 2;
t9 = t25 - t36;
t7 = (t15 * mrSges(5,1) + t17 * mrSges(5,2)) * qJD(2);
t5 = t23 * qJD(4) * qJD(2);
t4 = t34 + (Ifges(5,1) * t17 - t40) * qJD(2);
t3 = t33 + (-Ifges(5,2) * t15 + t39) * qJD(2);
t13 = [(qJD(2) * t7 + (-mrSges(3,2) + mrSges(4,3)) * t20 - t46 + m(4) * t30 + m(5) * (-t24 + t30)) * t18 + (t5 + (-mrSges(3,1) + mrSges(4,2)) * t20 + (t15 * t10 + t37) * qJD(2) + m(4) * (qJD(2) * t9 - t18 * t27 + t8) + m(5) * (qJD(2) * t45 + t8)) * t16; qJ(3) * t5 + qJD(3) * t7 + (t8 * mrSges(5,2) - t1 * mrSges(5,3)) * t17 + (qJD(2) * qJD(3) + t8) * mrSges(4,3) + m(4) * t21 + m(5) * (t19 * t42 + t21) + (-t3 / 0.2e1 + t12 * mrSges(5,1) + t19 * t10 - t33 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,4) * t32) * t31 + (-t16 * t37 + (-qJD(2) * mrSges(4,3) - t7) * t18 - m(5) * (t16 * t45 + t38) + (-t38 + (-t36 - t9) * t16) * m(4)) * qJD(1) + (-t10 * t29 + t8 * mrSges(5,1) + (m(5) * t19 - mrSges(5,3)) * t2 + (-t4 / 0.2e1 - t12 * mrSges(5,2) - t19 * t11 - t34 / 0.2e1 + (0.3e1 / 0.2e1 * t40 + (-0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,2)) * t17) * qJD(2)) * qJD(4)) * t15; m(5) * t24 - t20 * mrSges(4,3) + t46 + (-m(5) * t12 - t7 + (-t12 + t29) * m(4)) * qJD(2); t1 * mrSges(5,1) - t2 * mrSges(5,2) - t22 * t6 + (-t12 * t23 + t4 * t44 + t17 * t3 / 0.2e1 + ((-Ifges(5,1) * t15 - t39) * t43 + (-Ifges(5,2) * t17 - t40) * t44) * qJD(2) + (-Ifges(5,5) * t15 / 0.2e1 + Ifges(5,6) * t43) * qJD(4)) * qJD(2);];
tauc = t13(:);
