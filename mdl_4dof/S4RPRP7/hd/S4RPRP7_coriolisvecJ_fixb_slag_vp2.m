% Calculate vector of centrifugal and Coriolis load on the joints for
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:03
% EndTime: 2019-12-31 16:47:04
% DurationCPUTime: 0.54s
% Computational Cost: add. (281->112), mult. (655->157), div. (0->0), fcn. (216->2), ass. (0->42)
t49 = -mrSges(4,1) - mrSges(5,1);
t21 = sin(qJ(3));
t33 = qJD(3) * t21;
t48 = (((m(3) + m(4)) * qJ(2) + mrSges(3,3)) * qJD(1));
t35 = qJD(1) * t21;
t18 = -mrSges(5,2) * t35 + qJD(3) * mrSges(5,3);
t47 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t35 + t18;
t22 = cos(qJ(3));
t34 = qJD(1) * t22;
t46 = (mrSges(5,2) + mrSges(4,3)) * t34 + t49 * qJD(3);
t23 = -pkin(1) - pkin(5);
t41 = Ifges(5,1) * t22;
t40 = Ifges(4,4) * t21;
t39 = Ifges(4,4) * t22;
t38 = Ifges(5,5) * t21;
t37 = t22 * mrSges(4,2);
t19 = t23 * qJD(1) + qJD(2);
t36 = t22 * t19;
t32 = qJD(1) * qJD(2);
t31 = -Ifges(5,4) / 0.2e1 - Ifges(4,5) / 0.2e1;
t30 = 0.3e1 / 0.2e1 * Ifges(5,5) - 0.3e1 / 0.2e1 * Ifges(4,4);
t29 = Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1;
t4 = -qJD(3) * pkin(3) + qJD(4) - t36;
t28 = -t4 + t36;
t10 = qJD(3) * qJ(4) + t21 * t19;
t27 = t10 * t22 + t21 * t4;
t26 = pkin(3) * t22 + qJ(4) * t21;
t14 = t21 * pkin(3) - t22 * qJ(4) + qJ(2);
t25 = t46 * t21 + t47 * t22;
t2 = t26 * qJD(3) - t22 * qJD(4) + qJD(2);
t20 = Ifges(5,5) * t34;
t13 = qJD(1) * (t21 * mrSges(4,1) + t37);
t12 = (t21 * mrSges(5,1) - t22 * mrSges(5,3)) * qJD(1);
t11 = t26 * qJD(1);
t9 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t22 - t40) * qJD(1);
t8 = Ifges(5,4) * qJD(3) + (t38 + t41) * qJD(1);
t7 = Ifges(4,6) * qJD(3) + (-Ifges(4,2) * t21 + t39) * qJD(1);
t6 = Ifges(5,6) * qJD(3) + Ifges(5,3) * t35 + t20;
t5 = t14 * qJD(1);
t3 = (qJD(4) + t36) * qJD(3);
t1 = t2 * qJD(1);
t15 = [m(5) * (t1 * t14 + t5 * t2) + t2 * t12 + (t13 + (2 * t48)) * qJD(2) + (mrSges(4,2) * t32 - t1 * mrSges(5,3) + (t5 * mrSges(5,1) + t6 / 0.2e1 - t7 / 0.2e1 - t10 * mrSges(5,2) + t29 * qJD(3) + (m(5) * t10 + t47) * t23 + (0.2e1 * qJ(2) * mrSges(4,1) + t14 * mrSges(5,1) + t30 * t22) * qJD(1)) * qJD(3)) * t22 + (mrSges(4,1) * t32 + t1 * mrSges(5,1) + (m(5) * t23 - mrSges(5,2)) * t3 + (-t8 / 0.2e1 - t9 / 0.2e1 + t5 * mrSges(5,3) + t31 * qJD(3) + t28 * mrSges(5,2) + (-m(5) * t28 + t46) * t23 + (-0.2e1 * qJ(2) * mrSges(4,2) + t14 * mrSges(5,3) - t30 * t21 + (-0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,2)) * t22) * qJD(1)) * qJD(3)) * t21; m(5) * t3 * t21 + (m(5) * (-t21 * t36 + t27) + t25) * qJD(3) + (-m(5) * t5 - t12 - t13 - t48) * qJD(1); t3 * mrSges(5,3) + qJD(4) * t18 - t11 * t12 + ((t49 * t21 - t37) * qJD(3) - t25) * t19 + (t22 * t7 / 0.2e1 - t5 * (mrSges(5,1) * t22 + mrSges(5,3) * t21) + (-t21 * (Ifges(5,3) * t22 - t38) / 0.2e1 - qJ(2) * (mrSges(4,1) * t22 - mrSges(4,2) * t21)) * qJD(1) + t27 * mrSges(5,2) + ((-qJ(4) * mrSges(5,2) + t29) * t22 + (pkin(3) * mrSges(5,2) + t31) * t21) * qJD(3) - (t6 + t20 + (-Ifges(4,1) * t21 - t39) * qJD(1)) * t22 / 0.2e1 + (t8 + t9 + (-Ifges(4,2) * t22 - t40 + t41) * qJD(1)) * t21 / 0.2e1) * qJD(1) + (t3 * qJ(4) + t10 * qJD(4) - t5 * t11 + (-pkin(3) * t33 - t27) * t19) * m(5); -qJD(3) * t18 + (-mrSges(5,2) * t33 + t22 * t12) * qJD(1) + (-t10 * qJD(3) + t19 * t33 + t5 * t34) * m(5);];
tauc = t15(:);
