% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRP6
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:59
% EndTime: 2019-12-31 16:46:00
% DurationCPUTime: 0.52s
% Computational Cost: add. (294->119), mult. (699->170), div. (0->0), fcn. (246->2), ass. (0->51)
t30 = -pkin(1) - pkin(5);
t23 = t30 * qJD(1) + qJD(2);
t59 = -qJ(4) * qJD(1) + t23;
t58 = (((m(3) + m(4)) * qJ(2) + mrSges(3,3)) * qJD(1));
t28 = sin(qJ(3));
t52 = Ifges(4,4) * t28;
t29 = cos(qJ(3));
t51 = Ifges(4,4) * t29;
t50 = Ifges(5,4) * t28;
t49 = Ifges(5,4) * t29;
t48 = qJ(4) - t30;
t47 = qJD(1) * t28;
t46 = qJD(1) * t29;
t45 = qJD(3) * t28;
t44 = qJD(3) * t29;
t43 = t28 * qJD(4);
t42 = t29 * qJD(4);
t40 = qJ(4) * qJD(3);
t39 = qJD(1) * qJD(2);
t38 = 0.3e1 / 0.2e1 * Ifges(5,4) + 0.3e1 / 0.2e1 * Ifges(4,4);
t37 = -Ifges(5,5) / 0.2e1 - Ifges(4,5) / 0.2e1;
t36 = -Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1;
t12 = (t28 * mrSges(5,1) + t29 * mrSges(5,2)) * qJD(1);
t26 = t28 * pkin(3) + qJ(2);
t16 = t26 * qJD(1) + qJD(4);
t35 = -m(5) * t16 - t12;
t18 = t48 * t29;
t7 = t59 * t29;
t3 = qJD(3) * pkin(3) + t7;
t6 = t59 * t28;
t34 = -t28 * t3 + t29 * t6;
t24 = pkin(3) * t44 + qJD(2);
t33 = t28 * mrSges(4,1) + t29 * mrSges(4,2);
t20 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t47;
t22 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t46;
t32 = -t29 * t20 + t28 * t22;
t25 = qJD(1) * mrSges(5,1) * t44;
t21 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t46;
t19 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t47;
t17 = t48 * t28;
t15 = t24 * qJD(1);
t13 = qJD(1) * t33;
t11 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t29 - t52) * qJD(1);
t10 = Ifges(5,5) * qJD(3) + (Ifges(5,1) * t29 - t50) * qJD(1);
t9 = Ifges(4,6) * qJD(3) + (-Ifges(4,2) * t28 + t51) * qJD(1);
t8 = Ifges(5,6) * qJD(3) + (-Ifges(5,2) * t28 + t49) * qJD(1);
t5 = -qJD(3) * t18 - t43;
t4 = t48 * t45 - t42;
t2 = t23 * t44 + (-t29 * t40 - t43) * qJD(1);
t1 = -t23 * t45 + (t28 * t40 - t42) * qJD(1);
t14 = [t26 * t25 + m(5) * (-t1 * t18 + t15 * t26 + t16 * t24 - t2 * t17 + t3 * t4 + t6 * t5) + t5 * t19 + t4 * t21 + t24 * t12 + (t13 + (2 * t58)) * qJD(2) + (mrSges(4,2) * t39 + t15 * mrSges(5,2) - t1 * mrSges(5,3) + (-t6 * mrSges(5,3) + t30 * t20 - t8 / 0.2e1 - t9 / 0.2e1 + t16 * mrSges(5,1) + t36 * qJD(3) + (0.2e1 * qJ(2) * mrSges(4,1) + t17 * mrSges(5,3) - t38 * t29) * qJD(1)) * qJD(3)) * t29 + (mrSges(4,1) * t39 + t15 * mrSges(5,1) - t2 * mrSges(5,3) + (-t10 / 0.2e1 - t11 / 0.2e1 - t16 * mrSges(5,2) + t3 * mrSges(5,3) - t30 * t22 + t37 * qJD(3) + (-t18 * mrSges(5,3) - 0.2e1 * qJ(2) * mrSges(4,2) - t26 * mrSges(5,2) + t38 * t28 + (-0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(4,2)) * t29) * qJD(1)) * qJD(3)) * t28; m(5) * (t1 * t29 + t2 * t28) + (m(5) * t34 + t29 * t19 - t28 * t21 - t32) * qJD(3) + (-t13 + t35 - t58) * qJD(1); -t2 * mrSges(5,2) - t7 * t19 + (m(5) * pkin(3) + mrSges(5,1)) * t1 + (-m(5) * (-t3 + t7) + t21) * t6 + (-t33 * qJD(3) + t32) * t23 + (-t16 * (mrSges(5,1) * t29 - mrSges(5,2) * t28) + (-qJ(2) * (mrSges(4,1) * t29 - mrSges(4,2) * t28) - (-t49 - t51 + (-Ifges(4,1) - Ifges(5,1)) * t28) * t29 / 0.2e1) * qJD(1) + t35 * t29 * pkin(3) + t34 * mrSges(5,3) + (t36 * t29 + (pkin(3) * mrSges(5,3) + t37) * t28) * qJD(3) + (t10 + t11 + (-t50 - t52 + (-Ifges(4,2) - Ifges(5,2)) * t29) * qJD(1)) * t28 / 0.2e1 + (t8 + t9) * t29 / 0.2e1) * qJD(1); t25 + m(5) * t15 + (-mrSges(5,2) * t45 + t28 * t19 + t29 * t21 - m(5) * (-t28 * t6 - t29 * t3)) * qJD(1);];
tauc = t14(:);
