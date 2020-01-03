% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:34
% EndTime: 2019-12-31 16:31:35
% DurationCPUTime: 0.32s
% Computational Cost: add. (338->80), mult. (718->132), div. (0->0), fcn. (299->4), ass. (0->42)
t28 = qJD(2) + qJD(3);
t30 = sin(qJ(3));
t51 = qJD(2) * pkin(2);
t23 = t28 * pkin(6) + t30 * t51;
t29 = sin(qJ(4));
t31 = cos(qJ(4));
t11 = t31 * qJD(1) - t29 * t23;
t32 = cos(qJ(3));
t45 = qJD(3) * t51;
t42 = t32 * t45;
t2 = t11 * qJD(4) + t31 * t42;
t12 = t29 * qJD(1) + t31 * t23;
t3 = -t12 * qJD(4) - t29 * t42;
t61 = t2 * t31 - t3 * t29;
t58 = Ifges(5,4) * t29;
t56 = t28 * t29;
t55 = t28 * t31;
t54 = t29 * t32;
t53 = t31 * t32;
t52 = t32 * mrSges(4,2);
t50 = Ifges(5,5) * qJD(4);
t49 = Ifges(5,6) * qJD(4);
t48 = qJD(4) * t29;
t47 = qJD(4) * t31;
t46 = -qJD(2) - t28;
t43 = t47 / 0.2e1;
t40 = -t31 * mrSges(5,1) + t29 * mrSges(5,2);
t39 = (Ifges(5,2) * t31 + t58) * t28;
t38 = (mrSges(5,1) * t29 + mrSges(5,2) * t31) * qJD(4);
t21 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t56;
t22 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t55;
t35 = m(5) * (-t11 * t29 + t12 * t31) + t31 * t22 - t29 * t21;
t34 = -t22 * t48 - t21 * t47 + m(5) * (-t11 * t47 - t12 * t48 + t61);
t13 = t39 + t49;
t25 = Ifges(5,4) * t55;
t14 = Ifges(5,1) * t56 + t25 + t50;
t24 = -t28 * pkin(3) - t32 * t51;
t33 = t30 * t40 * t45 + qJD(4) ^ 2 * (Ifges(5,5) * t31 - Ifges(5,6) * t29) / 0.2e1 + t14 * t43 + t24 * t38 - (t39 + t13) * t48 / 0.2e1 + ((-t11 * t31 - t12 * t29) * qJD(4) + t61) * mrSges(5,3) + ((Ifges(5,1) * t31 - t58) * t48 + (0.3e1 * Ifges(5,4) * t31 + (Ifges(5,1) - 0.2e1 * Ifges(5,2)) * t29) * t43) * t28;
t27 = -t32 * pkin(2) - pkin(3);
t18 = t40 * t28;
t15 = t28 * t38;
t1 = [m(5) * (t2 * t29 + t3 * t31) + ((-t29 ^ 2 - t31 ^ 2) * t28 * mrSges(5,3) + t35) * qJD(4); t33 + ((t18 + m(5) * (qJD(2) * t27 + t24) + t46 * mrSges(4,1)) * t30 + (t46 * mrSges(4,2) + t35) * t32) * qJD(3) * pkin(2) + t34 * (t30 * pkin(2) + pkin(6)) + t27 * t15; t33 + (-t30 * t18 - m(5) * (-t11 * t54 + t12 * t53 + t24 * t30) + t21 * t54 - t22 * t53 + (t30 * mrSges(4,1) + t52) * t28 + (-t52 + (-m(5) * pkin(3) - mrSges(4,1)) * t30) * qJD(3)) * t51 + t34 * pkin(6) - pkin(3) * t15; t3 * mrSges(5,1) - t2 * mrSges(5,2) - t11 * t22 + t12 * t21 + ((t50 / 0.2e1 - t24 * mrSges(5,2) - t14 / 0.2e1 - t25 / 0.2e1 + t11 * mrSges(5,3)) * t31 + (-t49 / 0.2e1 - t24 * mrSges(5,1) + t13 / 0.2e1 + t12 * mrSges(5,3) + (t58 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t31) * t28) * t29) * t28;];
tauc = t1(:);
