% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:46
% EndTime: 2019-12-31 16:27:47
% DurationCPUTime: 0.47s
% Computational Cost: add. (259->112), mult. (756->154), div. (0->0), fcn. (291->2), ass. (0->50)
t28 = cos(qJ(3));
t37 = qJD(2) * t28;
t60 = -t37 / 0.2e1;
t59 = mrSges(4,1) + mrSges(5,1);
t51 = mrSges(5,2) + mrSges(4,3);
t58 = -Ifges(4,1) / 0.2e1;
t57 = Ifges(5,6) / 0.2e1;
t56 = Ifges(4,4) * t60;
t27 = sin(qJ(3));
t16 = pkin(5) * t37 + t27 * qJD(1);
t10 = qJD(3) * qJ(4) + t16;
t55 = m(5) * t10;
t54 = mrSges(5,2) * t37;
t53 = -qJD(2) / 0.2e1;
t52 = -qJD(3) / 0.2e1;
t38 = qJD(2) * t27;
t35 = pkin(5) * t38;
t36 = t28 * qJD(1);
t15 = -t35 + t36;
t50 = -t15 + qJD(4);
t49 = 0.2e1 * pkin(5);
t48 = m(4) / 0.2e1;
t47 = m(5) / 0.2e1;
t46 = (pkin(2) * mrSges(4,1));
t45 = (pkin(2) * mrSges(4,2));
t44 = Ifges(4,4) * t27;
t43 = Ifges(5,5) * t28;
t12 = t16 * qJD(3);
t42 = t12 * t28;
t41 = t59 * qJD(3) - t51 * t38;
t21 = qJD(3) * mrSges(5,3) + t54;
t40 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t37 + t21;
t34 = Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t33 = 0.3e1 / 0.2e1 * Ifges(5,5) - 0.3e1 / 0.2e1 * Ifges(4,4);
t32 = t57 - Ifges(4,6) / 0.2e1;
t31 = pkin(3) * t27 - qJ(4) * t28;
t17 = -t28 * pkin(3) - t27 * qJ(4) - pkin(2);
t4 = -qJD(3) * pkin(3) + t50;
t9 = t17 * qJD(2);
t30 = t15 * mrSges(4,3) + t9 * mrSges(5,3) + (Ifges(5,1) * t27 - t43) * t53 + t38 * t58 + t56 - t4 * mrSges(5,2) + (Ifges(5,4) + Ifges(4,5)) * t52;
t23 = Ifges(5,5) * t38;
t29 = t9 * mrSges(5,1) + qJD(3) * t57 + Ifges(5,3) * t60 + t23 / 0.2e1 + Ifges(4,6) * t52 + (Ifges(4,2) * t28 + t44) * t53 - t10 * mrSges(5,2) - t16 * mrSges(4,3);
t3 = t31 * qJD(3) - t27 * qJD(4);
t22 = qJD(3) * t36;
t14 = t31 * qJD(2);
t13 = (-t28 * mrSges(5,1) - t27 * mrSges(5,3)) * qJD(2);
t11 = -qJD(3) * t35 + t22;
t2 = t22 + (qJD(4) - t35) * qJD(3);
t1 = t3 * qJD(2);
t5 = [m(4) * (t11 * t27 - t42) + m(5) * (t2 * t27 - t42) + (t40 * t28 - t41 * t27 + m(4) * (-t15 * t27 + t16 * t28) + m(5) * (t10 * t28 + t27 * t4) + t51 * qJD(2) * (-t27 ^ 2 - t28 ^ 2)) * qJD(3); m(5) * (t1 * t17 + t9 * t3) + t3 * t13 + (-t1 * mrSges(5,1) + t2 * mrSges(5,2) + t11 * mrSges(4,3) + (t11 * t48 + t2 * t47) * t49 + (t34 * qJD(3) + (-t17 * mrSges(5,3) - t33 * t28 - (2 * t45)) * qJD(2) + (-m(4) * t15 + m(5) * t4 - t41) * pkin(5) - t30) * qJD(3)) * t28 + (-t1 * mrSges(5,3) + ((t48 + t47) * t49 + t51) * t12 + (t32 * qJD(3) + (-m(4) * t16 - t40 - t55) * pkin(5) + (-(2 * t46) + t17 * mrSges(5,1) + t33 * t27 + (-0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t28) * qJD(2) + t29) * qJD(3)) * t27; -t11 * mrSges(4,2) + t2 * mrSges(5,3) + qJD(4) * t21 - t14 * t13 + t41 * t16 - t40 * t15 - t59 * t12 + ((t56 + (t45 + t43 / 0.2e1) * qJD(2) + (-pkin(3) * mrSges(5,2) + t34) * qJD(3) + t30) * t28 + (-t23 / 0.2e1 + (t46 + t44 / 0.2e1) * qJD(2) + (-Ifges(5,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + t58) * t37 + (-qJ(4) * mrSges(5,2) + t32) * qJD(3) - t29) * t27) * qJD(2) + (-t12 * pkin(3) + t2 * qJ(4) + t50 * t10 - t9 * t14 - t4 * t16) * m(5); m(5) * t12 + (-t21 + t54 - t55) * qJD(3) + (m(5) * t9 + t13) * t38;];
tauc = t5(:);
