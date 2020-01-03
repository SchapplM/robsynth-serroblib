% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRP3
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
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:44
% EndTime: 2019-12-31 16:26:45
% DurationCPUTime: 0.51s
% Computational Cost: add. (271->115), mult. (793->157), div. (0->0), fcn. (320->2), ass. (0->50)
t63 = Ifges(4,4) + Ifges(5,4);
t59 = qJD(3) / 0.2e1;
t34 = cos(qJ(3));
t47 = qJD(2) * t34;
t58 = -t63 * t47 / 0.2e1;
t57 = m(4) * pkin(5);
t56 = m(5) * pkin(3);
t27 = -t34 * pkin(3) - pkin(2);
t17 = t27 * qJD(2) + qJD(4);
t55 = m(5) * t17;
t54 = pkin(2) * mrSges(4,1);
t53 = pkin(2) * mrSges(4,2);
t52 = -qJ(4) - pkin(5);
t33 = sin(qJ(3));
t44 = qJD(3) * qJD(2);
t40 = t33 * t44;
t51 = t34 * mrSges(5,2) * t44 + mrSges(5,1) * t40;
t49 = qJD(3) * pkin(3);
t48 = qJD(2) * t33;
t46 = t33 * qJD(1);
t45 = qJD(1) * qJD(3);
t43 = 0.3e1 / 0.2e1 * Ifges(5,4) + 0.3e1 / 0.2e1 * Ifges(4,4);
t42 = Ifges(5,5) / 0.2e1 + Ifges(4,5) / 0.2e1;
t41 = mrSges(4,3) + t57;
t23 = t52 * t34;
t22 = t52 * t33;
t39 = qJD(3) * t52;
t30 = t34 * qJD(1);
t4 = qJD(2) * t22 + t30;
t3 = t4 + t49;
t5 = -qJD(2) * t23 + t46;
t38 = -t3 * t33 + t34 * t5;
t37 = (-Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * qJD(3);
t16 = pkin(5) * t47 + t46;
t7 = -t33 * qJD(4) + t34 * t39;
t6 = t34 * qJD(4) + t33 * t39;
t36 = -t17 * mrSges(5,1) + t16 * mrSges(4,3) + t5 * mrSges(5,3) + ((Ifges(4,2) + Ifges(5,2)) * t34 + t63 * t33) * qJD(2) / 0.2e1 + (Ifges(5,6) + Ifges(4,6)) * t59;
t15 = -pkin(5) * t48 + t30;
t35 = t17 * mrSges(5,2) - t15 * mrSges(4,3) - t3 * mrSges(5,3) - t58 + (Ifges(4,1) + Ifges(5,1)) * t48 / 0.2e1 + (Ifges(4,5) + Ifges(5,5)) * t59;
t26 = t34 * t45;
t21 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t47;
t20 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t47;
t19 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t48;
t18 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t48;
t14 = (-t34 * mrSges(5,1) + t33 * mrSges(5,2)) * qJD(2);
t13 = t16 * qJD(3);
t12 = -pkin(5) * t40 + t26;
t2 = t7 * qJD(2) - t33 * t45;
t1 = t6 * qJD(2) + t26;
t8 = [m(4) * (t12 * t33 - t13 * t34) + m(5) * (t1 * t33 + t2 * t34) + ((t20 + t21) * t34 + (-t18 - t19) * t33 + m(4) * (-t15 * t33 + t16 * t34) + m(5) * t38 + (mrSges(4,3) + mrSges(5,3)) * qJD(2) * (-t33 ^ 2 - t34 ^ 2)) * qJD(3); m(5) * (-t1 * t23 + t2 * t22 + t3 * t7 + t5 * t6) + t7 * t18 + t6 * t20 + t27 * t51 + (t1 * mrSges(5,3) + t41 * t12 + (t42 * qJD(3) + (-m(4) * t15 - t19) * pkin(5) + (-t22 * mrSges(5,3) + t43 * t34 - 0.2e1 * t53) * qJD(2) + t35) * qJD(3)) * t34 + (-t2 * mrSges(5,3) + t41 * t13 + (pkin(3) * t14 - pkin(5) * t21 + t37 - t16 * t57 + pkin(3) * t55 + (t27 * t56 - 0.2e1 * t54 + t23 * mrSges(5,3) + (pkin(3) * mrSges(5,2) - t43) * t33 + (-pkin(3) * mrSges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t34) * qJD(2) - t36) * qJD(3)) * t33; -t13 * mrSges(4,1) - t12 * mrSges(4,2) - t1 * mrSges(5,2) - t15 * t21 + t16 * t19 - t4 * t20 + (mrSges(5,1) + t56) * t2 + (-m(5) * (-t3 + t4) + t18) * t5 + ((qJD(2) * t53 + (-pkin(3) * mrSges(5,3) + t42) * qJD(3) - t35 + t58) * t34 + (t37 + (t54 + (Ifges(5,4) / 0.2e1 + Ifges(4,4) / 0.2e1) * t33) * qJD(2) + (-t14 - t55) * pkin(3) + (Ifges(5,2) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(5,1) / 0.2e1 - Ifges(4,1) / 0.2e1) * t47 + t36) * t33) * qJD(2); (t33 * t18 - t34 * t20 + (t33 * t49 - t38) * m(5)) * qJD(2) + t51;];
tauc = t8(:);
