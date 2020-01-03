% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:32
% EndTime: 2019-12-31 16:42:33
% DurationCPUTime: 0.49s
% Computational Cost: add. (341->119), mult. (916->160), div. (0->0), fcn. (390->4), ass. (0->51)
t68 = Ifges(4,4) + Ifges(5,4);
t64 = qJD(3) / 0.2e1;
t39 = cos(qJ(3));
t55 = qJD(1) * t39;
t63 = -t68 * t55 / 0.2e1;
t29 = sin(pkin(6)) * pkin(1) + pkin(5);
t21 = t29 * qJD(1);
t44 = qJ(4) * qJD(1) + t21;
t38 = sin(qJ(3));
t53 = t38 * qJD(2);
t5 = t39 * t44 + t53;
t62 = m(5) * pkin(3);
t61 = m(4) * t29;
t47 = -cos(pkin(6)) * pkin(1) - pkin(2);
t20 = -pkin(3) * t39 + t47;
t12 = qJD(1) * t20 + qJD(4);
t60 = m(5) * t12;
t59 = (mrSges(5,1) * t38 + mrSges(5,2) * t39) * qJD(3) * qJD(1);
t57 = qJ(4) + t29;
t56 = qJD(1) * t38;
t54 = qJD(3) * t38;
t52 = t38 * qJD(4);
t33 = t39 * qJD(2);
t51 = t39 * qJD(4);
t49 = 0.3e1 / 0.2e1 * Ifges(5,4) + 0.3e1 / 0.2e1 * Ifges(4,4);
t48 = Ifges(5,5) / 0.2e1 + Ifges(4,5) / 0.2e1;
t46 = mrSges(4,3) + t61;
t45 = qJD(3) * t57;
t8 = qJD(3) * t33 - t21 * t54;
t4 = -t38 * t44 + t33;
t3 = qJD(3) * pkin(3) + t4;
t43 = -t3 * t38 + t39 * t5;
t42 = (-Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * qJD(3);
t11 = t21 * t39 + t53;
t24 = t47 * qJD(1);
t41 = t24 * mrSges(4,2) + t12 * mrSges(5,2) - t3 * mrSges(5,3) - t63 + (Ifges(4,1) + Ifges(5,1)) * t56 / 0.2e1 + (Ifges(4,5) + Ifges(5,5)) * t64;
t40 = -t24 * mrSges(4,1) - t12 * mrSges(5,1) + t11 * mrSges(4,3) + t5 * mrSges(5,3) + ((Ifges(4,2) + Ifges(5,2)) * t39 + t68 * t38) * qJD(1) / 0.2e1 + (Ifges(5,6) + Ifges(4,6)) * t64;
t26 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t55;
t25 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t55;
t23 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t56;
t22 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t56;
t19 = (-t39 * mrSges(5,1) + t38 * mrSges(5,2)) * qJD(1);
t18 = t57 * t39;
t17 = t57 * t38;
t10 = -t21 * t38 + t33;
t9 = t11 * qJD(3);
t7 = -t39 * t45 - t52;
t6 = -t38 * t45 + t51;
t2 = -qJD(1) * t52 - qJD(3) * t5;
t1 = (-qJ(4) * t54 + t51) * qJD(1) + t8;
t13 = [m(5) * (t1 * t18 - t2 * t17 + t3 * t7 + t5 * t6) + t7 * t22 + t6 * t25 + t20 * t59 + (t1 * mrSges(5,3) + t46 * t8 + (-t29 * t23 + t48 * qJD(3) - t46 * t10 + (mrSges(4,2) * t47 + t17 * mrSges(5,3) + t39 * t49) * qJD(1) + t41) * qJD(3)) * t39 + (-t2 * mrSges(5,3) + t46 * t9 + (pkin(3) * t19 - t29 * t26 + t42 - t11 * t61 + pkin(3) * t60 + (t20 * t62 + t47 * mrSges(4,1) - t18 * mrSges(5,3) + (pkin(3) * mrSges(5,2) - t49) * t38 + (-0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(4,2) - pkin(3) * mrSges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t39) * qJD(1) - t40) * qJD(3)) * t38; m(4) * (t8 * t38 - t9 * t39) + m(5) * (t1 * t38 + t2 * t39) + ((t25 + t26) * t39 + (-t22 - t23) * t38 + m(4) * (-t10 * t38 + t11 * t39) + m(5) * t43 + (mrSges(4,3) + mrSges(5,3)) * qJD(1) * (-t38 ^ 2 - t39 ^ 2)) * qJD(3); -t9 * mrSges(4,1) - t8 * mrSges(4,2) - t1 * mrSges(5,2) - t10 * t26 + t11 * t23 - t4 * t25 + (mrSges(5,1) + t62) * t2 + (-m(5) * (-t3 + t4) + t22) * t5 + ((t10 * mrSges(4,3) + (-mrSges(5,3) * pkin(3) + t48) * qJD(3) - t41 + t63) * t39 + ((Ifges(5,4) / 0.2e1 + Ifges(4,4) / 0.2e1) * t56 + t42 + (-t19 - t60) * pkin(3) + (Ifges(5,2) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(5,1) / 0.2e1 - Ifges(4,1) / 0.2e1) * t55 + t40) * t38) * qJD(1); (t38 * t22 - t39 * t25 + (pkin(3) * t54 - t43) * m(5)) * qJD(1) + t59;];
tauc = t13(:);
