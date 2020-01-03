% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRP5
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:51
% EndTime: 2019-12-31 16:28:53
% DurationCPUTime: 0.67s
% Computational Cost: add. (341->137), mult. (974->187), div. (0->0), fcn. (432->4), ass. (0->60)
t81 = Ifges(4,4) + Ifges(5,4);
t37 = sin(qJ(2));
t58 = t37 * qJD(1);
t27 = qJD(2) * pkin(5) + t58;
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t79 = t27 * (t36 ^ 2 + t38 ^ 2);
t78 = -qJD(3) / 0.2e1;
t54 = qJD(2) * qJD(3);
t15 = (mrSges(5,1) * t36 + mrSges(5,2) * t38) * t54;
t59 = qJD(3) * t36;
t17 = (pkin(3) * t59 + t58) * qJD(2);
t73 = m(5) * t17 + t15;
t60 = qJD(2) * t38;
t71 = -t81 * t60 / 0.2e1;
t31 = -t38 * pkin(3) - pkin(2);
t39 = cos(qJ(2));
t56 = t39 * qJD(1);
t10 = t31 * qJD(2) + qJD(4) - t56;
t18 = (-t38 * mrSges(5,1) + t36 * mrSges(5,2)) * qJD(2);
t28 = -qJD(2) * pkin(2) - t56;
t69 = m(5) * t10;
t46 = qJ(4) * qJD(2) + t27;
t7 = t46 * t38;
t70 = (t18 + t69) * pkin(3) + t28 * mrSges(4,1) + t10 * mrSges(5,1) - t7 * mrSges(5,3) - ((Ifges(4,2) + Ifges(5,2)) * t38 + t81 * t36) * qJD(2) / 0.2e1 + (Ifges(5,6) + Ifges(4,6)) * t78;
t66 = -qJ(4) - pkin(5);
t65 = t18 + (-t38 * mrSges(4,1) + t36 * mrSges(4,2)) * qJD(2);
t62 = qJD(2) * t36;
t20 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t62;
t21 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t62;
t64 = -t20 - t21;
t22 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t60;
t23 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t60;
t63 = t22 + t23;
t61 = qJD(2) * t37;
t57 = t38 * qJD(4);
t55 = qJ(4) * qJD(3);
t53 = qJD(3) * t38 * t27;
t52 = 0.3e1 / 0.2e1 * Ifges(5,4) + 0.3e1 / 0.2e1 * Ifges(4,4);
t51 = Ifges(5,5) / 0.2e1 + Ifges(4,5) / 0.2e1;
t50 = m(4) * pkin(5) + mrSges(4,3);
t48 = qJD(2) * t56;
t47 = qJD(3) * t66;
t3 = -t27 * t59 + t38 * t48;
t45 = (-Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * qJD(3);
t44 = t38 * t21 + t36 * t23;
t6 = t46 * t36;
t5 = qJD(3) * pkin(3) - t6;
t43 = m(5) * (-t36 * t5 + t38 * t7);
t42 = -t28 * mrSges(4,2) - t10 * mrSges(5,2) + t5 * mrSges(5,3) + t71 - (Ifges(4,1) + Ifges(5,1)) * t62 / 0.2e1 + (Ifges(4,5) + Ifges(5,5)) * t78;
t40 = qJD(2) ^ 2;
t25 = t66 * t38;
t24 = t66 * t36;
t16 = (mrSges(4,1) * t36 + mrSges(4,2) * t38) * t54;
t9 = -t36 * qJD(4) + t38 * t47;
t8 = t36 * t47 + t57;
t4 = -t36 * t48 - t53;
t2 = -t53 + (-t38 * t55 + (-qJD(4) - t56) * t36) * qJD(2);
t1 = (-t36 * t55 + t57) * qJD(2) + t3;
t11 = [(-t40 * mrSges(3,2) - t16 + (t79 * m(4) + t64 * t36 + t63 * t38 + t43) * qJD(2) - t73) * t39 + (-t40 * mrSges(3,1) + t65 * qJD(2) + m(4) * (qJD(2) * t28 + t3 * t38 - t36 * t4 - t48) + m(5) * (qJD(2) * t10 + t1 * t38 - t2 * t36) + (-t38 * t20 - t36 * t22 + m(5) * (-t36 * t7 - t38 * t5) - t44) * qJD(3)) * t37; m(5) * (-t1 * t25 + t17 * t31 + t2 * t24 + t5 * t9 + t7 * t8) - pkin(2) * t16 + t9 * t20 + t8 * t22 + t31 * t15 + ((-t65 - t69) * t37 + (-pkin(2) * t61 - t28 * t37 - t39 * t79) * m(4)) * qJD(1) + (t17 * mrSges(5,2) - t2 * mrSges(5,3) - t50 * t4 + (mrSges(4,2) * t61 + (m(5) * t5 - t64) * t39) * qJD(1) + (-pkin(5) * t23 + t45 + (t25 * mrSges(5,3) - t52 * t36) * qJD(2) + t70) * qJD(3)) * t36 + (-t17 * mrSges(5,1) + t1 * mrSges(5,3) + t50 * t3 + (-mrSges(4,1) * t61 + (-m(5) * t7 - t63) * t39) * qJD(1) + (-pkin(5) * t21 + t51 * qJD(3) + (-t24 * mrSges(5,3) + t52 * t38 + (0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t36) * qJD(2) - t42) * qJD(3)) * t38; t4 * mrSges(4,1) - t3 * mrSges(4,2) - t1 * mrSges(5,2) + t6 * t22 + t44 * t27 + (m(5) * pkin(3) + mrSges(5,1)) * t2 + (-m(5) * (-t5 - t6) + t20) * t7 + (((-pkin(3) * mrSges(5,3) + t51) * qJD(3) + t42 + t71) * t38 + ((Ifges(5,4) / 0.2e1 + Ifges(4,4) / 0.2e1) * t62 + t45 + (Ifges(5,2) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(5,1) / 0.2e1 - Ifges(4,1) / 0.2e1) * t60 - t70) * t36) * qJD(2); (t36 * t20 - t38 * t22 - t43) * qJD(2) + t73;];
tauc = t11(:);
