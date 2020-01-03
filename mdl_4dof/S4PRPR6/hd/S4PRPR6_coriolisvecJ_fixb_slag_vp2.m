% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:20
% EndTime: 2019-12-31 16:24:21
% DurationCPUTime: 0.68s
% Computational Cost: add. (460->115), mult. (1306->183), div. (0->0), fcn. (829->6), ass. (0->62)
t48 = cos(qJ(2));
t61 = qJD(1) * t48;
t54 = qJD(3) - t61;
t43 = sin(pkin(7));
t64 = pkin(5) + qJ(3);
t36 = t64 * t43;
t44 = cos(pkin(7));
t37 = t64 * t44;
t45 = sin(qJ(4));
t47 = cos(qJ(4));
t13 = -t36 * t47 - t37 * t45;
t51 = t43 * t45 - t44 * t47;
t72 = t13 * qJD(4) - t54 * t51;
t14 = -t36 * t45 + t37 * t47;
t33 = t43 * t47 + t44 * t45;
t50 = t33 * t48;
t71 = qJD(1) * t50 - t33 * qJD(3) - t14 * qJD(4);
t62 = t43 ^ 2 + t44 ^ 2;
t70 = t62 * mrSges(4,3);
t29 = t51 * qJD(4);
t69 = m(4) / 0.2e1;
t67 = -t29 / 0.2e1;
t30 = t33 * qJD(4);
t66 = -t30 / 0.2e1;
t28 = t33 * qJD(2);
t65 = Ifges(5,4) * t28;
t27 = t51 * qJD(2);
t52 = -mrSges(4,1) * t44 + mrSges(4,2) * t43;
t63 = mrSges(5,1) * t27 + mrSges(5,2) * t28 + qJD(2) * t52;
t46 = sin(qJ(2));
t60 = t46 * qJD(1);
t40 = -pkin(3) * t44 - pkin(2);
t35 = (qJD(3) + t61) * qJD(2);
t59 = t62 * t35;
t39 = qJD(2) * qJ(3) + t60;
t58 = t62 * t39;
t57 = qJD(2) * t60;
t22 = qJD(2) * t29;
t23 = qJD(2) * t30;
t9 = t23 * mrSges(5,1) - t22 * mrSges(5,2);
t56 = pkin(5) * qJD(2) + t39;
t55 = t62 * qJD(3);
t53 = t48 * t58;
t20 = t56 * t43;
t21 = t56 * t44;
t5 = -t20 * t47 - t21 * t45;
t6 = -t20 * t45 + t21 * t47;
t49 = qJD(2) ^ 2;
t38 = -qJD(2) * pkin(2) + t54;
t31 = t40 * qJD(2) + t54;
t26 = Ifges(5,4) * t27;
t25 = t51 * t46;
t24 = t33 * t46;
t16 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t28;
t15 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t27;
t11 = Ifges(5,1) * t28 + Ifges(5,5) * qJD(4) - t26;
t10 = -Ifges(5,2) * t27 + Ifges(5,6) * qJD(4) + t65;
t8 = -qJD(2) * t50 + t46 * t29;
t7 = -t48 * t27 - t46 * t30;
t2 = -t6 * qJD(4) - t33 * t35;
t1 = t5 * qJD(4) - t51 * t35;
t3 = [t7 * t15 + t8 * t16 - t48 * t9 + t63 * t46 * qJD(2) + (-t22 * t24 + t23 * t25) * mrSges(5,3) + m(5) * (-t1 * t25 - t2 * t24 + t5 * t8 + t6 * t7) + m(4) * t46 * t59 + 0.2e1 * (t53 * t69 + ((t38 - t61) * t69 + m(5) * (t31 - t61) / 0.2e1) * t46) * qJD(2) + (-t46 * mrSges(3,1) + (-mrSges(3,2) + t70) * t48) * t49; t40 * t9 + t10 * t66 + t31 * (t30 * mrSges(5,1) - t29 * mrSges(5,2)) + t11 * t67 + qJD(4) * (-Ifges(5,5) * t29 - Ifges(5,6) * t30) / 0.2e1 + t71 * t16 + t72 * t15 + (t23 * t51 - t27 * t66) * Ifges(5,2) + (-t22 * t33 + t28 * t67) * Ifges(5,1) + ((mrSges(5,1) * t51 + mrSges(5,2) * t33 + t52) * qJD(2) - t63) * t60 + (-t1 * t51 + t13 * t22 - t14 * t23 - t2 * t33 + t29 * t5 - t30 * t6) * mrSges(5,3) + (t59 + (-t62 * t61 + t55) * qJD(2)) * mrSges(4,3) + (t22 * t51 - t23 * t33 - t27 * t67 + t28 * t66) * Ifges(5,4) + (t1 * t14 + t2 * t13 - t31 * t60 + t40 * t57 + t71 * t5 + t72 * t6) * m(5) + (-pkin(2) * t57 + qJ(3) * t59 + t39 * t55 - (t38 * t46 + t53) * qJD(1)) * m(4); -m(5) * (-t27 * t6 - t28 * t5) + t27 * t15 + t28 * t16 - t49 * t70 + (m(5) * t60 + (-t58 + t60) * m(4)) * qJD(2) + t9; -Ifges(5,5) * t22 - Ifges(5,6) * t23 - t1 * mrSges(5,2) + t2 * mrSges(5,1) - t31 * (mrSges(5,1) * t28 - mrSges(5,2) * t27) - t28 * (-Ifges(5,1) * t27 - t65) / 0.2e1 + t28 * t10 / 0.2e1 - qJD(4) * (-Ifges(5,5) * t27 - Ifges(5,6) * t28) / 0.2e1 - t5 * t15 + t6 * t16 + (-t27 * t5 + t28 * t6) * mrSges(5,3) + (-Ifges(5,2) * t28 + t11 - t26) * t27 / 0.2e1;];
tauc = t3(:);
