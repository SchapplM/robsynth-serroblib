% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPPRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:26
% EndTime: 2019-12-05 14:59:27
% DurationCPUTime: 0.47s
% Computational Cost: add. (560->119), mult. (1397->195), div. (0->0), fcn. (1026->8), ass. (0->65)
t68 = -Ifges(6,1) / 0.2e1;
t40 = cos(qJ(5));
t54 = qJD(4) * t40;
t33 = Ifges(6,4) * t54;
t67 = -t33 / 0.2e1;
t66 = 2 * m(6);
t38 = sin(qJ(5));
t65 = Ifges(6,4) * t38;
t34 = sin(pkin(9));
t36 = cos(pkin(9));
t35 = sin(pkin(8));
t57 = qJD(1) * t35;
t27 = t34 * qJD(2) + t36 * t57;
t37 = cos(pkin(8));
t32 = -t37 * qJD(1) + qJD(3);
t39 = sin(qJ(4));
t41 = cos(qJ(4));
t16 = t41 * t27 + t39 * t32;
t10 = qJD(4) * t16;
t60 = t35 * t36;
t19 = t37 * t41 + t39 * t60;
t64 = t10 * t19;
t63 = t10 * t39;
t62 = t34 * t35;
t61 = t34 * t41;
t59 = Ifges(6,5) * qJD(5);
t58 = Ifges(6,6) * qJD(5);
t15 = -t39 * t27 + t41 * t32;
t9 = qJD(4) * t15;
t56 = qJD(4) * t38;
t55 = qJD(4) * t39;
t53 = qJD(4) * t41;
t52 = qJD(5) * qJD(4);
t51 = t41 * t60;
t50 = t34 * t55;
t49 = t59 / 0.2e1;
t48 = -t58 / 0.2e1;
t25 = (mrSges(6,1) * t38 + mrSges(6,2) * t40) * t52;
t42 = qJD(4) ^ 2;
t47 = t42 * mrSges(5,2) + t25;
t26 = -t36 * qJD(2) + t34 * t57;
t8 = qJD(4) * pkin(6) + t16;
t5 = t40 * t26 - t38 * t8;
t6 = t38 * t26 + t40 * t8;
t28 = (-t40 * mrSges(6,1) + t38 * mrSges(6,2)) * qJD(4);
t45 = -t42 * mrSges(5,1) + qJD(4) * t28;
t20 = -t37 * t39 + t51;
t11 = -t20 * t38 + t40 * t62;
t12 = t20 * t40 + t38 * t62;
t22 = -t38 * t36 + t40 * t61;
t21 = -t40 * t36 - t38 * t61;
t7 = -qJD(4) * pkin(4) - t15;
t44 = t5 * mrSges(6,3) + t56 * t68 + t67 - t59 / 0.2e1 - t7 * mrSges(6,2);
t43 = t6 * mrSges(6,3) + t58 / 0.2e1 + (Ifges(6,2) * t40 + t65) * qJD(4) / 0.2e1 - t7 * mrSges(6,1);
t30 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t54;
t29 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t56;
t18 = qJD(4) * t51 - t37 * t55;
t17 = t19 * qJD(4);
t14 = -qJD(5) * t22 + t38 * t50;
t13 = qJD(5) * t21 - t40 * t50;
t4 = qJD(5) * t11 - t17 * t40;
t3 = -qJD(5) * t12 + t17 * t38;
t2 = -qJD(5) * t6 - t38 * t9;
t1 = qJD(5) * t5 + t40 * t9;
t23 = [t18 * t28 + t19 * t25 + t3 * t29 + t4 * t30 + m(5) * (-t15 * t18 - t16 * t17 + t9 * t20 + t64) + m(6) * (t1 * t12 + t2 * t11 + t7 * t18 + t5 * t3 + t6 * t4 + t64) + (-t18 * mrSges(5,1) + t17 * mrSges(5,2) + (-t11 * t40 - t12 * t38) * qJD(5) * mrSges(6,3)) * qJD(4); m(6) * (t1 * t22 + t6 * t13 + t5 * t14 + t2 * t21) + t13 * t30 + t14 * t29 + (-t21 * t40 - t22 * t38) * mrSges(6,3) * t52 + (t45 * t41 + t47 * t39 + m(5) * (-t15 * t53 - t16 * t55 + t41 * t9 + t63) + m(6) * (t53 * t7 + t63)) * t34; ((-t38 * t29 + t40 * t30) * qJD(4) + m(6) * (-t5 * t56 + t54 * t6 - t10) - t47) * t41 + (m(6) * (qJD(4) * t7 + t1 * t40 - t2 * t38) + t45 + (-t40 * t29 - t38 * t30 + m(6) * (-t38 * t6 - t40 * t5)) * qJD(5)) * t39; -t10 * mrSges(5,1) + (-m(6) * t10 - t25) * pkin(4) + (-m(6) * t7 + qJD(4) * mrSges(5,1) - t28) * t16 + (-t10 * mrSges(6,1) + t1 * mrSges(6,3) - t15 * t30 + (-t15 * t6 / 0.2e1 + pkin(6) * t1 / 0.2e1) * t66 + (0.3e1 / 0.2e1 * t33 + t49 + (-m(6) * t5 - t29) * pkin(6) - t44) * qJD(5)) * t40 + (t10 * mrSges(6,2) - t2 * mrSges(6,3) + t15 * t29 + (t15 * t5 / 0.2e1 - pkin(6) * t2 / 0.2e1) * t66 + (t48 + (-m(6) * t6 - t30) * pkin(6) + (-0.3e1 / 0.2e1 * t65 + (0.3e1 / 0.2e1 * Ifges(6,1) - 0.3e1 / 0.2e1 * Ifges(6,2)) * t40) * qJD(4) - t43) * qJD(5)) * t38; t2 * mrSges(6,1) - t1 * mrSges(6,2) + t6 * t29 - t5 * t30 + ((t49 + t67 + t44) * t40 + (t48 + (t65 / 0.2e1 + (t68 + Ifges(6,2) / 0.2e1) * t40) * qJD(4) + t43) * t38) * qJD(4);];
tauc = t23(:);
