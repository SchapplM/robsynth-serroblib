% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPPRR1
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPPRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:50
% EndTime: 2019-12-05 14:57:53
% DurationCPUTime: 0.48s
% Computational Cost: add. (537->107), mult. (1418->177), div. (0->0), fcn. (1080->8), ass. (0->58)
t63 = -Ifges(6,1) / 0.2e1;
t40 = cos(qJ(5));
t49 = qJD(4) * t40;
t33 = Ifges(6,4) * t49;
t62 = -t33 / 0.2e1;
t34 = sin(pkin(9));
t36 = cos(pkin(9));
t39 = sin(qJ(4));
t41 = cos(qJ(4));
t61 = -t39 * t34 + t41 * t36;
t60 = 2 * m(6);
t38 = sin(qJ(5));
t59 = Ifges(6,4) * t38;
t35 = sin(pkin(8));
t52 = qJD(1) * t35;
t24 = t36 * qJD(2) - t34 * t52;
t25 = t34 * qJD(2) + t36 * t52;
t12 = t39 * t24 + t41 * t25;
t10 = t12 * qJD(4);
t27 = t41 * t34 + t39 * t36;
t17 = t27 * t35;
t58 = t10 * t17;
t57 = t10 * t61;
t54 = Ifges(6,5) * qJD(5);
t53 = Ifges(6,6) * qJD(5);
t51 = qJD(4) * t35;
t50 = qJD(4) * t38;
t11 = t41 * t24 - t39 * t25;
t9 = t11 * qJD(4);
t48 = t54 / 0.2e1;
t47 = -t53 / 0.2e1;
t28 = (-t40 * mrSges(6,1) + t38 * mrSges(6,2)) * qJD(4);
t46 = -qJD(4) * mrSges(5,1) + t28;
t37 = cos(pkin(8));
t32 = -t37 * qJD(1) + qJD(3);
t8 = qJD(4) * pkin(6) + t12;
t5 = t40 * t32 - t38 * t8;
t6 = t38 * t32 + t40 * t8;
t45 = -t38 * t5 + t40 * t6;
t18 = t61 * t35;
t14 = t40 * t18 - t37 * t38;
t13 = -t38 * t18 - t37 * t40;
t30 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t50;
t31 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t49;
t44 = -t38 * t30 + t40 * t31;
t7 = -qJD(4) * pkin(4) - t11;
t43 = t5 * mrSges(6,3) + t50 * t63 + t62 - t54 / 0.2e1 - t7 * mrSges(6,2);
t42 = t6 * mrSges(6,3) + t53 / 0.2e1 + (Ifges(6,2) * t40 + t59) * qJD(4) / 0.2e1 - t7 * mrSges(6,1);
t23 = (mrSges(6,1) * t38 + mrSges(6,2) * t40) * qJD(5) * qJD(4);
t20 = t27 * qJD(4);
t19 = t61 * qJD(4);
t16 = t61 * t51;
t15 = t27 * t51;
t4 = -t14 * qJD(5) + t38 * t15;
t3 = t13 * qJD(5) - t40 * t15;
t2 = -qJD(5) * t6 - t38 * t9;
t1 = qJD(5) * t5 + t40 * t9;
t21 = [t16 * t28 + t17 * t23 + t3 * t31 + t4 * t30 + m(5) * (-t11 * t16 - t12 * t15 + t9 * t18 + t58) + m(6) * (t1 * t14 + t2 * t13 + t7 * t16 + t6 * t3 + t5 * t4 + t58) + (-t16 * mrSges(5,1) + t15 * mrSges(5,2) + (-t13 * t40 - t14 * t38) * qJD(5) * mrSges(6,3)) * qJD(4); -t61 * t23 + t46 * t20 + (-t30 * t40 - t31 * t38) * t27 * qJD(5) + (-qJD(4) * mrSges(5,2) + t44) * t19 + m(5) * (-t11 * t20 + t12 * t19 + t9 * t27 - t57) + m(6) * (-t57 + t7 * t20 + t45 * t19 + (t1 * t40 - t2 * t38 + (-t6 * t38 - t5 * t40) * qJD(5)) * t27); m(6) * (t1 * t38 + t2 * t40) + (m(6) * t45 + (-t38 ^ 2 - t40 ^ 2) * qJD(4) * mrSges(6,3) + t44) * qJD(5); -t10 * mrSges(5,1) + (-m(6) * t10 - t23) * pkin(4) + (-m(6) * t7 - t46) * t12 + (-t10 * mrSges(6,1) + t1 * mrSges(6,3) - t11 * t31 + (-t11 * t6 / 0.2e1 + pkin(6) * t1 / 0.2e1) * t60 + (0.3e1 / 0.2e1 * t33 + t48 + (-m(6) * t5 - t30) * pkin(6) - t43) * qJD(5)) * t40 + (t10 * mrSges(6,2) - t2 * mrSges(6,3) + t11 * t30 + (t11 * t5 / 0.2e1 - pkin(6) * t2 / 0.2e1) * t60 + (t47 + (-m(6) * t6 - t31) * pkin(6) + (-0.3e1 / 0.2e1 * t59 + (-0.3e1 / 0.2e1 * Ifges(6,2) + 0.3e1 / 0.2e1 * Ifges(6,1)) * t40) * qJD(4) - t42) * qJD(5)) * t38; t2 * mrSges(6,1) - t1 * mrSges(6,2) + t6 * t30 - t5 * t31 + ((t48 + t62 + t43) * t40 + (t47 + (t59 / 0.2e1 + (t63 + Ifges(6,2) / 0.2e1) * t40) * qJD(4) + t42) * t38) * qJD(4);];
tauc = t21(:);
