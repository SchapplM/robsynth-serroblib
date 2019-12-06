% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:18
% EndTime: 2019-12-05 15:26:21
% DurationCPUTime: 0.61s
% Computational Cost: add. (500->115), mult. (1177->166), div. (0->0), fcn. (680->6), ass. (0->64)
t31 = sin(pkin(8));
t34 = sin(qJ(2));
t55 = qJD(1) * t34;
t27 = t31 * t55;
t32 = cos(pkin(8));
t36 = cos(qJ(2));
t53 = t36 * qJD(1);
t51 = t32 * t53;
t15 = -t27 + t51;
t52 = -t15 + qJD(4);
t30 = t31 * pkin(2) + qJ(4);
t26 = qJD(2) * pkin(2) + t53;
t10 = t31 * t26 + t32 * t55;
t7 = qJD(2) * qJ(4) + t10;
t25 = qJD(2) * t27;
t8 = -t25 + (qJD(4) + t51) * qJD(2);
t69 = t8 * t30 + t52 * t7;
t68 = -t31 * t34 + t32 * t36;
t35 = cos(qJ(5));
t67 = t35 ^ 2;
t33 = sin(qJ(5));
t66 = t33 / 0.2e1;
t65 = -t35 / 0.2e1;
t63 = Ifges(6,4) * t33;
t62 = Ifges(6,4) * t35;
t21 = t31 * t36 + t32 * t34;
t14 = t21 * qJD(2);
t11 = qJD(1) * t14;
t61 = t11 * t68;
t58 = -mrSges(5,2) + mrSges(4,1);
t57 = Ifges(6,5) * qJD(5);
t56 = Ifges(6,6) * qJD(5);
t54 = qJD(2) * mrSges(6,3);
t50 = -t32 * pkin(2) - pkin(3);
t9 = t32 * t26 - t27;
t45 = qJD(4) - t9;
t5 = (-pkin(3) - pkin(6)) * qJD(2) + t45;
t3 = -t33 * qJD(3) + t35 * t5;
t1 = t3 * qJD(5) + t33 * t11;
t4 = t35 * qJD(3) + t33 * t5;
t2 = -t4 * qJD(5) + t35 * t11;
t49 = -t1 * t33 - t2 * t35;
t16 = t68 * qJD(2);
t48 = t7 * t16 + t8 * t21;
t47 = t3 * t35 + t33 * t4;
t46 = -t3 * t33 + t35 * t4;
t44 = mrSges(6,1) * t35 - mrSges(6,2) * t33;
t43 = t33 * mrSges(6,1) + t35 * mrSges(6,2);
t23 = -qJD(5) * mrSges(6,2) - t33 * t54;
t24 = qJD(5) * mrSges(6,1) - t35 * t54;
t41 = t33 * t23 + t35 * t24;
t40 = (t35 * t23 - t33 * t24) * qJD(5);
t39 = t46 * qJD(5) - t49;
t38 = m(6) * t39;
t37 = qJD(2) ^ 2;
t29 = -pkin(6) + t50;
t22 = t43 * qJD(2);
t19 = t44 * qJD(5) * qJD(2);
t18 = t57 + (Ifges(6,1) * t35 - t63) * qJD(2);
t17 = t56 + (-Ifges(6,2) * t33 + t62) * qJD(2);
t13 = t21 * qJD(1);
t12 = qJD(2) * t51 - t25;
t6 = -qJD(2) * pkin(3) + t45;
t20 = [t16 * t22 + t21 * t19 + (-t34 * mrSges(3,1) - t36 * mrSges(3,2)) * t37 + t41 * t14 - t68 * t40 + m(4) * (t10 * t16 + t12 * t21 - t9 * t14 - t61) + m(5) * (t6 * t14 + t48 - t61) + m(6) * (t47 * t14 - t39 * t68 + t48) + ((-mrSges(4,2) + mrSges(5,3)) * t16 - t58 * t14) * qJD(2); -t12 * mrSges(4,2) + t30 * t19 + (mrSges(5,3) + t43) * t8 + t52 * t22 - t41 * t13 - t58 * t11 + t49 * mrSges(6,3) + ((t7 * mrSges(6,1) - t17 / 0.2e1 + t29 * t23 - t4 * mrSges(6,3) - t56 / 0.2e1) * t35 + (-t7 * mrSges(6,2) - t18 / 0.2e1 - t29 * t24 + t3 * mrSges(6,3) - t57 / 0.2e1) * t33) * qJD(5) + t29 * t38 + (t15 * mrSges(4,2) + t52 * mrSges(5,3) + t58 * t13 + (-0.3e1 / 0.2e1 * t67 * Ifges(6,4) + (0.3e1 / 0.2e1 * t63 + (-0.3e1 / 0.2e1 * Ifges(6,1) + 0.3e1 / 0.2e1 * Ifges(6,2)) * t35) * t33) * qJD(5)) * qJD(2) + (-t47 * t13 + t69) * m(6) + (t11 * t50 - t6 * t13 + t69) * m(5) + ((-t11 * t32 + t12 * t31) * pkin(2) - t10 * t15 + t9 * t13) * m(4); m(6) * (t1 * t35 - t2 * t33) + (-m(6) * t47 + (-t33 ^ 2 - t67) * t54 - t41) * qJD(5); -t37 * mrSges(5,3) + t40 + m(5) * t11 + t38 + (-t22 + (-m(5) - m(6)) * t7) * qJD(2); t2 * mrSges(6,1) - t1 * mrSges(6,2) - t3 * t23 + t4 * t24 + (-t7 * t44 + t18 * t66 + t35 * t17 / 0.2e1 + ((-Ifges(6,1) * t33 - t62) * t65 + (-Ifges(6,2) * t35 - t63) * t66) * qJD(2) + t46 * mrSges(6,3) + (-Ifges(6,5) * t33 / 0.2e1 + Ifges(6,6) * t65) * qJD(5)) * qJD(2);];
tauc = t20(:);
