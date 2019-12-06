% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:00
% EndTime: 2019-12-05 15:03:02
% DurationCPUTime: 0.55s
% Computational Cost: add. (434->95), mult. (1091->137), div. (0->0), fcn. (658->6), ass. (0->51)
t27 = sin(pkin(8));
t29 = sin(qJ(3));
t49 = cos(pkin(8));
t55 = cos(qJ(3));
t65 = -t29 * t27 + t55 * t49;
t64 = t65 * qJD(1);
t19 = t27 * t55 + t29 * t49;
t14 = t19 * qJD(3);
t10 = qJD(1) * t14;
t28 = sin(qJ(5));
t30 = cos(qJ(5));
t31 = -pkin(3) - pkin(6);
t36 = -t64 + qJD(4);
t5 = qJD(3) * t31 + t36;
t3 = -t28 * qJD(2) + t30 * t5;
t1 = qJD(5) * t3 + t28 * t10;
t4 = t30 * qJD(2) + t28 * t5;
t2 = -qJD(5) * t4 + t30 * t10;
t63 = t1 * t28 + t2 * t30;
t41 = mrSges(6,1) * t30 - mrSges(6,2) * t28;
t43 = t28 * t3 - t30 * t4;
t53 = Ifges(6,4) * t30;
t54 = Ifges(6,4) * t28;
t59 = -t30 / 0.2e1;
t60 = -t28 / 0.2e1;
t12 = t19 * qJD(1);
t8 = qJD(3) * qJ(4) + t12;
t62 = t43 * mrSges(6,3) + (Ifges(6,6) * qJD(5) + (-Ifges(6,2) * t28 + t53) * qJD(3)) * t59 + (Ifges(6,5) * qJD(5) + (Ifges(6,1) * t30 - t54) * qJD(3)) * t60 + t41 * t8;
t61 = t30 ^ 2;
t56 = t8 * t64;
t52 = t10 * t65;
t50 = -mrSges(5,2) + mrSges(4,1);
t48 = qJD(3) * mrSges(6,3);
t13 = t65 * qJD(3);
t9 = t64 * qJD(3);
t6 = -qJD(3) * qJD(4) - t9;
t45 = t8 * t13 - t6 * t19;
t44 = t28 * t4 + t3 * t30;
t40 = t28 * mrSges(6,1) + t30 * mrSges(6,2);
t21 = -qJD(5) * mrSges(6,2) - t28 * t48;
t22 = qJD(5) * mrSges(6,1) - t30 * t48;
t39 = t30 * t21 - t28 * t22;
t38 = t28 * t21 + t30 * t22;
t37 = -t6 * qJ(4) + t8 * qJD(4);
t35 = t39 * qJD(5);
t34 = (Ifges(6,5) * t60 + Ifges(6,6) * t59) * qJD(5);
t33 = -t43 * qJD(5) + t63;
t20 = t40 * qJD(3);
t17 = t41 * qJD(5) * qJD(3);
t7 = -qJD(3) * pkin(3) + t36;
t11 = [t13 * t20 + t19 * t17 + t38 * t14 - t65 * t35 + m(4) * (t12 * t13 - t14 * t64 + t9 * t19 - t52) + m(5) * (t7 * t14 + t45 - t52) + m(6) * (t14 * t44 - t33 * t65 + t45) + (-t50 * t14 - (mrSges(4,2) - mrSges(5,3)) * t13) * qJD(3); m(6) * (t1 * t30 - t2 * t28) + (-m(6) * t44 + (-t28 ^ 2 - t61) * t48 - t38) * qJD(5); -t9 * mrSges(4,2) + qJ(4) * t17 + (-mrSges(5,3) - t40) * t6 + t36 * t20 - t38 * t12 - t50 * t10 - t63 * mrSges(6,3) - m(6) * (t12 * t44 + t56) + m(6) * (t63 * t31 + t37) + (t34 + (-m(6) * t43 + t39) * t31 + t62) * qJD(5) + (t64 * mrSges(4,2) + t36 * mrSges(5,3) + t50 * t12 + (-0.3e1 / 0.2e1 * t61 * Ifges(6,4) + (0.3e1 / 0.2e1 * t54 + (-0.3e1 / 0.2e1 * Ifges(6,1) + 0.3e1 / 0.2e1 * Ifges(6,2)) * t30) * t28) * qJD(5)) * qJD(3) + (-t10 * pkin(3) - t7 * t12 + t37 - t56) * m(5); t35 + m(5) * t10 + m(6) * t33 + (-qJD(3) * mrSges(5,3) - t20 + (-m(5) - m(6)) * t8) * qJD(3); t2 * mrSges(6,1) - t1 * mrSges(6,2) - t3 * t21 + t4 * t22 + (((-Ifges(6,1) * t28 - t53) * t59 + t28 * (-Ifges(6,2) * t30 - t54) / 0.2e1) * qJD(3) + t34 - t62) * qJD(3);];
tauc = t11(:);
