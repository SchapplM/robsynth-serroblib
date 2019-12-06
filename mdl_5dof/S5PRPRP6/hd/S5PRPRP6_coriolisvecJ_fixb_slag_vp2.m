% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:21
% EndTime: 2019-12-05 15:40:25
% DurationCPUTime: 0.94s
% Computational Cost: add. (522->165), mult. (1175->226), div. (0->0), fcn. (473->4), ass. (0->69)
t85 = mrSges(5,1) + mrSges(6,1);
t84 = -mrSges(6,2) - mrSges(5,3);
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t58 = qJD(2) * t35;
t63 = t85 * qJD(4) + t84 * t58;
t60 = qJD(2) * t33;
t26 = -mrSges(6,2) * t60 + qJD(4) * mrSges(6,3);
t64 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t60 + t26;
t83 = t33 * t63 - t35 * t64;
t37 = -pkin(2) - pkin(6);
t36 = cos(qJ(2));
t54 = t36 * qJD(1);
t48 = qJD(3) - t54;
t16 = qJD(2) * t37 + t48;
t80 = t16 * (t33 ^ 2 + t35 ^ 2);
t79 = 2 * m(6);
t77 = -m(5) / 0.2e1;
t7 = qJD(4) * qJ(5) + t16 * t33;
t76 = -t7 / 0.2e1;
t20 = pkin(4) * t33 - qJ(5) * t35 + qJ(3);
t34 = sin(qJ(2));
t62 = qJD(1) * t34;
t6 = qJD(2) * t20 + t62;
t75 = m(6) * t6;
t59 = qJD(2) * t34;
t49 = qJD(1) * t59;
t57 = qJD(4) * t33;
t3 = t16 * t57 - t35 * t49;
t72 = t3 * t35;
t56 = qJD(4) * t35;
t4 = t16 * t56 + t33 * t49;
t71 = Ifges(6,1) * t35;
t70 = Ifges(5,4) * t33;
t69 = Ifges(5,4) * t35;
t68 = Ifges(6,5) * t33;
t27 = qJD(2) * qJ(3) + t62;
t66 = t27 * t36;
t18 = (t33 * mrSges(6,1) - t35 * mrSges(6,3)) * qJD(2);
t19 = (t33 * mrSges(5,1) + t35 * mrSges(5,2)) * qJD(2);
t65 = t18 + t19;
t55 = t27 * qJD(2);
t53 = qJD(2) * qJD(4);
t52 = -Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1;
t51 = 0.3e1 / 0.2e1 * Ifges(6,5) - 0.3e1 / 0.2e1 * Ifges(5,4);
t50 = -Ifges(5,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t5 = -qJD(4) * pkin(4) - t16 * t35 + qJD(5);
t47 = m(6) * t5 - t63;
t2 = qJD(4) * qJD(5) + t4;
t46 = t2 * t33 - t72;
t45 = t33 * t4 - t72;
t44 = t33 * t5 + t35 * t7;
t43 = mrSges(5,1) * t35 - mrSges(5,2) * t33;
t42 = mrSges(6,1) * t35 + mrSges(6,3) * t33;
t41 = pkin(4) * t35 + qJ(5) * t33;
t8 = qJD(4) * t41 - qJD(5) * t35 + qJD(3);
t39 = m(6) * t44 - t83;
t38 = qJD(2) ^ 2;
t30 = Ifges(6,5) * t58;
t22 = -qJD(2) * pkin(2) + t48;
t21 = (qJD(3) + t54) * qJD(2);
t15 = t43 * t53;
t14 = t42 * t53;
t13 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t35 - t70) * qJD(2);
t12 = Ifges(6,4) * qJD(4) + (t68 + t71) * qJD(2);
t11 = Ifges(5,6) * qJD(4) + (-Ifges(5,2) * t33 + t69) * qJD(2);
t10 = Ifges(6,6) * qJD(4) + Ifges(6,3) * t60 + t30;
t1 = (t8 + t54) * qJD(2);
t9 = [(t14 + t15 + (-mrSges(3,1) + mrSges(4,2)) * t38 + (t33 * t64 + t35 * t63) * qJD(2) + m(4) * (qJD(2) * t22 + t21) + m(5) * (qJD(2) * t80 + t21) + m(6) * (-t5 * t58 + t60 * t7 + t1)) * t34 + ((-mrSges(3,2) + mrSges(4,3)) * t38 + t65 * qJD(2) + t83 * qJD(4) + m(4) * (-t49 + t55) + m(5) * (-t45 + t55) + m(6) * (qJD(2) * t6 - t5 * t57 - t56 * t7 - t46)) * t36; qJ(3) * t15 + qJD(3) * t19 + t20 * t14 + t8 * t18 - t65 * t54 + (qJD(2) * t48 + t21) * mrSges(4,3) + m(6) * (t1 * t20 + t6 * t8) + 0.2e1 * (-t36 * t75 / 0.2e1 + (t34 * t80 + t66) * t77 - (pkin(2) * t59 + t22 * t34 + t66) * m(4) / 0.2e1) * qJD(1) + (t21 * mrSges(5,2) - t1 * mrSges(6,3) + t47 * t62 + (0.2e1 * (-m(6) / 0.2e1 + t77) * t37 - t84) * t3 + (t6 * mrSges(6,1) + t27 * mrSges(5,1) + t10 / 0.2e1 - t11 / 0.2e1 - t7 * mrSges(6,2) + t50 * qJD(4) + t51 * t58 + (m(6) * t7 + t64) * t37) * qJD(4)) * t35 + (t21 * mrSges(5,1) + t1 * mrSges(6,1) - t2 * mrSges(6,2) - t64 * t62 + (t2 * t37 / 0.2e1 + t62 * t76) * t79 + (-t12 / 0.2e1 - t13 / 0.2e1 - t27 * mrSges(5,2) + t6 * mrSges(6,3) - t5 * mrSges(6,2) + t52 * qJD(4) + t47 * t37 + (-t51 * t33 + (-0.3e1 / 0.2e1 * Ifges(6,1) - 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(6,3) + 0.3e1 / 0.2e1 * Ifges(5,2)) * t35) * qJD(2)) * qJD(4) + (m(5) * t37 - mrSges(5,3)) * t4) * t33 + (m(4) + m(5)) * (qJ(3) * t21 + qJD(3) * t27); -t38 * mrSges(4,3) + m(5) * t45 + m(6) * t46 + t39 * qJD(4) + (-m(5) * t27 - t75 + (-t27 + t62) * m(4) - t65) * qJD(2); -t4 * mrSges(5,2) + t2 * mrSges(6,3) + qJD(5) * t26 - t85 * t3 + m(6) * (-t3 * pkin(4) + t2 * qJ(5) + t7 * qJD(5)) - t39 * t16 + (t35 * t11 / 0.2e1 - t27 * t43 - t6 * t42 - (Ifges(6,3) * t35 - t68) * t60 / 0.2e1 + t44 * mrSges(6,2) + ((-qJ(5) * mrSges(6,2) + t50) * t35 + (pkin(4) * mrSges(6,2) + t52) * t33) * qJD(4) - (t30 + t10 + (-Ifges(5,1) * t33 - t69) * qJD(2)) * t35 / 0.2e1 + (t12 + t13 + (-Ifges(5,2) * t35 - t70 + t71) * qJD(2)) * t33 / 0.2e1) * qJD(2) + (-t18 - t75) * t41 * qJD(2); -qJD(4) * t26 + (-mrSges(6,2) * t57 + t35 * t18) * qJD(2) + (t3 / 0.2e1 + t6 * t58 / 0.2e1 + qJD(4) * t76) * t79;];
tauc = t9(:);
