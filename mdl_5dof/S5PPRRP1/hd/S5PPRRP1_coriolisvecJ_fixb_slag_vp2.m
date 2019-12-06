% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:39
% EndTime: 2019-12-05 15:06:41
% DurationCPUTime: 0.96s
% Computational Cost: add. (633->168), mult. (1775->222), div. (0->0), fcn. (1076->6), ass. (0->71)
t94 = -qJD(4) / 0.2e1;
t93 = Ifges(5,4) + Ifges(6,4);
t49 = sin(pkin(8));
t50 = cos(pkin(8));
t52 = sin(qJ(3));
t54 = cos(qJ(3));
t29 = t52 * t49 - t50 * t54;
t19 = t29 * qJD(1);
t53 = cos(qJ(4));
t43 = -pkin(4) * t53 - pkin(3);
t10 = qJD(3) * t43 + qJD(5) + t19;
t13 = -qJD(3) * pkin(3) + t19;
t51 = sin(qJ(4));
t31 = (-t53 * mrSges(6,1) + t51 * mrSges(6,2)) * qJD(3);
t30 = t49 * t54 + t52 * t50;
t20 = t30 * qJD(1);
t14 = qJD(3) * pkin(6) + t20;
t61 = qJ(5) * qJD(3) + t14;
t73 = qJD(2) * t51;
t7 = t53 * t61 + t73;
t83 = m(6) * t10;
t9 = t14 * t53 + t73;
t91 = (t31 + t83) * pkin(4) + t13 * mrSges(5,1) + t10 * mrSges(6,1) - t9 * mrSges(5,3) - t7 * mrSges(6,3) - ((Ifges(5,2) + Ifges(6,2)) * t53 + t93 * t51) * qJD(3) / 0.2e1 + (Ifges(6,6) + Ifges(5,6)) * t94;
t71 = qJD(3) * t53;
t85 = -t93 * t71 / 0.2e1;
t84 = m(5) * t9;
t22 = t30 * qJD(3);
t16 = qJD(1) * t22;
t82 = t16 * t29;
t21 = t29 * qJD(3);
t81 = t21 * t51;
t80 = t21 * t53;
t77 = -qJ(5) - pkin(6);
t72 = qJD(3) * t51;
t33 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t72;
t34 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t72;
t76 = -t33 - t34;
t35 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t71;
t36 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t71;
t75 = t35 + t36;
t67 = qJD(3) * qJD(4);
t27 = (mrSges(6,1) * t51 + mrSges(6,2) * t53) * t67;
t70 = qJD(4) * t51;
t69 = qJD(4) * t53;
t68 = qJD(5) * t53;
t46 = t53 * qJD(2);
t66 = 0.3e1 / 0.2e1 * Ifges(6,4) + 0.3e1 / 0.2e1 * Ifges(5,4);
t65 = Ifges(5,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t64 = m(5) * pkin(6) + mrSges(5,3);
t8 = -t14 * t51 + t46;
t63 = -m(5) * t8 - t34;
t62 = qJD(4) * t77;
t60 = -t31 + (t53 * mrSges(5,1) - t51 * mrSges(5,2) + mrSges(4,1)) * qJD(3);
t59 = (-Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * qJD(4);
t15 = qJD(1) * t21;
t3 = qJD(4) * t46 - t14 * t70 - t53 * t15;
t6 = -t51 * t61 + t46;
t5 = qJD(4) * pkin(4) + t6;
t58 = m(6) * (-t5 * t51 + t53 * t7);
t57 = t51 * t76 + t53 * t75;
t56 = -t13 * mrSges(5,2) - t10 * mrSges(6,2) + t8 * mrSges(5,3) + t5 * mrSges(6,3) + t85 - (Ifges(5,1) + Ifges(6,1)) * t72 / 0.2e1 + (Ifges(5,5) + Ifges(6,5)) * t94;
t38 = t77 * t53;
t37 = t77 * t51;
t28 = (mrSges(5,1) * t51 + mrSges(5,2) * t53) * t67;
t18 = -qJD(5) * t51 + t53 * t62;
t17 = t51 * t62 + t68;
t11 = (pkin(4) * t70 + t20) * qJD(3);
t4 = -qJD(4) * t9 + t15 * t51;
t2 = (-qJD(3) * qJD(5) + t15) * t51 - t7 * qJD(4);
t1 = (-qJ(5) * t70 + t68) * qJD(3) + t3;
t12 = [(t27 + t28) * t29 - t60 * t22 - (-qJD(3) * mrSges(4,2) + t57) * t21 + m(4) * (t19 * t22 - t20 * t21 + t82) + m(5) * (t13 * t22 + t8 * t81 - t9 * t80 + t82) + m(6) * (t10 * t22 + t11 * t29 + t5 * t81 - t7 * t80) + ((-t51 * t75 + t53 * t76) * qJD(4) - m(4) * t15 + m(5) * (t3 * t53 - t4 * t51 - t8 * t69 - t9 * t70) + m(6) * (t1 * t53 - t2 * t51 - t5 * t69 - t7 * t70)) * t30; m(5) * (t3 * t51 + t4 * t53) + m(6) * (t1 * t51 + t2 * t53) + (m(5) * (-t51 * t8 + t53 * t9) + t58 + t57 + (mrSges(5,3) + mrSges(6,3)) * qJD(3) * (-t51 ^ 2 - t53 ^ 2)) * qJD(4); -t16 * mrSges(4,1) + t17 * t35 + t18 * t33 + t43 * t27 + (-t19 * qJD(3) + t15) * mrSges(4,2) + m(6) * (-t1 * t38 + t11 * t43 + t17 * t7 + t18 * t5 + t2 * t37) + (t16 * mrSges(5,2) + t11 * mrSges(6,2) - t2 * mrSges(6,3) - t64 * t4 - (m(6) * t5 + t33 - t63) * t19 + (t59 + (t38 * mrSges(6,3) - t51 * t66) * qJD(3) + (-t36 - t84) * pkin(6) + t91) * qJD(4)) * t51 + (-t16 * mrSges(5,1) - t11 * mrSges(6,1) + t1 * mrSges(6,3) + t64 * t3 - (-m(6) * t7 - t75 - t84) * t19 + (t65 * qJD(4) + t63 * pkin(6) + (-t37 * mrSges(6,3) + t66 * t53 + (0.3e1 / 0.2e1 * Ifges(6,1) + 0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(6,2) - 0.3e1 / 0.2e1 * Ifges(5,2)) * t51) * qJD(3) - t56) * qJD(4)) * t53 + (-m(5) * t16 - t28) * pkin(3) + (-m(5) * t13 + t60 - t83) * t20; t4 * mrSges(5,1) - t3 * mrSges(5,2) - t1 * mrSges(6,2) + t9 * t34 - t6 * t35 - t8 * t36 + (m(6) * pkin(4) + mrSges(6,1)) * t2 + (t33 - m(6) * (-t5 + t6)) * t7 + (((-pkin(4) * mrSges(6,3) + t65) * qJD(4) + t56 + t85) * t53 + ((Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1) * t72 + t59 + (-Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t71 - t91) * t51) * qJD(3); m(6) * t11 + (t51 * t33 - t53 * t35 - t58) * qJD(3) + t27;];
tauc = t12(:);
