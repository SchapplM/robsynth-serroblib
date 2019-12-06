% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRRP2
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:32
% EndTime: 2019-12-05 15:08:36
% DurationCPUTime: 0.87s
% Computational Cost: add. (628->164), mult. (1746->242), div. (0->0), fcn. (1059->6), ass. (0->71)
t85 = mrSges(5,1) + mrSges(6,1);
t82 = mrSges(6,2) + mrSges(5,3);
t42 = sin(pkin(8));
t43 = cos(pkin(8));
t45 = sin(qJ(3));
t47 = cos(qJ(3));
t26 = t45 * t42 - t43 * t47;
t44 = sin(qJ(4));
t62 = qJD(3) * t44;
t53 = t62 / 0.2e1;
t84 = -qJD(4) / 0.2e1;
t83 = qJD(4) / 0.2e1;
t46 = cos(qJ(4));
t27 = t42 * t47 + t45 * t43;
t17 = t27 * qJD(1);
t12 = qJD(3) * pkin(6) + t17;
t73 = t12 * t44;
t8 = qJD(2) * t46 - t73;
t63 = -t8 + qJD(5);
t81 = 2 * m(5);
t80 = 2 * m(6);
t79 = pkin(6) / 0.2e1;
t18 = t26 * qJD(3);
t13 = qJD(1) * t18;
t9 = qJD(2) * t44 + t12 * t46;
t3 = t9 * qJD(4) - t13 * t44;
t78 = t3 / 0.2e1;
t16 = t26 * qJD(1);
t77 = t16 / 0.2e1;
t76 = -t17 / 0.2e1;
t75 = t3 * t44;
t74 = t3 * t46;
t19 = t27 * qJD(3);
t14 = qJD(1) * t19;
t72 = t14 * t26;
t71 = t18 * t44;
t70 = t18 * t46;
t59 = qJD(4) * t46;
t67 = qJD(2) * t59 - t46 * t13;
t66 = t85 * qJD(4) - t82 * t62;
t61 = qJD(3) * t46;
t34 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t61;
t35 = mrSges(6,2) * t61 + qJD(4) * mrSges(6,3);
t65 = t34 + t35;
t60 = qJD(4) * t44;
t58 = qJD(3) * qJD(4);
t57 = pkin(6) * t78;
t56 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t55 = 0.3e1 / 0.2e1 * Ifges(6,5) - 0.3e1 / 0.2e1 * Ifges(5,4);
t54 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t28 = (-t46 * mrSges(6,1) - t44 * mrSges(6,3)) * qJD(3);
t52 = -t28 + (t46 * mrSges(5,1) - t44 * mrSges(5,2) + mrSges(4,1)) * qJD(3);
t51 = pkin(4) * t44 - qJ(5) * t46;
t31 = -pkin(4) * t46 - qJ(5) * t44 - pkin(3);
t50 = -t66 * t44 + t65 * t46;
t15 = t51 * qJD(4) - qJD(5) * t44;
t11 = -qJD(3) * pkin(3) + t16;
t38 = Ifges(6,5) * t62;
t6 = qJD(4) * qJ(5) + t9;
t7 = t31 * qJD(3) + t16;
t49 = t11 * mrSges(5,1) + t7 * mrSges(6,1) + Ifges(6,6) * t83 - Ifges(6,3) * t61 / 0.2e1 + t38 / 0.2e1 + Ifges(5,6) * t84 - (Ifges(5,4) * t44 + Ifges(5,2) * t46) * qJD(3) / 0.2e1 - t6 * mrSges(6,2) - t9 * mrSges(5,3);
t39 = Ifges(5,4) * t61;
t5 = -qJD(4) * pkin(4) + t63;
t48 = t11 * mrSges(5,2) + t5 * mrSges(6,2) + (Ifges(6,1) * t44 - Ifges(6,5) * t46) * qJD(3) / 0.2e1 + Ifges(5,1) * t53 + t39 / 0.2e1 - t7 * mrSges(6,3) - t8 * mrSges(5,3) + (Ifges(6,4) + Ifges(5,5)) * t83;
t30 = t51 * qJD(3);
t25 = (mrSges(5,1) * t44 + mrSges(5,2) * t46) * t58;
t24 = (mrSges(6,1) * t44 - mrSges(6,3) * t46) * t58;
t4 = (t17 + t15) * qJD(3);
t2 = -t12 * t60 + t67;
t1 = (qJD(5) - t73) * qJD(4) + t67;
t10 = [(t24 + t25) * t26 - t52 * t19 - (-qJD(3) * mrSges(4,2) + t50) * t18 + m(4) * (t16 * t19 - t17 * t18 + t72) + m(5) * (t11 * t19 - t9 * t70 + t8 * t71 + t72) + m(6) * (t19 * t7 + t26 * t4 - t5 * t71 - t6 * t70) + ((-t65 * t44 - t66 * t46) * qJD(4) - m(4) * t13 + m(5) * (t2 * t46 - t8 * t59 - t9 * t60 + t75) + m(6) * (t1 * t46 + t5 * t59 - t6 * t60 + t75)) * t27; m(5) * (t2 * t44 - t74) + m(6) * (t1 * t44 - t74) + (m(5) * (-t44 * t8 + t46 * t9) + m(6) * (t44 * t5 + t46 * t6) + t50 + t82 * qJD(3) * (-t44 ^ 2 - t46 ^ 2)) * qJD(4); -t14 * mrSges(4,1) - pkin(3) * t25 + t15 * t28 + t31 * t24 + (-t16 * qJD(3) + t13) * mrSges(4,2) + t52 * t17 + (t31 * t4 / 0.2e1 + (t15 / 0.2e1 + t76) * t7) * t80 + (-t14 * pkin(3) / 0.2e1 + t11 * t76) * t81 + (-t14 * mrSges(5,1) - t4 * mrSges(6,1) + t1 * mrSges(6,2) + t2 * mrSges(5,3) + t65 * t16 + (t1 * t79 + t6 * t77) * t80 + (t2 * t79 + t77 * t9) * t81 + (t56 * qJD(4) - t55 * t61 + (-m(5) * t8 + m(6) * t5 - t66) * pkin(6) + t48) * qJD(4)) * t46 + (t14 * mrSges(5,2) - t4 * mrSges(6,3) + t82 * t3 - t66 * t16 + (t5 * t77 + t57) * t80 + (t57 - t16 * t8 / 0.2e1) * t81 + (t54 * qJD(4) + (-m(5) * t9 - m(6) * t6 - t65) * pkin(6) + (t55 * t44 + (-0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(6,3) + 0.3e1 / 0.2e1 * Ifges(6,1) + 0.3e1 / 0.2e1 * Ifges(5,1)) * t46) * qJD(3) + t49) * qJD(4)) * t44; -t2 * mrSges(5,2) + t1 * mrSges(6,3) - t30 * t28 - t8 * t34 + t66 * t9 + t63 * t35 - t85 * t3 + ((-t39 / 0.2e1 + Ifges(6,5) * t61 / 0.2e1 + (-pkin(4) * mrSges(6,2) + t56) * qJD(4) - t48) * t46 + (-t38 / 0.2e1 + Ifges(5,4) * t53 + (-qJ(5) * mrSges(6,2) + t54) * qJD(4) + (-Ifges(6,1) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 - Ifges(5,1) / 0.2e1) * t61 - t49) * t44) * qJD(3) + (-t3 * pkin(4) + t1 * qJ(5) - t7 * t30 - t5 * t9 + t63 * t6) * m(6); -qJD(4) * t35 + (mrSges(6,2) * t59 + t28 * t44) * qJD(3) + (t7 * t53 + t6 * t84 + t78) * t80;];
tauc = t10(:);
