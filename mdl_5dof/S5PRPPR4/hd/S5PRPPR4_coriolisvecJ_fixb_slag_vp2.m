% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPPR4
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:46
% EndTime: 2019-12-31 17:36:48
% DurationCPUTime: 0.75s
% Computational Cost: add. (478->122), mult. (1248->182), div. (0->0), fcn. (777->4), ass. (0->60)
t46 = sin(pkin(8));
t47 = cos(pkin(8));
t80 = t46 ^ 2 + t47 ^ 2;
t82 = t80 * (mrSges(5,2) + mrSges(4,3));
t48 = sin(qJ(5));
t49 = cos(qJ(5));
t65 = t46 * qJD(2);
t68 = qJD(2) * t47;
t25 = -t48 * t65 - t49 * t68;
t79 = -t25 / 0.2e1;
t54 = t46 * t48 + t47 * t49;
t23 = t54 * qJD(5);
t16 = qJD(2) * t23;
t78 = t54 * t16;
t30 = t46 * t49 - t47 * t48;
t24 = t30 * qJD(5);
t17 = qJD(2) * t24;
t77 = t30 * t17;
t59 = t46 * qJ(4) + pkin(2);
t75 = (t47 * pkin(3) + t59) * qJD(2);
t26 = -t48 * t68 + t49 * t65;
t73 = t26 / 0.2e1;
t72 = Ifges(6,4) * t26;
t71 = -pkin(6) + qJ(3);
t64 = qJD(2) * qJ(3);
t33 = t46 * qJD(1) + t47 * t64;
t27 = (pkin(3) + pkin(4)) * t47 + t59;
t69 = qJD(2) * t27;
t67 = qJD(3) * t46;
t63 = qJD(2) * qJD(3);
t62 = t80 * qJ(3) * t63 + t33 * t47 * qJD(3);
t60 = qJD(4) * t65;
t32 = t47 * qJD(1) - t46 * t64;
t55 = t17 * mrSges(6,1) - t16 * mrSges(6,2);
t28 = qJD(4) - t32;
t15 = -pkin(6) * t65 + t28;
t18 = -pkin(6) * t68 + t33;
t5 = t49 * t15 - t48 * t18;
t6 = t48 * t15 + t49 * t18;
t34 = t71 * t46;
t35 = t71 * t47;
t10 = t49 * t34 - t48 * t35;
t11 = t48 * t34 + t49 * t35;
t31 = (-t47 * mrSges(5,1) - t46 * mrSges(5,3)) * qJD(2);
t52 = t30 * qJD(3);
t51 = t54 * qJD(3);
t22 = t33 * t68;
t20 = qJD(3) - t75;
t19 = Ifges(6,4) * t25;
t14 = qJD(5) * mrSges(6,1) - t26 * mrSges(6,3);
t13 = -qJD(5) * mrSges(6,2) + t25 * mrSges(6,3);
t12 = -qJD(3) + t69;
t9 = -t25 * mrSges(6,1) + t26 * mrSges(6,2);
t8 = Ifges(6,1) * t26 + Ifges(6,5) * qJD(5) + t19;
t7 = Ifges(6,2) * t25 + Ifges(6,6) * qJD(5) + t72;
t4 = -t11 * qJD(5) + t52;
t3 = t10 * qJD(5) + t51;
t2 = qJD(2) * t52 - t6 * qJD(5);
t1 = qJD(2) * t51 + t5 * qJD(5);
t21 = [m(6) * (t1 * t30 - t2 * t54 - t6 * t23 - t5 * t24) - t23 * t13 - t24 * t14 + (-t77 - t78) * mrSges(6,3); -t23 * t8 / 0.2e1 + t12 * (t24 * mrSges(6,1) - t23 * mrSges(6,2)) + qJD(5) * (-Ifges(6,5) * t23 - Ifges(6,6) * t24) / 0.2e1 - t24 * t7 / 0.2e1 + t3 * t13 + t4 * t14 + t27 * t55 + m(4) * (-t32 * t67 + t62) + m(5) * (t28 * t67 + t62) + m(6) * (t1 * t11 + t2 * t10 + t6 * t3 + t5 * t4) + (mrSges(6,1) * t54 + t30 * mrSges(6,2)) * t60 + (t9 - 0.2e1 * t31 + m(5) * (-t20 + t75) + m(6) * (t12 + t69)) * qJD(4) * t46 + (t17 * t54 + t24 * t79) * Ifges(6,2) + (-t16 * t30 - t23 * t73) * Ifges(6,1) + (-t1 * t54 + t10 * t16 - t11 * t17 - t2 * t30 + t5 * t23 - t6 * t24) * mrSges(6,3) + (t23 * t79 - t24 * t73 - t77 + t78) * Ifges(6,4) + 0.2e1 * t63 * t82; t25 * t13 - t26 * t14 + (-m(5) - m(6)) * t60 - m(4) * (-t32 * t65 + t22) - m(5) * (t28 * t65 + t22) - m(6) * (-t6 * t25 + t5 * t26) - t55 - qJD(2) ^ 2 * t82; m(6) * (t1 * t48 + t2 * t49) + (t49 * t16 - t48 * t17) * mrSges(6,3) + (-m(6) * t12 + t31 - t9 + (qJD(3) + t20) * m(5)) * t65 + (m(6) * (-t48 * t5 + t49 * t6) + t49 * t13 - t48 * t14) * qJD(5); -Ifges(6,5) * t16 - Ifges(6,6) * t17 - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t12 * (t26 * mrSges(6,1) + t25 * mrSges(6,2)) - t26 * (Ifges(6,1) * t25 - t72) / 0.2e1 + t7 * t73 - qJD(5) * (Ifges(6,5) * t25 - Ifges(6,6) * t26) / 0.2e1 - t5 * t13 + t6 * t14 + (t5 * t25 + t6 * t26) * mrSges(6,3) + (-Ifges(6,2) * t26 + t19 + t8) * t79;];
tauc = t21(:);
