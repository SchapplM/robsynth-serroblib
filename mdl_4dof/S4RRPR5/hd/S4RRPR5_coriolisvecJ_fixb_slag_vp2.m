% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:28
% EndTime: 2019-12-31 17:03:29
% DurationCPUTime: 0.38s
% Computational Cost: add. (393->94), mult. (704->135), div. (0->0), fcn. (244->4), ass. (0->59)
t83 = m(4) + m(5);
t28 = qJD(1) + qJD(2);
t35 = -pkin(2) - pkin(6);
t34 = cos(qJ(2));
t67 = pkin(1) * qJD(1);
t50 = -t34 * t67 + qJD(3);
t16 = t35 * t28 + t50;
t33 = cos(qJ(4));
t32 = sin(qJ(2));
t66 = pkin(1) * qJD(2);
t57 = qJD(1) * t66;
t53 = t32 * t57;
t31 = sin(qJ(4));
t62 = qJD(4) * t31;
t1 = -t16 * t62 + t33 * t53;
t61 = qJD(4) * t33;
t2 = t16 * t61 + t31 * t53;
t51 = t1 * t33 + t2 * t31;
t82 = m(5) * t51;
t71 = t31 * mrSges(5,3);
t21 = -qJD(4) * mrSges(5,2) - t28 * t71;
t69 = t33 * mrSges(5,3);
t23 = qJD(4) * mrSges(5,1) - t28 * t69;
t81 = t31 * t21 + t33 * t23;
t48 = mrSges(5,1) * t31 + mrSges(5,2) * t33;
t17 = t48 * t28;
t80 = mrSges(4,3) * t28 + t17;
t79 = (t31 ^ 2 + t33 ^ 2) * t16;
t77 = Ifges(5,4) * t31;
t76 = Ifges(5,4) * t33;
t75 = Ifges(5,5) * t31;
t74 = Ifges(5,6) * t33;
t22 = qJ(3) * t28 + t32 * t67;
t72 = t22 * t34;
t65 = qJD(2) * t32;
t64 = qJD(3) * t28;
t60 = pkin(1) * t65;
t59 = t34 * t66;
t58 = -pkin(1) * t34 - pkin(2);
t54 = t28 * t60;
t52 = t34 * t57;
t49 = mrSges(5,1) * t33 - mrSges(5,2) * t31;
t46 = t21 * t33 - t23 * t31;
t19 = t52 + t64;
t45 = t19 * qJ(3) + t22 * qJD(3);
t44 = t31 * (-Ifges(5,2) * t33 - t77);
t43 = t33 * (-Ifges(5,1) * t31 - t76);
t42 = (Ifges(5,1) * t33 - t77) * t28;
t41 = (-Ifges(5,2) * t31 + t76) * t28;
t40 = t49 * qJD(4);
t39 = t46 * qJD(4);
t11 = Ifges(5,6) * qJD(4) + t41;
t12 = Ifges(5,5) * qJD(4) + t42;
t36 = t19 * t48 + t22 * t40 + qJD(4) ^ 2 * (-t74 - t75) / 0.2e1 + mrSges(4,2) * t53 + (t43 - t44) * qJD(4) * t28 - (t12 + t42) * t62 / 0.2e1 - (t11 + t41) * t61 / 0.2e1;
t27 = pkin(1) * t32 + qJ(3);
t25 = qJD(3) + t59;
t20 = -pkin(2) * t28 + t50;
t13 = t28 * t40;
t3 = [mrSges(4,2) * t54 + t19 * mrSges(4,3) - t1 * t69 + t27 * t13 - t2 * t71 + t36 + t83 * (t19 * t27 + t22 * t25) + (m(4) * (t58 * qJD(1) + t20) + m(5) * t79 + t81) * t60 + (t21 * t61 - t23 * t62 + t82) * (-pkin(6) + t58) + t80 * t25 + (-t28 * t59 - t52) * mrSges(3,2) + (-t53 - t54) * mrSges(3,1); t36 + ((-qJD(2) * mrSges(3,2) - t17) * t34 + (-qJD(2) * mrSges(3,1) - t81) * t32 - m(5) * (t32 * t79 + t72) + ((mrSges(3,2) - mrSges(4,3)) * t34 + (mrSges(3,1) - mrSges(4,2)) * t32) * t28 + (-pkin(2) * t65 - t20 * t32 - t72) * m(4)) * t67 + m(5) * (t51 * t35 + t45) + m(4) * t45 - t51 * mrSges(5,3) + (t19 + t64) * mrSges(4,3) + t35 * t39 + qJ(3) * t13 + qJD(3) * t17; m(4) * t53 + t82 + t39 + (-t83 * t22 - t80) * t28; t1 * mrSges(5,1) - t2 * mrSges(5,2) - t46 * t16 + (-t22 * t49 + t31 * t12 / 0.2e1 + t33 * t11 / 0.2e1 + (-t43 / 0.2e1 + t44 / 0.2e1) * t28 + (-t75 / 0.2e1 - t74 / 0.2e1) * qJD(4)) * t28;];
tauc = t3(:);
