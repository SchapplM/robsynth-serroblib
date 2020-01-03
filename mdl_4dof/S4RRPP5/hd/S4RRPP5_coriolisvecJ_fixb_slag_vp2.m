% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:10
% EndTime: 2019-12-31 17:00:13
% DurationCPUTime: 0.96s
% Computational Cost: add. (393->168), mult. (1086->216), div. (0->0), fcn. (390->2), ass. (0->68)
t83 = -mrSges(3,1) + mrSges(4,2);
t48 = cos(qJ(2));
t63 = qJD(1) * t48;
t82 = pkin(3) * t63 + qJD(4);
t47 = sin(qJ(2));
t64 = qJD(1) * t47;
t41 = pkin(5) * t64;
t22 = -pkin(3) * t64 - t41;
t81 = qJD(3) - t22;
t43 = pkin(5) * t63;
t80 = -t43 - t82;
t79 = -Ifges(3,1) / 0.2e1;
t78 = -Ifges(5,6) * t64 / 0.2e1;
t77 = -Ifges(3,4) * t63 / 0.2e1;
t76 = t63 / 0.2e1;
t68 = pkin(3) + pkin(5);
t75 = t68 * t47;
t74 = -qJD(1) / 0.2e1;
t73 = qJD(1) / 0.2e1;
t72 = -qJD(2) / 0.2e1;
t71 = qJD(2) / 0.2e1;
t27 = -qJD(2) * pkin(2) + qJD(3) + t41;
t70 = m(4) * t27 + (mrSges(4,1) + mrSges(3,3)) * t64 + t83 * qJD(2);
t60 = qJD(2) * qJ(3);
t30 = -t43 - t60;
t33 = -mrSges(4,1) * t63 - qJD(2) * mrSges(4,3);
t69 = m(4) * t30 + qJD(2) * mrSges(3,2) - mrSges(3,3) * t63 + t33;
t67 = (pkin(1) * mrSges(3,1));
t66 = pkin(1) * mrSges(3,2);
t46 = -pkin(2) - qJ(4);
t34 = mrSges(5,1) * t63 + qJD(2) * mrSges(5,2);
t65 = t33 - t34;
t62 = qJD(2) * t47;
t61 = t47 * qJD(3);
t42 = pkin(2) * t64;
t37 = t68 * t48;
t58 = m(4) * pkin(5) + mrSges(4,1);
t57 = -t47 * qJ(3) - pkin(1);
t56 = Ifges(3,5) / 0.2e1 + Ifges(5,5) / 0.2e1 - Ifges(4,4) / 0.2e1;
t55 = -Ifges(3,6) / 0.2e1 + Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t54 = 0.3e1 / 0.2e1 * Ifges(4,6) + 0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * Ifges(5,6);
t53 = -qJ(3) * t48 + qJ(4) * t47;
t28 = -pkin(2) * t48 + t57;
t52 = -t48 * t60 - t61;
t18 = t46 * t48 + t57;
t51 = qJD(2) * t53 - qJD(4) * t48 - t61;
t17 = t28 * qJD(1);
t3 = t18 * qJD(1);
t5 = qJD(2) * t46 + t81;
t50 = t17 * mrSges(4,3) + t3 * mrSges(5,2) + t64 * t79 + t77 + (-Ifges(5,6) * t48 + Ifges(5,3) * t47) * t74 + Ifges(4,4) * t71 + (-Ifges(4,2) * t47 - Ifges(4,6) * t48) * t73 - t27 * mrSges(4,1) - t5 * mrSges(5,1) + (Ifges(3,5) + Ifges(5,5)) * t72;
t10 = -t30 + t82;
t49 = t10 * mrSges(5,1) + t17 * mrSges(4,2) + Ifges(3,6) * t71 + (Ifges(3,4) * t47 + Ifges(3,2) * t48) * t73 + (-Ifges(4,6) * t47 - Ifges(4,3) * t48) * t74 + Ifges(5,2) * t76 + t78 - t3 * mrSges(5,3) - t30 * mrSges(4,1) + (Ifges(4,5) + Ifges(5,4)) * t72;
t45 = pkin(2) * t62;
t38 = qJD(2) * t42;
t32 = mrSges(5,1) * t64 - qJD(2) * mrSges(5,3);
t26 = (-qJD(3) + t41) * qJD(2);
t25 = qJD(2) * t37;
t24 = qJD(2) * t75;
t21 = (-t47 * mrSges(5,2) - t48 * mrSges(5,3)) * qJD(1);
t19 = (t48 * mrSges(4,2) - t47 * mrSges(4,3)) * qJD(1);
t9 = t45 + t52;
t8 = qJD(1) * t53 + t42;
t7 = (qJD(1) * t37 - qJD(4)) * qJD(2);
t6 = (-qJD(1) * t75 + qJD(3)) * qJD(2);
t4 = qJD(1) * t52 + t38;
t2 = t45 + t51;
t1 = qJD(1) * t51 + t38;
t11 = [t9 * t19 + t2 * t21 - t24 * t34 + t25 * t32 + m(4) * (t17 * t9 + t4 * t28) + m(5) * (t1 * t18 - t10 * t24 + t2 * t3 + t25 * t5 + t37 * t6 + t7 * t75) + (t7 * mrSges(5,1) - t1 * mrSges(5,2) - t4 * mrSges(4,3) + (t55 * qJD(2) + t69 * pkin(5) + (-t37 * mrSges(5,1) - t28 * mrSges(4,2) + t18 * mrSges(5,3) - t47 * t54 - (2 * t67)) * qJD(1) - t49) * qJD(2)) * t47 + (t6 * mrSges(5,1) + t4 * mrSges(4,2) - t1 * mrSges(5,3) - t58 * t26 + (t56 * qJD(2) + t70 * pkin(5) + (t54 * t48 - 0.2e1 * t66 - t18 * mrSges(5,2) - t28 * mrSges(4,3) + t75 * mrSges(5,1) + (0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(3,1) + 0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(5,2) + t58 * pkin(5)) * t47) * qJD(1) - t50) * qJD(2)) * t48; t6 * mrSges(5,2) - t26 * mrSges(4,3) - t7 * mrSges(5,3) - t8 * t21 - t22 * t34 + t80 * t32 - t65 * qJD(3) + m(4) * (-t26 * qJ(3) - t30 * qJD(3)) + ((t77 + (t66 + (Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * t48) * qJD(1) + t50) * t48 + (t78 + (t67 + (Ifges(3,4) / 0.2e1 + Ifges(4,6) / 0.2e1) * t47) * qJD(1) + (-Ifges(5,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + t79 + Ifges(5,2) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t63 + t49) * t47 + (-t69 * t47 - t70 * t48) * pkin(5) + ((pkin(5) * mrSges(3,2) + (-mrSges(4,1) - mrSges(5,1)) * qJ(3) + t55) * t47 + (-pkin(2) * mrSges(4,1) + t46 * mrSges(5,1) + (-m(4) * pkin(2) + t83) * pkin(5) + t56) * t48) * qJD(2)) * qJD(1) + (-m(4) * t17 - t19) * (-qJ(3) * t63 + t42) + (t6 * qJ(3) + t81 * t10 - t3 * t8 + t7 * t46 + t80 * t5) * m(5); t65 * qJD(2) + ((t19 + t21) * t47 + (mrSges(5,1) + t58) * t48 * qJD(2)) * qJD(1) - m(4) * (-qJD(2) * t30 - t17 * t64) + (-t10 * qJD(2) + t3 * t64 + t7) * m(5); (-mrSges(5,1) * t62 + t48 * t21) * qJD(1) + qJD(2) * t32 + (t5 * qJD(2) + 0.2e1 * t3 * t76 + t6) * m(5);];
tauc = t11(:);
