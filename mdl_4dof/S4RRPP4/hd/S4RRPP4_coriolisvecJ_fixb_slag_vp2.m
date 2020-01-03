% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPP4
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:53
% EndTime: 2019-12-31 16:58:55
% DurationCPUTime: 0.80s
% Computational Cost: add. (390->160), mult. (1092->207), div. (0->0), fcn. (399->2), ass. (0->64)
t83 = Ifges(5,4) + Ifges(4,5);
t81 = -mrSges(3,1) - mrSges(4,1);
t78 = -qJD(2) / 0.2e1;
t77 = qJD(2) / 0.2e1;
t46 = sin(qJ(2));
t65 = qJD(1) * t46;
t42 = pkin(5) * t65;
t60 = qJ(4) * qJD(1);
t21 = t46 * t60 - t42;
t76 = -t21 + qJD(3);
t75 = -t83 * t65 / 0.2e1;
t47 = cos(qJ(2));
t64 = qJD(1) * t47;
t43 = pkin(5) * t64;
t25 = -t47 * t60 + t43;
t45 = qJD(2) * qJ(3);
t12 = t25 + t45;
t58 = qJ(3) * t46 + pkin(1);
t28 = -pkin(2) * t47 - t58;
t19 = t28 * qJD(1);
t32 = t43 + t45;
t35 = mrSges(4,2) * t64 + qJD(2) * mrSges(4,3);
t72 = pkin(2) + pkin(3);
t20 = t72 * t47 + t58;
t4 = t20 * qJD(1) + qJD(4);
t69 = Ifges(3,4) * t46;
t74 = -(m(4) * t32 - qJD(2) * mrSges(3,2) + mrSges(3,3) * t64 + t35) * pkin(5) + t12 * mrSges(5,3) + t19 * mrSges(4,1) + Ifges(4,6) * t77 - (Ifges(3,2) * t47 + t69) * qJD(1) / 0.2e1 - t32 * mrSges(4,2) - t4 * mrSges(5,1) - t75 - (Ifges(4,3) + Ifges(5,2)) * t64 / 0.2e1 + (Ifges(5,6) + Ifges(3,6)) * t78;
t27 = -qJD(2) * pkin(2) + qJD(3) + t42;
t41 = Ifges(3,4) * t64;
t7 = -t72 * qJD(2) + t76;
t73 = (m(4) * t27 + (mrSges(4,2) + mrSges(3,3)) * t65 + t81 * qJD(2)) * pkin(5) + t27 * mrSges(4,2) + t4 * mrSges(5,2) + Ifges(5,5) * t78 + Ifges(3,1) * t65 / 0.2e1 + t41 / 0.2e1 - t19 * mrSges(4,3) - t7 * mrSges(5,3) + (-t83 * t47 + (Ifges(4,1) + Ifges(5,1)) * t46) * qJD(1) / 0.2e1 + (Ifges(4,4) + Ifges(3,5)) * t77;
t71 = pkin(1) * mrSges(3,1);
t70 = pkin(1) * mrSges(3,2);
t68 = pkin(5) - qJ(4);
t33 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t64;
t67 = t33 + t35;
t66 = qJ(3) * t47;
t63 = qJD(2) * t46;
t62 = qJD(2) * t47;
t61 = qJD(3) * t46;
t59 = m(4) * pkin(5) + mrSges(4,2);
t37 = t68 * t47;
t55 = 0.3e1 / 0.2e1 * Ifges(5,4) + 0.3e1 / 0.2e1 * Ifges(4,5) - 0.3e1 / 0.2e1 * Ifges(3,4);
t54 = -Ifges(5,5) / 0.2e1 + Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t53 = -Ifges(5,6) / 0.2e1 + Ifges(4,6) / 0.2e1 - Ifges(3,6) / 0.2e1;
t52 = pkin(2) * t46 - t66;
t51 = -t72 * t46 + t66;
t11 = qJD(2) * t37 - qJD(4) * t46;
t9 = -qJD(4) * t47 - t68 * t63;
t10 = t52 * qJD(2) - t61;
t3 = t51 * qJD(2) + t61;
t44 = qJD(2) * qJD(3);
t38 = qJD(1) * mrSges(5,2) * t62;
t36 = t68 * t46;
t29 = -qJD(2) * mrSges(5,1) - mrSges(5,3) * t65;
t26 = -qJD(2) * t42 + t44;
t23 = (t47 * mrSges(5,1) + t46 * mrSges(5,2)) * qJD(1);
t22 = (-t47 * mrSges(4,1) - t46 * mrSges(4,3)) * qJD(1);
t8 = t51 * qJD(1);
t6 = t11 * qJD(1);
t5 = t10 * qJD(1);
t2 = t9 * qJD(1) + t44;
t1 = t3 * qJD(1);
t13 = [t10 * t22 + t11 * t29 + t20 * t38 + t3 * t23 + t9 * t33 + m(4) * (t10 * t19 + t28 * t5) + m(5) * (t1 * t20 + t11 * t7 + t12 * t9 + t2 * t37 + t3 * t4 + t6 * t36) + (t1 * mrSges(5,2) - t5 * mrSges(4,3) - t6 * mrSges(5,3) + (t53 * qJD(2) + (t28 * mrSges(4,1) - t20 * mrSges(5,1) + t37 * mrSges(5,3) + t55 * t46 - 0.2e1 * t71) * qJD(1) + t74) * qJD(2)) * t46 + (-t5 * mrSges(4,1) + t1 * mrSges(5,1) - t2 * mrSges(5,3) + t59 * t26 + (t54 * qJD(2) + (-t36 * mrSges(5,3) - 0.2e1 * t70 - t28 * mrSges(4,3) - t55 * t47 + (0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + t59 * pkin(5)) * t46) * qJD(1) + t73) * qJD(2)) * t47; -t6 * mrSges(5,1) + t2 * mrSges(5,2) + t26 * mrSges(4,3) - t21 * t33 - t8 * t23 - t25 * t29 + t67 * qJD(3) + m(4) * (t26 * qJ(3) + t32 * qJD(3)) + (((t71 + t69 / 0.2e1) * qJD(1) + (pkin(5) * mrSges(3,2) + (-mrSges(4,2) + mrSges(5,3)) * qJ(3) + t53) * qJD(2) - t74 + t75) * t46 + (-t41 / 0.2e1 + (t70 + (Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t47) * qJD(1) + (-Ifges(4,1) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1) * t65 + (-pkin(2) * mrSges(4,2) + t72 * mrSges(5,3) + (-m(4) * pkin(2) + t81) * pkin(5) + t54) * qJD(2) - t73) * t47) * qJD(1) + (-m(4) * t19 - t22) * t52 * qJD(1) + (t2 * qJ(3) + t76 * t12 - t7 * t25 - t4 * t8 - t6 * t72) * m(5); -t67 * qJD(2) + ((t22 - t23) * t46 + (-mrSges(5,3) + t59) * t62) * qJD(1) - m(4) * (qJD(2) * t32 - t19 * t65) + (-t12 * qJD(2) - t4 * t65 + t6) * m(5); t38 + m(5) * t1 + (-mrSges(5,1) * t63 + t46 * t29 + t47 * t33 - m(5) * (-t12 * t47 - t46 * t7)) * qJD(1);];
tauc = t13(:);
