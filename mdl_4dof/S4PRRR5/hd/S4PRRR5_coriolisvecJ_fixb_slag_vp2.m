% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:35
% EndTime: 2019-12-31 16:33:36
% DurationCPUTime: 0.43s
% Computational Cost: add. (557->95), mult. (1249->156), div. (0->0), fcn. (726->6), ass. (0->59)
t41 = sin(qJ(4));
t44 = cos(qJ(4));
t86 = mrSges(5,1) * t44 - mrSges(5,2) * t41 + mrSges(4,1);
t46 = cos(qJ(2));
t34 = qJD(2) * pkin(2) + qJD(1) * t46;
t42 = sin(qJ(3));
t45 = cos(qJ(3));
t43 = sin(qJ(2));
t71 = qJD(1) * t43;
t20 = t34 * t42 + t45 * t71;
t38 = qJD(2) + qJD(3);
t14 = pkin(6) * t38 + t20;
t30 = t42 * t43 - t45 * t46;
t52 = t30 * qJD(2);
t69 = qJD(3) * t45;
t70 = qJD(3) * t42;
t6 = t34 * t69 + (-t43 * t70 - t52) * qJD(1);
t67 = qJD(4) * t44;
t3 = -t14 * t67 - t41 * t6;
t68 = qJD(4) * t41;
t2 = -t14 * t68 + t44 * t6;
t83 = t2 * t44;
t61 = -t3 * t41 + t83;
t66 = t14 * (t41 ^ 2 + t44 ^ 2);
t85 = t41 / 0.2e1;
t31 = t42 * t46 + t43 * t45;
t53 = t31 * qJD(2);
t7 = t34 * t70 + (t43 * t69 + t53) * qJD(1);
t81 = t30 * t7;
t80 = Ifges(5,1) * t41;
t79 = Ifges(5,4) * t41;
t77 = Ifges(5,5) * t44;
t76 = Ifges(5,2) * t44;
t75 = Ifges(5,6) * t41;
t73 = t38 * t44;
t72 = t41 * mrSges(5,3);
t64 = t67 / 0.2e1;
t63 = t86 * t38;
t59 = mrSges(5,1) * t41 + mrSges(5,2) * t44;
t32 = qJD(4) * mrSges(5,1) - t38 * t72;
t33 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t73;
t58 = t32 * t44 + t33 * t41;
t19 = t34 * t45 - t42 * t71;
t57 = t41 * (Ifges(5,1) * t44 - t79);
t56 = (t76 + t79) * t38;
t55 = t59 * qJD(4);
t54 = t58 * qJD(4);
t51 = t38 * mrSges(4,2) + t41 * t32 - t44 * t33;
t13 = -pkin(3) * t38 - t19;
t21 = Ifges(5,6) * qJD(4) + t56;
t35 = Ifges(5,4) * t73;
t22 = Ifges(5,5) * qJD(4) + t38 * t80 + t35;
t48 = -t6 * mrSges(4,2) + mrSges(5,3) * t83 + t13 * t55 + t22 * t64 - t3 * t72 + qJD(4) ^ 2 * (-t75 + t77) / 0.2e1 - t86 * t7 - (t21 + t56) * t68 / 0.2e1 + (qJD(4) * t57 + (0.3e1 * Ifges(5,4) * t44 - 0.2e1 * Ifges(5,2) * t41 + t80) * t64) * t38;
t29 = t30 * qJD(1);
t28 = t31 * qJD(1);
t23 = t38 * t55;
t10 = t31 * qJD(3) + t53;
t9 = -t30 * qJD(3) - t52;
t1 = [t30 * t23 + (-mrSges(3,1) * t43 - mrSges(3,2) * t46) * qJD(2) ^ 2 - t63 * t10 - t31 * t54 - t51 * t9 + m(4) * (-t10 * t19 + t20 * t9 + t31 * t6 + t81) + m(5) * (t10 * t13 + t61 * t31 + t9 * t66 + t81); t48 - t51 * t29 + t63 * t28 - m(4) * (-t19 * t28 - t20 * t29) - m(5) * (t13 * t28 - t29 * t66) + (m(4) * (t42 * t6 - t45 * t7) + (-t63 * t42 - t51 * t45 + m(4) * (-t19 * t42 + t20 * t45) + m(5) * (t13 * t42 + t45 * t66)) * qJD(3)) * pkin(2) + (m(5) * t7 + t23) * (-pkin(2) * t45 - pkin(3)) + (m(5) * t61 - t54) * (pkin(2) * t42 + pkin(6)); -pkin(3) * t23 - pkin(6) * t54 + t51 * t19 + t63 * t20 + t48 + (-t7 * pkin(3) + t61 * pkin(6) - t13 * t20 - t19 * t66) * m(5); t3 * mrSges(5,1) - t2 * mrSges(5,2) + t58 * t14 + (-t13 * t59 + t21 * t85 + (-t57 / 0.2e1 + t76 * t85) * t38 + (t77 / 0.2e1 - t75 / 0.2e1) * qJD(4) - (t22 + t35) * t44 / 0.2e1) * t38;];
tauc = t1(:);
