% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRR4
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:30
% EndTime: 2019-12-31 16:32:32
% DurationCPUTime: 0.75s
% Computational Cost: add. (731->151), mult. (1953->231), div. (0->0), fcn. (1147->4), ass. (0->79)
t66 = Ifges(4,5) * qJD(3) / 0.2e1;
t57 = qJD(3) + qJD(4);
t58 = sin(qJ(4));
t59 = sin(qJ(3));
t60 = cos(qJ(4));
t61 = cos(qJ(3));
t42 = -t58 * t59 + t60 * t61;
t38 = t42 * qJD(2);
t91 = t38 / 0.2e1;
t43 = t58 * t61 + t60 * t59;
t39 = t43 * qJD(2);
t90 = t39 / 0.2e1;
t89 = -pkin(6) - pkin(5);
t54 = -pkin(3) * t61 - pkin(2);
t50 = qJD(2) * t54;
t88 = m(5) * t50;
t87 = pkin(2) * mrSges(4,1);
t86 = pkin(2) * mrSges(4,2);
t85 = mrSges(5,3) * t38;
t84 = Ifges(4,4) * t59;
t83 = Ifges(5,4) * t39;
t52 = t89 * t61;
t71 = t59 * qJD(1);
t35 = -qJD(2) * t52 + t71;
t82 = t35 * t58;
t81 = t35 * t60;
t80 = t39 * mrSges(5,3);
t24 = t57 * t42;
t20 = t24 * qJD(2);
t79 = t42 * t20;
t25 = t57 * t43;
t21 = t25 * qJD(2);
t78 = t43 * t21;
t76 = Ifges(4,6) * qJD(3);
t75 = qJD(2) * t59;
t74 = qJD(2) * t61;
t73 = qJD(4) * t58;
t72 = qJD(4) * t60;
t56 = t61 * qJD(1);
t70 = pkin(5) * t75;
t69 = m(4) * pkin(5) + mrSges(4,3);
t68 = qJD(2) * t89;
t67 = qJD(3) * t89;
t65 = -t76 / 0.2e1;
t46 = t56 - t70;
t55 = Ifges(4,4) * t74;
t64 = Ifges(4,1) * t75 / 0.2e1 + t55 / 0.2e1 + t66 - t46 * mrSges(4,3);
t47 = pkin(5) * t74 + t71;
t63 = t47 * mrSges(4,3) + t76 / 0.2e1 + (Ifges(4,2) * t61 + t84) * qJD(2) / 0.2e1;
t44 = t59 * t67;
t34 = t59 * t68 + t56;
t32 = qJD(3) * pkin(3) + t34;
t11 = t32 * t60 - t82;
t12 = t32 * t58 + t81;
t51 = t89 * t59;
t26 = t51 * t60 + t52 * t58;
t27 = t51 * t58 - t52 * t60;
t16 = Ifges(5,2) * t38 + Ifges(5,6) * t57 + t83;
t33 = Ifges(5,4) * t38;
t17 = Ifges(5,1) * t39 + Ifges(5,5) * t57 + t33;
t53 = qJD(3) * t56;
t30 = qJD(2) * t44 + t53;
t31 = (t61 * t68 - t71) * qJD(3);
t2 = qJD(4) * t11 + t30 * t60 + t31 * t58;
t3 = -qJD(4) * t12 - t30 * t58 + t31 * t60;
t62 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + t11 * t85 + t16 * t90 - t50 * (mrSges(5,1) * t39 + mrSges(5,2) * t38) - t39 * (Ifges(5,1) * t38 - t83) / 0.2e1 - t57 * (Ifges(5,5) * t38 - Ifges(5,6) * t39) / 0.2e1 - Ifges(5,6) * t21 + Ifges(5,5) * t20 - (-Ifges(5,2) * t39 + t17 + t33) * t38 / 0.2e1;
t49 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t74;
t48 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t75;
t45 = t61 * t67;
t41 = qJD(3) * t47;
t40 = -qJD(3) * t70 + t53;
t29 = mrSges(5,1) * t57 - t80;
t28 = -mrSges(5,2) * t57 + t85;
t23 = -mrSges(5,1) * t38 + mrSges(5,2) * t39;
t14 = t34 * t60 - t82;
t13 = -t34 * t58 - t81;
t6 = -t27 * qJD(4) - t44 * t58 + t45 * t60;
t5 = t26 * qJD(4) + t44 * t60 + t45 * t58;
t1 = [t24 * t28 - t25 * t29 + (-t78 - t79) * mrSges(5,3) + m(4) * (t40 * t59 - t41 * t61) + m(5) * (-t11 * t25 + t12 * t24 + t2 * t43 + t3 * t42) + (t61 * t49 - t59 * t48 + m(4) * (-t46 * t59 + t47 * t61) + (-t59 ^ 2 - t61 ^ 2) * qJD(2) * mrSges(4,3)) * qJD(3); m(5) * (t11 * t6 + t12 * t5 + t2 * t27 + t26 * t3) + t50 * (t25 * mrSges(5,1) + t24 * mrSges(5,2)) + t54 * (t21 * mrSges(5,1) + t20 * mrSges(5,2)) + t57 * (Ifges(5,5) * t24 - Ifges(5,6) * t25) / 0.2e1 - t25 * t16 / 0.2e1 + t5 * t28 + t6 * t29 + t24 * t17 / 0.2e1 + (-t42 * t21 - t25 * t91) * Ifges(5,2) + (t20 * t43 + t24 * t90) * Ifges(5,1) + (-t11 * t24 - t12 * t25 + t2 * t42 - t20 * t26 - t21 * t27 - t3 * t43) * mrSges(5,3) + (t69 * t40 + (t66 + (-0.2e1 * t86 + 0.3e1 / 0.2e1 * Ifges(4,4) * t61) * qJD(2) + (-m(4) * t46 - t48) * pkin(5) + t64) * qJD(3)) * t61 + (t24 * t91 - t25 * t90 - t78 + t79) * Ifges(5,4) + (t69 * t41 + (t65 + (-m(4) * t47 - t49) * pkin(5) + (-0.2e1 * t87 - 0.3e1 / 0.2e1 * t84 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t61) * qJD(2) + (0.2e1 * t88 + t23 + qJD(2) * (-t42 * mrSges(5,1) + t43 * mrSges(5,2))) * pkin(3) - t63) * qJD(3)) * t59; -m(5) * (t11 * t13 + t12 * t14) + t12 * t80 + t47 * t48 - t46 * t49 - t40 * mrSges(4,2) - t41 * mrSges(4,1) - t14 * t28 - t13 * t29 + t62 + (m(5) * (-t11 * t73 + t12 * t72 + t2 * t58 + t3 * t60) + t28 * t72 - t29 * t73 + (-t20 * t60 - t21 * t58) * mrSges(5,3)) * pkin(3) + ((-t55 / 0.2e1 + t66 + qJD(2) * t86 - t64) * t61 + (t65 + (t87 + t84 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t61) * qJD(2) + (-t23 - t88) * pkin(3) + t63) * t59) * qJD(2); t62 - t11 * t28 + (t29 + t80) * t12;];
tauc = t1(:);
