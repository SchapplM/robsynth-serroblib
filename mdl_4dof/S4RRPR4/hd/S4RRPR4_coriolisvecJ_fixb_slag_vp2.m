% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:26
% EndTime: 2019-12-31 17:02:28
% DurationCPUTime: 0.67s
% Computational Cost: add. (946->132), mult. (1791->200), div. (0->0), fcn. (1043->6), ass. (0->74)
t72 = cos(qJ(2));
t90 = pkin(1) * qJD(1);
t84 = t72 * t90;
t75 = qJD(3) - t84;
t67 = sin(pkin(7));
t68 = cos(pkin(7));
t69 = sin(qJ(4));
t71 = cos(qJ(4));
t52 = -t69 * t67 + t71 * t68;
t47 = t52 * qJD(4);
t106 = t47 / 0.2e1;
t53 = t71 * t67 + t69 * t68;
t48 = t53 * qJD(4);
t105 = -t48 / 0.2e1;
t66 = qJD(1) + qJD(2);
t89 = pkin(1) * qJD(2);
t82 = qJD(1) * t89;
t76 = t72 * t82;
t51 = t66 * qJD(3) + t76;
t91 = t67 ^ 2 + t68 ^ 2;
t80 = t91 * t51;
t104 = t91 * mrSges(4,3) * t66;
t57 = (-pkin(6) - qJ(3)) * t67;
t62 = t68 * pkin(6);
t58 = t68 * qJ(3) + t62;
t28 = t71 * t57 - t69 * t58;
t102 = t28 * qJD(4) + t75 * t52;
t29 = t69 * t57 + t71 * t58;
t101 = -t29 * qJD(4) - t75 * t53;
t41 = t52 * t66;
t42 = t53 * t66;
t74 = -t68 * mrSges(4,1) + t67 * mrSges(4,2);
t100 = -t41 * mrSges(5,1) + t42 * mrSges(5,2) + t74 * t66;
t70 = sin(qJ(2));
t85 = t70 * t90;
t55 = t66 * qJ(3) + t85;
t81 = pkin(6) * t66 + t55;
t33 = t81 * t67;
t34 = t81 * t68;
t16 = -t71 * t33 - t69 * t34;
t17 = -t69 * t33 + t71 * t34;
t3 = -t17 * qJD(4) - t53 * t51;
t99 = -t16 * t47 - t3 * t53;
t97 = t42 / 0.2e1;
t95 = t72 * pkin(1);
t93 = Ifges(5,4) * t42;
t88 = qJD(2) * t72;
t83 = mrSges(3,2) * t88;
t61 = -t68 * pkin(3) - pkin(2);
t79 = t91 * t55;
t37 = t66 * t47;
t38 = t66 * t48;
t18 = t38 * mrSges(5,1) + t37 * mrSges(5,2);
t78 = t91 * qJD(3);
t77 = t70 * t82;
t60 = t70 * pkin(1) + qJ(3);
t49 = (-pkin(6) - t60) * t67;
t50 = t68 * t60 + t62;
t26 = t71 * t49 - t69 * t50;
t27 = t69 * t49 + t71 * t50;
t19 = Ifges(5,2) * t41 + Ifges(5,6) * qJD(4) + t93;
t2 = t16 * qJD(4) + t52 * t51;
t39 = Ifges(5,4) * t41;
t20 = Ifges(5,1) * t42 + Ifges(5,5) * qJD(4) + t39;
t40 = t61 * t66 + t75;
t73 = t20 * t106 + t40 * (t48 * mrSges(5,1) + t47 * mrSges(5,2)) + qJD(4) * (Ifges(5,5) * t47 - Ifges(5,6) * t48) / 0.2e1 + t19 * t105 + mrSges(4,3) * t80 + (-t52 * mrSges(5,1) + t53 * mrSges(5,2) + t74) * t77 + (t41 * t105 - t38 * t52) * Ifges(5,2) + (t53 * t37 + t47 * t97) * Ifges(5,1) + (-t17 * t48 + t2 * t52) * mrSges(5,3) + (t41 * t106 + t37 * t52 - t53 * t38 - t48 * t97) * Ifges(5,4);
t59 = pkin(1) * t88 + qJD(3);
t56 = t61 - t95;
t54 = -t66 * pkin(2) + t75;
t31 = qJD(4) * mrSges(5,1) - t42 * mrSges(5,3);
t30 = -qJD(4) * mrSges(5,2) + t41 * mrSges(5,3);
t6 = -t27 * qJD(4) - t53 * t59;
t5 = t26 * qJD(4) + t52 * t59;
t1 = [m(4) * (t59 * t79 + t60 * t80) + m(5) * (t16 * t6 + t17 * t5 + t2 * t27 + t3 * t26) + t73 + t56 * t18 + t5 * t30 + t6 * t31 - mrSges(3,1) * t77 - mrSges(3,2) * t76 - pkin(1) * t66 * t83 + t59 * t104 + (m(4) * (t54 + (-pkin(2) - t95) * qJD(1)) + m(5) * (qJD(1) * t56 + t40) - t66 * mrSges(3,1) + t100) * t70 * t89 + (-t26 * t37 - t27 * t38 + t99) * mrSges(5,3); t73 + ((mrSges(3,1) * t70 + mrSges(3,2) * t72) * t90 + (-t91 * t84 + t78) * mrSges(4,3)) * t66 + t101 * t31 + t102 * t30 + (-t83 + (-mrSges(3,1) * qJD(2) - t100) * t70) * t90 + (-t28 * t37 - t29 * t38 + t99) * mrSges(5,3) + t61 * t18 + (t101 * t16 + t102 * t17 + t2 * t29 + t3 * t28 - t40 * t85 + t61 * t77) * m(5) + (-pkin(2) * t77 + qJ(3) * t80 + t55 * t78 - (t54 * t70 + t72 * t79) * t90) * m(4); -t41 * t30 + t42 * t31 + (m(4) + m(5)) * t77 - m(5) * (-t16 * t42 + t17 * t41) + t18 + (-m(4) * t79 - t104) * t66; Ifges(5,5) * t37 - Ifges(5,6) * t38 - t2 * mrSges(5,2) + t3 * mrSges(5,1) - t40 * (t42 * mrSges(5,1) + t41 * mrSges(5,2)) - t42 * (Ifges(5,1) * t41 - t93) / 0.2e1 + t19 * t97 - qJD(4) * (Ifges(5,5) * t41 - Ifges(5,6) * t42) / 0.2e1 - t16 * t30 + t17 * t31 + (t16 * t41 + t17 * t42) * mrSges(5,3) - (-Ifges(5,2) * t42 + t20 + t39) * t41 / 0.2e1;];
tauc = t1(:);
