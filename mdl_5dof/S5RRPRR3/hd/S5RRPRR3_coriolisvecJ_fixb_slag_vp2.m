% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:21
% EndTime: 2020-01-03 12:00:23
% DurationCPUTime: 0.82s
% Computational Cost: add. (2004->157), mult. (4024->236), div. (0->0), fcn. (2242->8), ass. (0->89)
t68 = cos(pkin(9));
t71 = sin(qJ(2));
t111 = t68 * t71;
t67 = sin(pkin(9));
t74 = cos(qJ(2));
t85 = pkin(1) * (-t67 * t74 - t111);
t49 = qJD(1) * t85;
t112 = t67 * t71;
t84 = pkin(1) * (t68 * t74 - t112);
t51 = qJD(1) * t84;
t70 = sin(qJ(4));
t73 = cos(qJ(4));
t125 = pkin(2) * t67;
t63 = pkin(2) * t68 + pkin(3);
t86 = -t125 * t70 + t63 * t73;
t109 = -t86 * qJD(4) + t49 * t70 + t51 * t73;
t105 = qJD(1) * pkin(1);
t66 = qJD(1) + qJD(2);
t58 = t66 * pkin(2) + t105 * t74;
t97 = t71 * t105;
t26 = t68 * t58 - t67 * t97;
t21 = pkin(3) * t66 + t26;
t27 = t58 * t67 + t68 * t97;
t16 = t21 * t70 + t27 * t73;
t65 = qJD(4) + t66;
t14 = pkin(8) * t65 + t16;
t69 = sin(qJ(5));
t72 = cos(qJ(5));
t11 = qJD(3) * t72 - t14 * t69;
t12 = qJD(3) * t69 + t14 * t72;
t89 = -t11 * t69 + t12 * t72;
t127 = m(6) * t89;
t64 = pkin(1) * t74 + pkin(2);
t93 = -pkin(1) * t112 + t68 * t64;
t46 = pkin(3) + t93;
t53 = pkin(1) * t111 + t64 * t67;
t108 = t70 * t46 + t73 * t53;
t106 = t73 * t125 + t70 * t63;
t107 = t106 * qJD(4) + t73 * t49 - t51 * t70;
t126 = mrSges(3,1) * t71 + mrSges(3,2) * t74;
t15 = t21 * t73 - t27 * t70;
t50 = qJD(2) * t85;
t38 = qJD(1) * t50;
t52 = qJD(2) * t84;
t39 = qJD(1) * t52;
t6 = qJD(4) * t15 + t38 * t70 + t39 * t73;
t2 = qJD(5) * t11 + t6 * t72;
t124 = t2 * t72;
t102 = qJD(5) * t12;
t3 = -t6 * t69 - t102;
t123 = t3 * t69;
t120 = Ifges(6,1) * t69;
t119 = Ifges(6,4) * t69;
t117 = t11 * mrSges(6,3);
t114 = t65 * t72;
t113 = t66 * mrSges(4,1);
t110 = t69 * mrSges(6,3);
t104 = Ifges(6,5) * qJD(5);
t103 = Ifges(6,6) * qJD(5);
t101 = qJD(5) * t69;
t100 = qJD(5) * t72;
t95 = t100 / 0.2e1;
t91 = -mrSges(6,1) * t72 + mrSges(6,2) * t69;
t45 = t91 * t65;
t94 = t65 * mrSges(5,1) - t45;
t90 = -t11 * t72 - t12 * t69;
t88 = t46 * t73 - t53 * t70;
t56 = qJD(5) * mrSges(6,1) - t110 * t65;
t57 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t114;
t87 = -t69 * t56 + t72 * t57;
t83 = (Ifges(6,2) * t72 + t119) * t65;
t82 = (mrSges(6,1) * t69 + mrSges(6,2) * t72) * qJD(5);
t81 = t65 * mrSges(5,2) - t87;
t80 = qJD(5) * t90 - t123;
t77 = m(6) * (t80 + t124);
t13 = -pkin(4) * t65 - t15;
t29 = t83 + t103;
t61 = Ifges(6,4) * t114;
t30 = t120 * t65 + t104 + t61;
t7 = qJD(4) * t16 - t73 * t38 + t39 * t70;
t76 = -t6 * mrSges(5,2) + mrSges(6,3) * t124 + t13 * t82 + t30 * t95 + qJD(5) ^ 2 * (Ifges(6,5) * t72 - Ifges(6,6) * t69) / 0.2e1 - (t29 + t83) * t101 / 0.2e1 + (-mrSges(5,1) + t91) * t7 + ((Ifges(6,1) * t72 - t119) * t101 + (0.3e1 * Ifges(6,4) * t72 - 0.2e1 * Ifges(6,2) * t69 + t120) * t95) * t65;
t75 = t38 * mrSges(4,1) + t76;
t48 = pkin(8) + t106;
t47 = -pkin(4) - t86;
t32 = t65 * t82;
t18 = pkin(8) + t108;
t17 = -pkin(4) - t88;
t9 = t108 * qJD(4) - t73 * t50 + t52 * t70;
t1 = [t18 * t77 + m(6) * (t13 * t9 + t17 * t7) + ((-t72 * t56 - t69 * t57) * t18 + t90 * mrSges(6,3)) * qJD(5) + t75 - t94 * t9 + m(4) * (t26 * t50 + t27 * t52 + t38 * t93 + t39 * t53) + m(5) * (t108 * t6 - t15 * t9 - t7 * t88) + (-t52 * t66 - t39) * mrSges(4,2) - t3 * t110 + t50 * t113 + t17 * t32 + t126 * pkin(1) * qJD(2) * (-qJD(1) - t66) + (m(5) * t16 + t127 - t81) * (qJD(4) * t88 + t50 * t70 + t52 * t73); (t51 * t66 - t39) * mrSges(4,2) + t75 + (-qJD(5) * t48 * t57 + t109 * t56 + (-t3 - t102) * mrSges(6,3)) * t69 + (-mrSges(5,1) * t107 + mrSges(5,2) * t109) * t65 + (-t109 * t57 + (-t48 * t56 - t117) * qJD(5)) * t72 - t49 * t113 + t47 * t32 + t107 * t45 + t48 * t77 + t126 * t105 * (-qJD(2) + t66) + (t107 * t13 - t109 * t89 + t47 * t7) * m(6) + (t106 * t6 - t107 * t15 - t109 * t16 - t7 * t86) * m(5) + (-t26 * t49 - t27 * t51 + (t38 * t68 + t39 * t67) * pkin(2)) * m(4); m(6) * (t2 * t69 + t3 * t72) + (t127 + (-t69 ^ 2 - t72 ^ 2) * t65 * mrSges(6,3) + t87) * qJD(5); t76 + t80 * mrSges(6,3) + t94 * t16 + t81 * t15 - pkin(4) * t32 + (-t56 * t100 - t57 * t101) * pkin(8) + (-t7 * pkin(4) + (-t100 * t11 - t101 * t12 - t123 + t124) * pkin(8) - t13 * t16 - t89 * t15) * m(6); t3 * mrSges(6,1) - t2 * mrSges(6,2) - t11 * t57 + t12 * t56 + ((t104 / 0.2e1 - t13 * mrSges(6,2) - t30 / 0.2e1 - t61 / 0.2e1 + t117) * t72 + (-t103 / 0.2e1 - t13 * mrSges(6,1) + t29 / 0.2e1 + t12 * mrSges(6,3) + (t119 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t72) * t65) * t69) * t65;];
tauc = t1(:);
