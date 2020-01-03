% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:07
% EndTime: 2020-01-03 11:45:11
% DurationCPUTime: 1.38s
% Computational Cost: add. (1239->204), mult. (2587->269), div. (0->0), fcn. (1256->6), ass. (0->98)
t154 = Ifges(5,1) + Ifges(6,1);
t90 = sin(qJ(4));
t92 = cos(qJ(4));
t153 = -mrSges(5,1) * t92 + mrSges(5,2) * t90 - mrSges(4,1);
t148 = pkin(1) * sin(pkin(8));
t79 = cos(pkin(8)) * pkin(1) + pkin(2);
t91 = sin(qJ(3));
t93 = cos(qJ(3));
t127 = t93 * t148 + t91 * t79;
t114 = -t91 * t148 + t79 * t93;
t120 = qJD(1) * t148;
t68 = t79 * qJD(1);
t39 = t120 * t93 + t91 * t68;
t85 = qJD(1) + qJD(3);
t22 = t85 * pkin(7) + t39;
t113 = qJ(5) * t85 + t22;
t124 = qJD(2) * t90;
t13 = t113 * t92 + t124;
t134 = t85 * t90;
t64 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t134;
t65 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t134;
t133 = t85 * t92;
t66 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t133;
t132 = t92 * mrSges(5,3);
t67 = -qJD(4) * mrSges(5,2) + t132 * t85;
t149 = (t66 + t67) * t92 - (t64 + t65) * t90;
t17 = t22 * t92 + t124;
t38 = -t91 * t120 + t68 * t93;
t31 = t38 * qJD(3);
t6 = -qJD(4) * t17 - t31 * t90;
t147 = t6 * t90;
t146 = t92 * pkin(4);
t145 = Ifges(5,4) * t90;
t143 = Ifges(6,4) * t90;
t83 = t92 * qJD(2);
t16 = -t22 * t90 + t83;
t141 = t16 * mrSges(5,3);
t140 = t16 * t90;
t139 = t17 * mrSges(5,3);
t44 = t114 * qJD(3);
t138 = t44 * t92;
t57 = pkin(7) + t127;
t137 = t57 * t92;
t135 = t85 * mrSges(4,2);
t131 = -qJ(5) - pkin(7);
t130 = qJD(4) * t83 + t92 * t31;
t123 = qJD(4) * t90;
t118 = t85 * t123;
t122 = qJD(4) * t92;
t50 = t85 * mrSges(6,2) * t122 + mrSges(6,1) * t118;
t125 = -qJ(5) - t57;
t82 = t92 * qJD(5);
t119 = pkin(4) * t123;
t81 = -pkin(3) - t146;
t115 = t153 * t85;
t112 = qJD(4) * t131;
t111 = qJD(4) * t125;
t56 = -pkin(3) - t114;
t107 = -mrSges(6,1) * t92 + mrSges(6,2) * t90;
t103 = t113 * t90;
t12 = t83 - t103;
t10 = qJD(4) * pkin(4) + t12;
t106 = -t10 * t90 + t13 * t92;
t105 = -t16 * t92 - t17 * t90;
t104 = t17 * t92 - t140;
t102 = m(6) * t106;
t101 = (Ifges(5,2) * t92 + t145) * t85;
t100 = (Ifges(6,2) * t92 + t143) * t85;
t99 = (mrSges(5,1) * t90 + mrSges(5,2) * t92) * qJD(4);
t45 = t127 * qJD(3);
t32 = t39 * qJD(3);
t15 = t81 * t85 + qJD(5) - t38;
t19 = pkin(4) * t118 + t32;
t2 = -qJD(4) * t103 + t82 * t85 + t130;
t21 = -t85 * pkin(3) - t38;
t3 = (-qJD(5) * t85 - t31) * t90 - t13 * qJD(4);
t46 = Ifges(6,6) * qJD(4) + t100;
t47 = Ifges(5,6) * qJD(4) + t101;
t75 = Ifges(6,4) * t133;
t48 = Ifges(6,1) * t134 + Ifges(6,5) * qJD(4) + t75;
t76 = Ifges(5,4) * t133;
t49 = Ifges(5,1) * t134 + Ifges(5,5) * qJD(4) + t76;
t5 = -t123 * t22 + t130;
t94 = t15 * (mrSges(6,1) * t90 + mrSges(6,2) * t92) * qJD(4) + t19 * t107 + (-t6 * mrSges(5,3) - t3 * mrSges(6,3)) * t90 + t5 * t132 + t21 * t99 + t2 * t92 * mrSges(6,3) - t31 * mrSges(4,2) + t153 * t32 + (t154 * t92 - t143 - t145) * t118 + ((Ifges(5,5) + Ifges(6,5)) * t92 + (-Ifges(5,6) - Ifges(6,6)) * t90) * qJD(4) ^ 2 / 0.2e1 - (t100 + t101 + t47 + t46) * t123 / 0.2e1 + (t49 + t48 + ((-0.2e1 * Ifges(5,2) - 0.2e1 * Ifges(6,2) + t154) * t90 + 0.3e1 * (Ifges(5,4) + Ifges(6,4)) * t92) * t85) * t122 / 0.2e1;
t84 = t92 * qJ(5);
t70 = pkin(7) * t92 + t84;
t69 = t131 * t90;
t60 = t107 * t85;
t59 = -qJD(5) * t90 + t112 * t92;
t58 = t112 * t90 + t82;
t51 = t85 * t99;
t37 = t56 - t146;
t26 = t45 + t119;
t24 = t84 + t137;
t23 = t125 * t90;
t8 = (-qJD(5) - t44) * t90 + t92 * t111;
t7 = t111 * t90 + t138 + t82;
t1 = [t115 * t45 + (-t90 * t65 + t92 * t67 - t135) * t44 + m(5) * (t137 * t5 + t138 * t17 - t140 * t44 - t147 * t57 + t21 * t45 + t32 * t56) + m(4) * (-t32 * t114 + t127 * t31 - t38 * t45 + t39 * t44) + m(6) * (t10 * t8 + t13 * t7 + t15 * t26 + t19 * t37 + t2 * t24 + t23 * t3) + t94 + t7 * t66 + t37 * t50 + t56 * t51 + t26 * t60 + t8 * t64 + (t105 * mrSges(5,3) + (m(5) * t105 - t92 * t65 - t90 * t67) * t57 + (-t10 * t92 - t13 * t90 + (-t23 * t92 - t24 * t90) * t85) * mrSges(6,3)) * qJD(4); m(5) * (t5 * t90 + t6 * t92) + m(6) * (t2 * t90 + t3 * t92) + (m(5) * t104 + t102 + (mrSges(5,3) + mrSges(6,3)) * t85 * (-t90 ^ 2 - t92 ^ 2) + t149) * qJD(4); ((-t141 - pkin(7) * t65 + (-t69 * t85 - t10) * mrSges(6,3)) * t92 + (-t139 + pkin(4) * t60 - pkin(7) * t67 + (-t70 * t85 - t13) * mrSges(6,3)) * t90) * qJD(4) + (t135 - t149) * t38 + (-t115 - t60) * t39 + t94 + t58 * t66 + t81 * t50 - pkin(3) * t51 + t59 * t64 + (t10 * t59 - t106 * t38 + t13 * t58 + t19 * t81 + t2 * t70 + t3 * t69 + (t119 - t39) * t15) * m(6) + ((qJD(4) * t105 + t5 * t92 - t147) * pkin(7) - t32 * pkin(3) - t104 * t38 - t21 * t39) * m(5); t6 * mrSges(5,1) - t5 * mrSges(5,2) - t2 * mrSges(6,2) - t12 * t66 - t16 * t67 + t17 * t65 + (m(6) * pkin(4) + mrSges(6,1)) * t3 + (t64 - m(6) * (-t10 + t12)) * t13 + ((-t75 / 0.2e1 - t76 / 0.2e1 - t48 / 0.2e1 - t49 / 0.2e1 + t141 + t10 * mrSges(6,3) - t15 * mrSges(6,2) - t21 * mrSges(5,2) + (Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1 - pkin(4) * mrSges(6,3)) * qJD(4)) * t92 + (t139 + t13 * mrSges(6,3) - t15 * mrSges(6,1) - t21 * mrSges(5,1) + t46 / 0.2e1 + t47 / 0.2e1 + (Ifges(5,4) / 0.2e1 + Ifges(6,4) / 0.2e1) * t134 + (-Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * qJD(4) + (-m(6) * t15 - t60) * pkin(4) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1 - Ifges(6,1) / 0.2e1) * t133) * t90) * t85; m(6) * t19 + (t90 * t64 - t92 * t66 - t102) * t85 + t50;];
tauc = t1(:);
