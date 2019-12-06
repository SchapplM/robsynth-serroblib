% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:40
% EndTime: 2019-12-05 17:07:44
% DurationCPUTime: 1.28s
% Computational Cost: add. (1987->215), mult. (3586->319), div. (0->0), fcn. (2063->6), ass. (0->126)
t104 = qJD(4) + qJD(5);
t106 = sin(qJ(5));
t107 = sin(qJ(4));
t109 = cos(qJ(5));
t110 = cos(qJ(4));
t84 = -t106 * t107 + t109 * t110;
t48 = t104 * t84;
t178 = t48 / 0.2e1;
t105 = qJD(2) + qJD(3);
t72 = t84 * t105;
t177 = -t72 / 0.2e1;
t36 = t48 * t105;
t176 = t36 * t84;
t85 = t106 * t110 + t107 * t109;
t49 = t104 * t85;
t37 = t49 * t105;
t175 = t37 * t85;
t111 = cos(qJ(3));
t157 = pkin(2) * qJD(2);
t137 = t111 * t157;
t167 = -pkin(8) - pkin(7);
t95 = t167 * t107;
t103 = t110 * pkin(8);
t96 = pkin(7) * t110 + t103;
t56 = -t106 * t96 + t109 * t95;
t133 = qJD(4) * t167;
t86 = t107 * t133;
t87 = t110 * t133;
t174 = qJD(5) * t56 + t106 * t87 + t109 * t86 - t84 * t137;
t57 = t106 * t95 + t109 * t96;
t173 = -qJD(5) * t57 - t106 * t86 + t109 * t87 + t85 * t137;
t143 = qJD(4) * t107;
t102 = t110 * qJD(1);
t156 = pkin(2) * qJD(3);
t134 = qJD(2) * t156;
t126 = t111 * t134;
t160 = qJD(4) * t102 + t110 * t126;
t108 = sin(qJ(3));
t91 = pkin(7) * t105 + t108 * t157;
t45 = -t91 * t143 + t160;
t120 = t107 * t126;
t145 = qJD(1) * t107;
t69 = t110 * t91 + t145;
t46 = -qJD(4) * t69 - t120;
t172 = -t46 * t107 + t110 * t45;
t73 = t85 * t105;
t39 = -mrSges(6,1) * t72 + mrSges(6,2) * t73;
t101 = -t110 * pkin(4) - pkin(3);
t75 = t101 * t105 - t137;
t171 = -m(6) * t75 - t39;
t68 = -t107 * t91 + t102;
t170 = -t107 * t68 + t110 * t69;
t131 = pkin(8) * t105 + t91;
t51 = t110 * t131 + t145;
t168 = t73 / 0.2e1;
t99 = pkin(2) * t108 + pkin(7);
t165 = -pkin(8) - t99;
t164 = mrSges(6,3) * t72;
t163 = Ifges(6,4) * t73;
t162 = pkin(2) * t111;
t161 = t73 * mrSges(6,3);
t159 = Ifges(5,4) * t107;
t155 = t106 * t51;
t153 = t109 * t51;
t149 = Ifges(5,5) * qJD(4);
t148 = Ifges(5,6) * qJD(4);
t147 = t105 * t107;
t146 = t105 * t110;
t144 = qJD(3) * t108;
t142 = qJD(4) * t110;
t141 = qJD(5) * t106;
t140 = qJD(5) * t109;
t139 = -qJD(2) - t105;
t138 = -qJD(3) + t105;
t136 = t111 * t156;
t135 = pkin(4) * t143;
t132 = t105 * t143;
t130 = qJD(4) * t165;
t127 = t108 * t134;
t123 = t131 * t107;
t50 = t102 - t123;
t47 = qJD(4) * pkin(4) + t50;
t22 = t109 * t47 - t155;
t23 = t106 * t47 + t153;
t40 = -qJD(4) * t123 + t160;
t41 = -qJD(4) * t51 - t120;
t4 = -qJD(5) * t23 - t106 * t40 + t109 * t41;
t124 = -t22 * t48 - t4 * t85;
t122 = -mrSges(5,1) * t110 + mrSges(5,2) * t107;
t82 = t165 * t107;
t83 = t110 * t99 + t103;
t43 = -t106 * t83 + t109 * t82;
t44 = t106 * t82 + t109 * t83;
t89 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t147;
t90 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t146;
t121 = t107 * t89 - t110 * t90;
t119 = (Ifges(5,2) * t110 + t159) * t105;
t118 = (mrSges(5,1) * t107 + mrSges(5,2) * t110) * qJD(4);
t115 = m(5) * t170 - t121;
t3 = qJD(5) * t22 + t106 * t41 + t109 * t40;
t31 = Ifges(6,2) * t72 + Ifges(6,6) * t104 + t163;
t65 = Ifges(6,4) * t72;
t32 = Ifges(6,1) * t73 + Ifges(6,5) * t104 + t65;
t114 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + t22 * t164 + t31 * t168 - t75 * (mrSges(6,1) * t73 + mrSges(6,2) * t72) - t73 * (Ifges(6,1) * t72 - t163) / 0.2e1 - t104 * (Ifges(6,5) * t72 - Ifges(6,6) * t73) / 0.2e1 - Ifges(6,6) * t37 + Ifges(6,5) * t36 + (-Ifges(6,2) * t73 + t32 + t65) * t177;
t113 = m(5) * (-t68 * t142 - t69 * t143 + t172) - t89 * t142 - t90 * t143;
t70 = t119 + t148;
t97 = Ifges(5,4) * t146;
t71 = Ifges(5,1) * t147 + t149 + t97;
t78 = pkin(4) * t132 + t127;
t92 = -t105 * pkin(3) - t137;
t112 = t122 * t127 + qJD(4) ^ 2 * (Ifges(5,5) * t110 - Ifges(5,6) * t107) / 0.2e1 + t104 * (Ifges(6,5) * t48 - Ifges(6,6) * t49) / 0.2e1 + t78 * (-mrSges(6,1) * t84 + mrSges(6,2) * t85) + t75 * (mrSges(6,1) * t49 + mrSges(6,2) * t48) + t32 * t178 - t49 * t31 / 0.2e1 + (Ifges(5,1) * t110 - t159) * t132 + t92 * t118 - (t119 + t70) * t143 / 0.2e1 + (t49 * t177 - t84 * t37) * Ifges(6,2) + (t48 * t168 + t36 * t85) * Ifges(6,1) + (-t23 * t49 + t3 * t84) * mrSges(6,3) + ((-t107 * t69 - t110 * t68) * qJD(4) + t172) * mrSges(5,3) + (-t49 * t168 + t72 * t178 - t175 + t176) * Ifges(6,4) + (t71 + (0.3e1 * Ifges(5,4) * t110 + (Ifges(5,1) - 0.2e1 * Ifges(5,2)) * t107) * t105) * t142 / 0.2e1;
t100 = -pkin(3) - t162;
t93 = t101 - t162;
t88 = pkin(2) * t144 + t135;
t79 = t122 * t105;
t74 = t105 * t118;
t59 = -t107 * t136 + t110 * t130;
t58 = t107 * t130 + t110 * t136;
t53 = mrSges(6,1) * t104 - t161;
t52 = -mrSges(6,2) * t104 + t164;
t27 = t109 * t50 - t155;
t26 = -t106 * t50 - t153;
t10 = mrSges(6,1) * t37 + mrSges(6,2) * t36;
t9 = -qJD(5) * t44 - t106 * t58 + t109 * t59;
t8 = qJD(5) * t43 + t106 * t59 + t109 * t58;
t1 = [t48 * t52 - t49 * t53 + (-t175 - t176) * mrSges(6,3) + m(5) * (t107 * t45 + t110 * t46) + m(6) * (-t22 * t49 + t23 * t48 + t3 * t85 + t4 * t84) + ((-t107 ^ 2 - t110 ^ 2) * t105 * mrSges(5,3) + t115) * qJD(4); m(6) * (t22 * t9 + t23 * t8 + t3 * t44 + t4 * t43 + t75 * t88 + t78 * t93) + (-t36 * t43 - t37 * t44 + t124) * mrSges(6,3) + t113 * t99 + ((m(5) * (qJD(2) * t100 + t92) + t79 + t139 * mrSges(4,1)) * t108 + (t139 * mrSges(4,2) + t115) * t111) * t156 + t88 * t39 + t93 * t10 + t100 * t74 + t8 * t52 + t9 * t53 + t112; (-t56 * t36 - t57 * t37 + t124) * mrSges(6,3) + t113 * pkin(7) + ((t138 * mrSges(4,2) + t121) * t111 + (t138 * mrSges(4,1) + t171 - t79) * t108 + (-pkin(3) * t144 - t108 * t92 - t111 * t170) * m(5)) * t157 + t101 * t10 - pkin(3) * t74 + t39 * t135 + t174 * t52 + t173 * t53 + t112 + (t101 * t78 + t135 * t75 + t173 * t22 + t174 * t23 + t3 * t57 + t4 * t56) * m(6); -m(6) * (t22 * t26 + t23 * t27) + (m(6) * (t106 * t3 + t109 * t4 + t23 * t140 - t22 * t141) + t52 * t140 - t53 * t141 + (-t106 * t37 - t109 * t36) * mrSges(6,3)) * pkin(4) + ((-t92 * mrSges(5,2) + t149 / 0.2e1 - t97 / 0.2e1 - t71 / 0.2e1 + t68 * mrSges(5,3)) * t110 + (-t92 * mrSges(5,1) - t148 / 0.2e1 + t70 / 0.2e1 + t69 * mrSges(5,3) + (t159 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t110) * t105 + t171 * pkin(4)) * t107) * t105 + t23 * t161 + t69 * t89 - t68 * t90 - t27 * t52 - t26 * t53 - t45 * mrSges(5,2) + t46 * mrSges(5,1) + t114; -t22 * t52 + (t53 + t161) * t23 + t114;];
tauc = t1(:);
