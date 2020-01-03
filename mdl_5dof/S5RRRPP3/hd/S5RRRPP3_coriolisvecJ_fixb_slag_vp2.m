% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:08
% EndTime: 2019-12-31 20:53:12
% DurationCPUTime: 1.73s
% Computational Cost: add. (1309->273), mult. (2383->341), div. (0->0), fcn. (942->4), ass. (0->136)
t203 = mrSges(4,1) - mrSges(5,2);
t111 = (qJD(1) + qJD(2));
t202 = 2 * t111;
t115 = sin(qJ(3));
t117 = cos(qJ(3));
t116 = sin(qJ(2));
t176 = pkin(1) * qJD(1);
t91 = (pkin(7) * t111) + t116 * t176;
t201 = (t115 ^ 2 + t117 ^ 2) * t91;
t118 = cos(qJ(2));
t175 = pkin(1) * qJD(2);
t152 = qJD(1) * t175;
t144 = t118 * t152;
t136 = t115 * t144;
t160 = qJD(3) * t117;
t22 = t160 * t91 + t136;
t173 = t115 * t22;
t161 = qJD(3) * t115;
t94 = t117 * t144;
t21 = -t161 * t91 + t94;
t200 = t117 * t21 + t173;
t78 = t115 * t91;
t17 = -t94 + (-qJD(4) + t78) * qJD(3);
t199 = -t117 * t17 + t173;
t114 = -pkin(3) - qJ(5);
t167 = t111 * t115;
t41 = -pkin(4) * t167 - t78;
t193 = -t41 + qJD(4);
t14 = qJD(3) * t114 + t193;
t50 = -qJD(3) * pkin(3) + qJD(4) + t78;
t198 = t50 * mrSges(5,1) + t14 * mrSges(6,1);
t166 = t111 * t117;
t86 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t166;
t88 = -mrSges(5,1) * t166 - qJD(3) * mrSges(5,3);
t184 = t86 - t88;
t185 = (-mrSges(5,1) - mrSges(4,3)) * t167 + t203 * qJD(3);
t197 = t184 * t115 + t185 * t117;
t79 = t117 * t91;
t42 = pkin(4) * t166 + t79;
t192 = -t42 - qJD(5);
t191 = pkin(1) * t118;
t162 = qJD(3) * qJ(4);
t67 = -t79 - t162;
t188 = t67 * mrSges(5,1);
t187 = mrSges(5,1) + mrSges(6,1);
t139 = mrSges(5,2) * t117 - mrSges(5,3) * t115;
t69 = t139 * t111;
t138 = -mrSges(6,2) * t115 - mrSges(6,3) * t117;
t72 = t138 * t111;
t186 = t69 + t72;
t89 = mrSges(6,1) * t166 + qJD(3) * mrSges(6,2);
t183 = t88 - t89;
t182 = Ifges(4,4) * t115;
t180 = Ifges(5,6) * t115;
t179 = Ifges(5,6) * t117;
t178 = Ifges(6,6) * t115;
t177 = Ifges(6,6) * t117;
t145 = -qJD(5) - t162;
t18 = -t145 + t42;
t174 = t115 * t18;
t102 = t116 * t152;
t106 = pkin(3) * t161;
t169 = t111 * t106 + t102;
t168 = mrSges(6,1) * qJD(3);
t165 = t115 * t118;
t164 = t117 * t118;
t163 = qJD(2) * t116;
t159 = qJD(4) * t115;
t158 = -qJD(1) - t111;
t157 = -qJD(2) + t111;
t155 = t118 * t176;
t154 = t118 * t175;
t153 = mrSges(6,1) * t167;
t151 = pkin(4) * t111 + t91;
t149 = t161 / 0.2e1;
t146 = -t115 * qJ(4) - pkin(2);
t143 = t111 * t149;
t140 = -mrSges(4,1) * t117 + mrSges(4,2) * t115;
t137 = t115 * t67 + t117 * t50;
t135 = -qJ(4) * t117 + qJ(5) * t115;
t93 = -t117 * pkin(3) + t146;
t134 = (Ifges(4,2) * t117 + t182) * t111;
t133 = (-Ifges(5,2) * t115 - t179) * t111;
t132 = (-Ifges(5,3) * t117 - t180) * t111;
t131 = (Ifges(6,3) * t115 - t177) * t111;
t130 = -qJ(4) * t160 - t159;
t129 = (mrSges(4,1) * t115 + mrSges(4,2) * t117) * qJD(3);
t128 = (-mrSges(5,2) * t115 - mrSges(5,3) * t117) * qJD(3);
t127 = (-mrSges(6,2) * t117 + mrSges(6,3) * t115) * qJD(3);
t77 = t114 * t117 + t146;
t68 = t106 + t130;
t19 = qJ(5) * t161 + t117 * t145 + t106 - t159;
t120 = m(5) * t199 + m(4) * t200 + (m(5) * t137 - t197) * qJD(3);
t11 = t111 * t77 - t155;
t12 = t111 * t130 + t169;
t2 = (qJD(3) * t135 - qJD(5) * t117 - t159) * t111 + t169;
t20 = t111 * t93 - t155;
t51 = Ifges(4,6) * qJD(3) + t134;
t99 = Ifges(4,4) * t166;
t52 = Ifges(4,1) * t167 + Ifges(4,5) * qJD(3) + t99;
t53 = Ifges(6,5) * qJD(3) + t131;
t54 = Ifges(5,5) * qJD(3) + t132;
t98 = Ifges(6,6) * t167;
t55 = Ifges(6,4) * qJD(3) - Ifges(6,2) * t166 + t98;
t56 = Ifges(5,4) * qJD(3) + t133;
t8 = t94 + (-t115 * t151 + qJD(4)) * qJD(3);
t9 = t136 + (t117 * t151 - qJD(5)) * qJD(3);
t92 = -t111 * pkin(2) - t155;
t119 = (-Ifges(6,2) * t117 + t178) * t143 + t2 * t138 + t12 * t139 + t11 * t127 + t20 * t128 + t92 * t129 + t161 * t188 + t140 * t102 + 0.2e1 * (t178 - t182 + (Ifges(4,1) + Ifges(6,3)) * t117) * t143 + t198 * t160 + t200 * mrSges(4,3) + (t9 * t115 + t8 * t117) * mrSges(6,1) + t199 * mrSges(5,1) + ((-Ifges(5,4) + Ifges(4,5) + Ifges(6,5)) * t117 + (Ifges(6,4) + Ifges(5,5) - Ifges(4,6)) * t115) * qJD(3) ^ 2 / 0.2e1 - ((-Ifges(5,2) * t117 + t180) * t202 + t134 + t51) * t161 / 0.2e1 + (t132 + t55 + t54) * t149 - ((t177 - t179 + (Ifges(6,2) + Ifges(5,3)) * t115) * t202 + t133 + t56) * t160 / 0.2e1 + ((0.3e1 * Ifges(4,4) * t117 + (Ifges(4,1) - 0.2e1 * Ifges(4,2)) * t115) * t111 + t131 + t53 + t52) * t160 / 0.2e1;
t110 = t117 * pkin(4);
t109 = t115 * pkin(4);
t108 = pkin(1) * t163;
t107 = pkin(4) * t160;
t105 = -pkin(2) - t191;
t104 = pkin(1) * t116 + pkin(7);
t100 = pkin(3) * t167;
t96 = pkin(7) * t117 + t110;
t95 = pkin(7) * t115 + t109;
t87 = -qJD(3) * mrSges(6,3) + t153;
t84 = pkin(7) * t160 + t107;
t83 = (-pkin(4) - pkin(7)) * t161;
t82 = t104 * t117 + t110;
t81 = t104 * t115 + t109;
t80 = t93 - t191;
t71 = -qJ(4) * t166 + t100;
t70 = t140 * t111;
t60 = t77 - t191;
t59 = t111 * t129;
t58 = t111 * t128;
t57 = t111 * t127;
t43 = t108 + t68;
t32 = t111 * t135 + t100;
t31 = t104 * t160 + t115 * t154 + t107;
t30 = t117 * t154 + (-pkin(4) - t104) * t161;
t13 = t108 + t19;
t1 = [m(6) * (t11 * t13 + t14 * t31 + t18 * t30 + t2 * t60 + t8 * t82 + t81 * t9) + t120 * t104 + ((m(4) * (qJD(1) * t105 + t92) + t70 + (t158 * mrSges(3,1))) * t116 + (t184 * t117 - t185 * t115 + (t158 * mrSges(3,2)) + m(5) * (t115 * t50 - t117 * t67) + m(4) * t201) * t118) * t175 + t119 + t30 * t89 + t105 * t59 + t13 * t72 + t80 * t58 + t31 * t87 + t60 * t57 + t43 * t69 + (-t174 + (-t115 * t82 + t117 * t81) * t111) * t168 + m(5) * (t12 * t80 + t20 * t43); t120 * pkin(7) + ((mrSges(3,1) * t157 - t186 - t70) * t116 + (t157 * mrSges(3,2) + (-t86 + t183) * t117 + (-t87 + t185) * t115) * t118 - m(5) * (t116 * t20 - t164 * t67 + t165 * t50) - m(6) * (t11 * t116 + t14 * t165 + t164 * t18) + (-pkin(2) * t163 - t116 * t92 - t118 * t201) * m(4)) * t176 + m(5) * (t12 * t93 + t20 * t68) + m(6) * (t11 * t19 + t14 * t84 + t18 * t83 + t2 * t77 + t8 * t96 + t9 * t95) + t119 + (-t174 + (-t115 * t96 + t117 * t95) * t111) * t168 + t83 * t89 + t93 * t58 + t19 * t72 + t77 * t57 + t84 * t87 - pkin(2) * t59 + t68 * t69; -t21 * mrSges(4,2) + t8 * mrSges(6,2) - t17 * mrSges(5,3) - t9 * mrSges(6,3) - t32 * t72 - t41 * t89 - t71 * t69 + t192 * t87 - t203 * t22 - t183 * qJD(4) + t197 * t91 + ((-t92 * mrSges(4,2) + t20 * mrSges(5,3) + t11 * mrSges(6,2) - t99 / 0.2e1 - t52 / 0.2e1 - t53 / 0.2e1 + t56 / 0.2e1 + (-Ifges(5,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t166 - t198) * t117 + (-t92 * mrSges(4,1) + t20 * mrSges(5,2) - t11 * mrSges(6,3) - t98 / 0.2e1 + t51 / 0.2e1 - t54 / 0.2e1 - t55 / 0.2e1 - t188 + t18 * mrSges(6,1) + (Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1) * t167 + (-Ifges(5,2) / 0.2e1 + Ifges(5,3) / 0.2e1 - Ifges(4,1) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(6,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t166) * t115 + ((Ifges(6,5) / 0.2e1 - Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1 - pkin(3) * mrSges(5,1) + t114 * mrSges(6,1)) * t117 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1 - t187 * qJ(4)) * t115) * qJD(3)) * t111 + (t8 * qJ(4) - t11 * t32 + t9 * t114 + t14 * t192 + t18 * t193) * m(6) + (-t22 * pkin(3) - t17 * qJ(4) - t67 * qJD(4) - t137 * t91 - t20 * t71) * m(5); t186 * t167 + (t166 * t187 + t183) * qJD(3) + (-t18 * qJD(3) + t11 * t167 + t9) * m(6) + (t67 * qJD(3) + t167 * t20 + t22) * m(5); m(6) * t8 + (m(6) * t11 + t72) * t166 + (m(6) * t14 - t153 + t87) * qJD(3);];
tauc = t1(:);
