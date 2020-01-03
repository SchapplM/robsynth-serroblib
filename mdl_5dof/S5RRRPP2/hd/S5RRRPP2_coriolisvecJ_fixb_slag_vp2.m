% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPP2
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:29
% EndTime: 2019-12-31 20:51:34
% DurationCPUTime: 1.88s
% Computational Cost: add. (1314->276), mult. (2388->349), div. (0->0), fcn. (952->4), ass. (0->138)
t206 = mrSges(4,1) + mrSges(5,1);
t205 = Ifges(6,2) + Ifges(5,3);
t113 = sin(qJ(3));
t115 = cos(qJ(3));
t108 = qJD(1) + qJD(2);
t114 = sin(qJ(2));
t178 = pkin(1) * qJD(1);
t89 = pkin(7) * t108 + t114 * t178;
t204 = (t113 ^ 2 + t115 ^ 2) * t89;
t109 = qJD(3) * qJD(4);
t162 = qJD(3) * t113;
t116 = cos(qJ(2));
t177 = pkin(1) * qJD(2);
t153 = qJD(1) * t177;
t144 = t116 * t153;
t92 = t115 * t144;
t21 = -t89 * t162 + t92;
t17 = t109 + t21;
t161 = qJD(3) * t115;
t22 = t113 * t144 + t89 * t161;
t176 = t113 * t22;
t203 = t115 * t17 + t176;
t202 = (Ifges(6,4) + Ifges(5,5)) * t113;
t168 = t108 * t115;
t88 = mrSges(5,2) * t168 + qJD(3) * mrSges(5,3);
t185 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t168 + t88;
t169 = t108 * t113;
t187 = (-mrSges(5,2) - mrSges(4,3)) * t169 + t206 * qJD(3);
t201 = t185 * t113 + t187 * t115;
t105 = t113 * qJ(4);
t107 = t115 * pkin(3);
t198 = t105 + t107;
t77 = t113 * t89;
t38 = qJ(5) * t169 - t77;
t196 = -t38 + qJD(4);
t110 = qJD(3) * qJ(4);
t78 = t115 * t89;
t65 = t78 + t110;
t175 = t113 * t65;
t48 = -qJD(3) * pkin(3) + qJD(4) + t77;
t134 = t115 * t48 - t175;
t173 = t115 * t21;
t195 = (m(5) * t134 - t201) * qJD(3) + m(5) * t203 + m(4) * (t173 + t176);
t194 = qJD(3) ^ 2;
t193 = pkin(3) + pkin(4);
t39 = -qJ(5) * t168 + t78;
t20 = t110 + t39;
t192 = t20 * mrSges(6,3);
t191 = t48 * mrSges(5,2);
t190 = mrSges(5,2) - mrSges(6,3);
t189 = pkin(7) - qJ(5);
t137 = -mrSges(5,1) * t115 - mrSges(5,3) * t113;
t69 = t137 * t108;
t136 = mrSges(6,1) * t115 + mrSges(6,2) * t113;
t70 = t136 * t108;
t188 = t69 - t70;
t86 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t168;
t186 = t86 + t88;
t184 = Ifges(4,4) * t113;
t181 = Ifges(6,4) * t115;
t179 = Ifges(5,5) * t115;
t12 = -t193 * qJD(3) + t196;
t171 = t12 * t115;
t170 = qJ(4) * t115;
t167 = t113 * t116;
t166 = t115 * t116;
t165 = t65 * qJD(3);
t101 = pkin(1) * t114 + pkin(7);
t164 = -qJ(5) + t101;
t163 = qJD(2) * t114;
t160 = qJD(5) * t115;
t104 = t113 * qJD(4);
t159 = -qJD(1) - t108;
t158 = -qJD(2) + t108;
t156 = pkin(2) + t198;
t155 = t116 * t178;
t154 = pkin(1) * t163;
t152 = mrSges(6,1) * t162;
t67 = pkin(3) * t162 - qJ(4) * t161 - t104;
t102 = -pkin(1) * t116 - pkin(2);
t151 = t108 * t161;
t94 = t189 * t115;
t149 = t162 / 0.2e1;
t147 = pkin(2) + t105;
t146 = qJ(5) * t108 - t89;
t82 = t164 * t115;
t145 = t114 * t153;
t143 = t108 * t149;
t140 = t116 * t177 - qJD(5);
t3 = -t108 * t160 + t146 * t162 + t109 + t92;
t5 = (-qJD(5) * t108 + t144) * t113 - t146 * t161;
t139 = -t5 * t113 - t3 * t115;
t138 = -mrSges(4,1) * t115 + mrSges(4,2) * t113;
t135 = pkin(3) * t113 - t170;
t79 = t102 - t198;
t131 = -t193 * t113 + t170;
t130 = (Ifges(5,1) * t113 - t179) * t108;
t129 = (Ifges(6,1) * t113 - t181) * t108;
t128 = (Ifges(4,2) * t115 + t184) * t108;
t40 = -pkin(4) * t162 - t67;
t127 = (mrSges(4,1) * t113 + mrSges(4,2) * t115) * qJD(3);
t126 = (mrSges(5,1) * t113 - mrSges(5,3) * t115) * qJD(3);
t19 = -t155 + (-t147 - t107) * t108;
t96 = Ifges(5,5) * t169;
t49 = Ifges(5,6) * qJD(3) - Ifges(5,3) * t168 + t96;
t97 = Ifges(6,4) * t169;
t50 = -Ifges(6,2) * t168 - Ifges(6,6) * qJD(3) + t97;
t51 = Ifges(4,6) * qJD(3) + t128;
t52 = -Ifges(6,5) * qJD(3) + t129;
t53 = Ifges(5,4) * qJD(3) + t130;
t98 = Ifges(4,4) * t168;
t54 = Ifges(4,1) * t169 + Ifges(4,5) * qJD(3) + t98;
t6 = -t145 + (t131 * qJD(3) + t104) * t108;
t8 = t155 + qJD(5) + (t193 * t115 + t147) * t108;
t9 = t145 + (t135 * qJD(3) - t104) * t108;
t90 = -pkin(2) * t108 - t155;
t118 = -t194 * (Ifges(6,5) * t115 + Ifges(6,6) * t113) / 0.2e1 + t6 * t136 + t9 * t137 + t8 * (-mrSges(6,1) * t113 + mrSges(6,2) * t115) * qJD(3) + t19 * t126 + t90 * t127 + mrSges(4,3) * t173 + t161 * t191 + t162 * t192 + t138 * t145 - (t205 * t113 + t179 + t181) * t151 + ((Ifges(5,4) + Ifges(4,5)) * t115 + (-Ifges(4,6) + Ifges(5,6)) * t113) * t194 / 0.2e1 - (t128 + t51) * t162 / 0.2e1 + (t50 + t49) * t149 + (-t205 * t115 + t202) * t143 + t203 * mrSges(5,2) + 0.2e1 * (-t184 + (Ifges(4,1) + Ifges(5,1) + Ifges(6,1)) * t115 + t202) * t143 + ((0.3e1 * Ifges(4,4) * t115 + (Ifges(4,1) - 0.2e1 * Ifges(4,2)) * t113) * t108 + t130 + t129 + t54 + t53 + t52) * t161 / 0.2e1;
t106 = t115 * pkin(4);
t99 = qJ(5) * t162;
t95 = mrSges(6,2) * t151;
t93 = t189 * t113;
t83 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t169;
t81 = t164 * t113;
t80 = t106 + t156;
t72 = t135 * t108;
t71 = t138 * t108;
t68 = qJD(3) * t94 - qJD(5) * t113;
t66 = -pkin(7) * t162 - t160 + t99;
t64 = t106 - t79;
t57 = t108 * t127;
t56 = -t108 * t152 + t95;
t55 = t108 * t126;
t41 = t67 + t154;
t31 = t131 * t108;
t18 = t40 - t154;
t16 = qJD(3) * t82 + t140 * t113;
t15 = -t101 * t162 + t140 * t115 + t99;
t1 = [t139 * mrSges(6,3) + m(6) * (t12 * t16 + t15 * t20 + t8 * t18 + t3 * t82 + t5 * t81 + t6 * t64) + (-mrSges(5,2) * t175 + (-t171 + (t113 * t82 - t115 * t81) * t108) * mrSges(6,3)) * qJD(3) + t118 + ((t71 + m(4) * (qJD(1) * t102 + t90) + t159 * mrSges(3,1)) * t114 + (t185 * t115 - t187 * t113 + t159 * mrSges(3,2) + m(5) * (t113 * t48 + t115 * t65) + m(4) * t204) * t116) * t177 + t102 * t57 + t16 * t83 + t15 * t86 + mrSges(4,3) * t176 + t64 * t56 + t41 * t69 + t18 * t70 + t79 * t55 + m(5) * (t19 * t41 + t79 * t9) + t195 * t101; t195 * pkin(7) + m(6) * (t12 * t68 + t20 * t66 + t3 * t94 + t40 * t8 + t5 * t93 + t6 * t80) + t118 + ((t158 * mrSges(3,1) - t188 - t71) * t114 + (t158 * mrSges(3,2) + (-t86 - t185) * t115 + (-t83 + t187) * t113) * t116 - m(5) * (t114 * t19 + t65 * t166 + t48 * t167) - m(6) * (-t114 * t8 + t12 * t167 + t20 * t166) + (-pkin(2) * t163 - t114 * t90 - t116 * t204) * m(4)) * t178 + ((-t171 + (t113 * t94 - t115 * t93) * t108) * qJD(3) + t139) * mrSges(6,3) + (-mrSges(5,2) * t165 + t22 * mrSges(4,3)) * t113 + m(5) * (-t156 * t9 + t19 * t67) + t68 * t83 + t66 * t86 - t156 * t55 + t67 * t69 + t40 * t70 + t80 * t56 - pkin(2) * t57; -t5 * mrSges(6,1) - t21 * mrSges(4,2) + t3 * mrSges(6,2) + t17 * mrSges(5,3) - t31 * t70 - t38 * t86 - t39 * t83 - t72 * t69 - t206 * t22 + t186 * qJD(4) + t201 * t89 + ((-t49 / 0.2e1 - t50 / 0.2e1 + t51 / 0.2e1 - t90 * mrSges(4,1) + t8 * mrSges(6,1) - t19 * mrSges(5,1) - t96 / 0.2e1 - t97 / 0.2e1 - t192 + t65 * mrSges(5,2) + Ifges(4,4) * t169 / 0.2e1 + (Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1 - Ifges(6,6) / 0.2e1 - t190 * qJ(4)) * qJD(3)) * t113 + (-t98 / 0.2e1 - t52 / 0.2e1 - t53 / 0.2e1 - t54 / 0.2e1 - t90 * mrSges(4,2) - t8 * mrSges(6,2) + t19 * mrSges(5,3) + t12 * mrSges(6,3) - t191 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t168 + (Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1 - Ifges(6,5) / 0.2e1 + t193 * mrSges(6,3) - pkin(3) * mrSges(5,2)) * qJD(3) + (-Ifges(5,1) / 0.2e1 - Ifges(6,1) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t169) * t115) * t108 + (t3 * qJ(4) - t12 * t39 - t193 * t5 + t196 * t20 - t8 * t31) * m(6) + (-t22 * pkin(3) + t17 * qJ(4) + t65 * qJD(4) - t134 * t89 - t19 * t72) * m(5); t188 * t169 + (t190 * t168 - t186) * qJD(3) + (-t20 * qJD(3) - t8 * t169 + t5) * m(6) + (t19 * t169 - t165 + t22) * m(5); t95 + m(6) * t6 + (-t152 + t113 * t83 + t115 * t86 - m(6) * (-t113 * t12 - t115 * t20)) * t108;];
tauc = t1(:);
