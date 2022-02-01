% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:01:41
% EndTime: 2022-01-20 12:01:47
% DurationCPUTime: 1.91s
% Computational Cost: add. (4159->275), mult. (6725->405), div. (0->0), fcn. (3843->8), ass. (0->157)
t135 = sin(qJ(5));
t136 = sin(qJ(4));
t139 = cos(qJ(5));
t140 = cos(qJ(4));
t110 = -t135 * t136 + t139 * t140;
t131 = qJD(4) + qJD(5);
t71 = t131 * t110;
t219 = t71 / 0.2e1;
t132 = qJD(1) + qJD(2);
t129 = qJD(3) + t132;
t92 = t110 * t129;
t218 = -t92 / 0.2e1;
t217 = -mrSges(5,1) * t140 + mrSges(5,2) * t136 - mrSges(4,1);
t111 = t135 * t140 + t136 * t139;
t209 = -pkin(9) - pkin(8);
t169 = qJD(4) * t209;
t112 = t136 * t169;
t113 = t140 * t169;
t118 = t209 * t136;
t130 = t140 * pkin(9);
t119 = pkin(8) * t140 + t130;
t84 = t118 * t135 + t119 * t139;
t142 = cos(qJ(2));
t196 = pkin(1) * qJD(1);
t171 = t142 * t196;
t115 = t132 * pkin(2) + t171;
t137 = sin(qJ(3));
t138 = sin(qJ(2));
t172 = t138 * t196;
t121 = t137 * t172;
t141 = cos(qJ(3));
t88 = t115 * t141 - t121;
t216 = -t84 * qJD(5) + t111 * t88 - t112 * t135 + t113 * t139;
t83 = t118 * t139 - t119 * t135;
t215 = t83 * qJD(5) - t110 * t88 + t112 * t139 + t113 * t135;
t104 = t141 * t171 - t121;
t124 = pkin(2) * t137 + pkin(8);
t201 = -pkin(9) - t124;
t106 = t201 * t136;
t187 = t124 * t140;
t107 = t130 + t187;
t70 = t106 * t135 + t107 * t139;
t163 = qJD(4) * t201;
t180 = qJD(3) * t141;
t170 = pkin(2) * t180;
t85 = t136 * t163 + t140 * t170;
t86 = -t136 * t170 + t140 * t163;
t214 = -t70 * qJD(5) + t111 * t104 - t135 * t85 + t139 * t86;
t69 = t106 * t139 - t107 * t135;
t213 = t69 * qJD(5) - t110 * t104 + t135 * t86 + t139 * t85;
t212 = mrSges(3,1) * t138 + mrSges(3,2) * t142;
t89 = t115 * t137 + t141 * t172;
t74 = pkin(8) * t129 + t89;
t162 = (t136 ^ 2 + t140 ^ 2) * t74;
t93 = t111 * t129;
t210 = t93 / 0.2e1;
t208 = mrSges(6,3) * t92;
t207 = Ifges(6,4) * t93;
t206 = pkin(2) * t141;
t205 = pkin(4) * t140;
t177 = qJD(4) * t140;
t181 = qJD(3) * t137;
t184 = t137 * t138;
t145 = (-t138 * t181 + (t141 * t142 - t184) * qJD(2)) * pkin(1);
t60 = qJD(1) * t145 + t115 * t180;
t193 = t136 * t60;
t28 = -t74 * t177 - t193;
t204 = t28 * mrSges(5,3);
t203 = t93 * mrSges(6,3);
t126 = pkin(1) * t142 + pkin(2);
t183 = t138 * t141;
t182 = pkin(1) * t183 + t137 * t126;
t102 = pkin(8) + t182;
t202 = -pkin(9) - t102;
t198 = Ifges(5,4) * t136;
t167 = pkin(9) * t129 + t74;
t63 = t167 * t140;
t195 = t135 * t63;
t194 = t136 * t28;
t192 = t139 * t63;
t178 = qJD(4) * t136;
t55 = t140 * t60;
t27 = -t74 * t178 + t55;
t191 = t140 * t27;
t189 = Ifges(5,5) * qJD(4);
t188 = Ifges(5,6) * qJD(4);
t186 = t129 * t136;
t185 = t129 * t140;
t108 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t186;
t179 = qJD(4) * t108;
t176 = qJD(5) * t135;
t175 = qJD(5) * t139;
t128 = pkin(4) * t178;
t127 = -pkin(3) - t205;
t168 = t129 * t178;
t164 = qJD(4) * t202;
t161 = t217 * t129;
t160 = -pkin(1) * t184 + t126 * t141;
t54 = -mrSges(6,1) * t92 + mrSges(6,2) * t93;
t158 = -t161 - t54;
t101 = -pkin(3) - t160;
t62 = t167 * t136;
t59 = qJD(4) * pkin(4) - t62;
t13 = t139 * t59 - t195;
t14 = t135 * t59 + t192;
t156 = qJD(4) * t167;
t15 = -t136 * t156 + t55;
t16 = -t140 * t156 - t193;
t4 = -t14 * qJD(5) - t135 * t15 + t139 * t16;
t157 = -t111 * t4 - t13 * t71;
t81 = t202 * t136;
t82 = t102 * t140 + t130;
t42 = -t135 * t82 + t139 * t81;
t43 = t135 * t81 + t139 * t82;
t154 = t191 - t194;
t109 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t185;
t153 = t108 * t140 + t109 * t136;
t152 = t137 * t142 + t183;
t151 = (Ifges(5,2) * t140 + t198) * t129;
t150 = (mrSges(5,1) * t136 + mrSges(5,2) * t140) * qJD(4);
t149 = t129 * mrSges(4,2) + t136 * t108 - t140 * t109;
t3 = t13 * qJD(5) + t135 * t16 + t139 * t15;
t47 = Ifges(6,2) * t92 + Ifges(6,6) * t131 + t207;
t87 = Ifges(6,4) * t92;
t48 = Ifges(6,1) * t93 + Ifges(6,5) * t131 + t87;
t51 = t71 * t129;
t72 = t131 * t111;
t52 = t72 * t129;
t66 = t127 * t129 - t88;
t146 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + t13 * t208 + t47 * t210 - t66 * (mrSges(6,1) * t93 + mrSges(6,2) * t92) - t93 * (Ifges(6,1) * t92 - t207) / 0.2e1 - t131 * (Ifges(6,5) * t92 - Ifges(6,6) * t93) / 0.2e1 - Ifges(6,6) * t52 + Ifges(6,5) * t51 + (-Ifges(6,2) * t93 + t48 + t87) * t218;
t144 = (t152 * qJD(2) + t138 * t180) * pkin(1);
t68 = t126 * t181 + t144;
t61 = qJD(1) * t144 + t115 * t181;
t44 = pkin(4) * t168 + t61;
t73 = -pkin(3) * t129 - t88;
t90 = t151 + t188;
t117 = Ifges(5,4) * t185;
t91 = Ifges(5,1) * t186 + t117 + t189;
t143 = t131 * (Ifges(6,5) * t71 - Ifges(6,6) * t72) / 0.2e1 + t44 * (-mrSges(6,1) * t110 + mrSges(6,2) * t111) + qJD(4) ^ 2 * (Ifges(5,5) * t140 - Ifges(5,6) * t136) / 0.2e1 + (Ifges(5,1) * t140 - t198) * t168 + t48 * t219 - t72 * t47 / 0.2e1 + t66 * (mrSges(6,1) * t72 + mrSges(6,2) * t71) + mrSges(5,3) * t191 + t73 * t150 + t217 * t61 - (t90 + t151) * t178 / 0.2e1 + (-t110 * t52 + t72 * t218) * Ifges(6,2) + (t51 * t111 + t71 * t210) * Ifges(6,1) + (t3 * t110 - t14 * t72) * mrSges(6,3) + (t110 * t51 - t52 * t111 - t72 * t210 + t92 * t219) * Ifges(6,4) + (t91 + (0.3e1 * Ifges(5,4) * t140 + (Ifges(5,1) - 0.2e1 * Ifges(5,2)) * t136) * t129) * t177 / 0.2e1;
t125 = -pkin(3) - t206;
t116 = t127 - t206;
t114 = pkin(2) * t181 + t128;
t103 = t152 * t196;
t97 = t101 - t205;
t94 = t129 * t150;
t76 = mrSges(6,1) * t131 - t203;
t75 = -mrSges(6,2) * t131 + t208;
t67 = t126 * t180 + t145;
t65 = t128 + t68;
t37 = -t136 * t67 + t140 * t164;
t36 = t136 * t164 + t140 * t67;
t20 = -t139 * t62 - t195;
t19 = t135 * t62 - t192;
t12 = mrSges(6,1) * t52 + mrSges(6,2) * t51;
t6 = -t43 * qJD(5) - t135 * t36 + t139 * t37;
t5 = t42 * qJD(5) + t135 * t37 + t139 * t36;
t1 = [(-t42 * t51 - t43 * t52 + t157) * mrSges(6,3) + (-t129 * t67 - t60) * mrSges(4,2) + t101 * t94 + t97 * t12 + t143 + t6 * t76 + t5 * t75 + t65 * t54 + t161 * t68 + (-t102 * t179 + t67 * t109) * t140 + (-qJD(4) * t102 * t109 - t67 * t108 - t204) * t136 + m(5) * (t101 * t61 + t154 * t102 + t67 * t162 + t68 * t73) + m(4) * (-t61 * t160 + t60 * t182 + t89 * t67 - t88 * t68) + m(6) * (t13 * t6 + t14 * t5 + t3 * t43 + t4 * t42 + t44 * t97 + t65 * t66) + t212 * qJD(2) * pkin(1) * (-qJD(1) - t132); -mrSges(5,3) * t194 + t114 * t54 + t116 * t12 + t125 * t94 + t143 + (-t51 * t69 - t52 * t70 + t157) * mrSges(6,3) + t149 * t104 + t158 * t103 - t60 * mrSges(4,2) - t153 * t124 * qJD(4) - m(4) * (-t103 * t88 + t104 * t89) + m(5) * (-t124 * t194 + t125 * t61 + t27 * t187) + t214 * t76 + (m(4) * (t137 * t60 - t141 * t61) + (t161 * t137 - t149 * t141 + m(5) * (t137 * t73 + t141 * t162) + m(4) * (-t137 * t88 + t141 * t89)) * qJD(3)) * pkin(2) + t213 * t75 - m(5) * (t103 * t73 + t104 * t162) + t212 * t196 * (-qJD(2) + t132) + (t116 * t44 + t3 * t70 + t4 * t69 + (t114 - t103) * t66 + t213 * t14 + t214 * t13) * m(6); (-t51 * t83 - t52 * t84 + t157) * mrSges(6,3) + (t129 * t88 - t60) * mrSges(4,2) + (-t204 + t88 * t108 + (pkin(4) * t54 - pkin(8) * t109) * qJD(4)) * t136 + t127 * t12 - pkin(3) * t94 + t143 + t158 * t89 + (-pkin(8) * t179 - t88 * t109) * t140 + t215 * t75 + t216 * t76 + (t127 * t44 + t3 * t84 + t4 * t83 + (t128 - t89) * t66 + t215 * t14 + t216 * t13) * m(6) + (-pkin(3) * t61 + t154 * pkin(8) - t88 * t162 - t73 * t89) * m(5); t153 * t74 + t146 - t20 * t75 - t19 * t76 - t27 * mrSges(5,2) + t28 * mrSges(5,1) + t14 * t203 - m(6) * (t13 * t19 + t14 * t20) + (m(6) * (-t13 * t176 + t135 * t3 + t139 * t4 + t14 * t175) + t75 * t175 - t76 * t176 + (-t135 * t52 - t139 * t51) * mrSges(6,3)) * pkin(4) + ((-t117 / 0.2e1 - t91 / 0.2e1 - t73 * mrSges(5,2) + t189 / 0.2e1) * t140 + (t90 / 0.2e1 - t73 * mrSges(5,1) - t188 / 0.2e1 + (t198 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t140) * t129 + (-m(6) * t66 - t54) * pkin(4)) * t136) * t129; t146 - t13 * t75 + (t76 + t203) * t14;];
tauc = t1(:);
