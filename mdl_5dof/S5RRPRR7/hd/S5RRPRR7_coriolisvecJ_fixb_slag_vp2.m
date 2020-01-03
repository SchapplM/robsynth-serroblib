% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:28
% EndTime: 2019-12-31 20:15:31
% DurationCPUTime: 1.46s
% Computational Cost: add. (2365->229), mult. (3562->320), div. (0->0), fcn. (1894->6), ass. (0->137)
t115 = qJD(1) + qJD(2);
t118 = sin(qJ(5));
t119 = sin(qJ(4));
t121 = cos(qJ(5));
t122 = cos(qJ(4));
t168 = t121 * t122;
t199 = -t118 * t119 + t168;
t74 = t199 * t115;
t114 = qJD(4) + qJD(5);
t142 = mrSges(5,1) * t122 - mrSges(5,2) * t119;
t130 = t142 * qJD(4);
t182 = Ifges(5,4) * t122;
t131 = (-Ifges(5,2) * t119 + t182) * t115;
t183 = Ifges(5,4) * t119;
t132 = (Ifges(5,1) * t122 - t183) * t115;
t133 = t122 * (-Ifges(5,1) * t119 - t182);
t134 = t119 * (-Ifges(5,2) * t122 - t183);
t136 = t118 * t122 + t121 * t119;
t141 = mrSges(5,1) * t119 + mrSges(5,2) * t122;
t120 = sin(qJ(2));
t178 = pkin(1) * qJD(2);
t155 = qJD(1) * t178;
t147 = t120 * t155;
t163 = qJD(4) * t122;
t164 = qJD(4) * t119;
t180 = Ifges(5,6) * t122;
t181 = Ifges(5,5) * t119;
t195 = t74 / 0.2e1;
t124 = -pkin(2) - pkin(7);
t123 = cos(qJ(2));
t179 = pkin(1) * qJD(1);
t143 = -t123 * t179 + qJD(3);
t79 = t124 * t115 + t143;
t98 = t122 * t147;
t44 = -t79 * t164 + t98;
t45 = t119 * t147 + t79 * t163;
t198 = -t119 * t45 - t122 * t44;
t171 = t115 * t119;
t49 = -pkin(8) * t171 + t119 * t79;
t174 = t121 * t49;
t170 = t115 * t122;
t159 = pkin(8) * t170;
t50 = t122 * t79 - t159;
t46 = qJD(4) * pkin(4) + t50;
t23 = t118 * t46 + t174;
t177 = t118 * t49;
t22 = t121 * t46 - t177;
t40 = t98 + (pkin(8) * t115 - t79) * t164;
t41 = -qJD(4) * t159 + t45;
t3 = t22 * qJD(5) + t118 * t40 + t121 * t41;
t162 = qJD(5) * t118;
t47 = t114 * t168 - t118 * t164 - t119 * t162;
t203 = -t3 * t136 - t23 * t47;
t37 = t114 * t74;
t205 = t37 * t136;
t48 = t114 * t136;
t36 = t48 * t115;
t206 = t199 * t36;
t73 = t136 * t115;
t209 = t73 / 0.2e1;
t190 = Ifges(6,4) * t74;
t32 = -Ifges(6,2) * t73 + Ifges(6,6) * t114 + t190;
t67 = Ifges(6,4) * t73;
t33 = Ifges(6,1) * t74 + Ifges(6,5) * t114 - t67;
t112 = pkin(4) * t163;
t102 = qJD(3) + t112;
t146 = t123 * t155;
t64 = t102 * t115 + t146;
t71 = Ifges(5,6) * qJD(4) + t131;
t72 = Ifges(5,5) * qJD(4) + t132;
t113 = t119 * pkin(4);
t108 = qJ(3) + t113;
t156 = t120 * t179;
t75 = t108 * t115 + t156;
t166 = qJD(3) * t115;
t87 = t146 + t166;
t92 = qJ(3) * t115 + t156;
t210 = t198 * mrSges(5,3) + mrSges(4,2) * t147 + qJD(4) ^ 2 * (-t180 - t181) / 0.2e1 + t114 * (-Ifges(6,5) * t48 - Ifges(6,6) * t47) / 0.2e1 + t64 * (mrSges(6,1) * t136 + mrSges(6,2) * t199) + t75 * (mrSges(6,1) * t47 - mrSges(6,2) * t48) - t47 * t32 / 0.2e1 - t48 * t33 / 0.2e1 + t92 * t130 + t87 * t141 + (t133 - t134) * qJD(4) * t115 - (t132 + t72) * t164 / 0.2e1 - (t131 + t71) * t163 / 0.2e1 + (t47 * t209 + t205) * Ifges(6,2) + (-t48 * t195 - t206) * Ifges(6,1) + t203 * mrSges(6,3) + (t136 * t36 - t47 * t195 - t199 * t37 + t48 * t209) * Ifges(6,4);
t185 = -pkin(8) + t124;
t94 = t185 * t119;
t95 = t185 * t122;
t52 = t118 * t95 + t121 * t94;
t111 = pkin(8) * t164;
t84 = -t124 * t164 + t111;
t85 = qJD(4) * t95;
t208 = -t52 * qJD(5) - t118 * t85 + t121 * t84 - t156 * t199;
t51 = -t118 * t94 + t121 * t95;
t207 = t51 * qJD(5) + t118 * t84 + t121 * t85 - t136 * t156;
t39 = mrSges(6,1) * t73 + mrSges(6,2) * t74;
t204 = -m(6) * t75 - t39;
t80 = t141 * t115;
t202 = mrSges(4,3) * t115 + t80;
t201 = -t202 + t204;
t154 = -pkin(1) * t123 - pkin(2);
t106 = -pkin(7) + t154;
t167 = qJD(2) * t120;
t158 = pkin(1) * t167;
t200 = -t106 * t164 + t122 * t158;
t4 = -t23 * qJD(5) - t118 * t41 + t121 * t40;
t144 = t199 * t4 - t22 * t48;
t197 = (t119 ^ 2 + t122 ^ 2) * t79;
t191 = mrSges(6,3) * t73;
t187 = t74 * mrSges(6,3);
t186 = -pkin(8) + t106;
t172 = t123 * t92;
t161 = qJD(5) * t121;
t160 = -qJD(2) + t115;
t157 = t123 * t178;
t83 = t186 * t122;
t107 = pkin(1) * t120 + qJ(3);
t150 = t115 * t158;
t149 = t119 * t158;
t101 = qJD(3) + t157;
t140 = t92 * t101 + t87 * t107;
t82 = t186 * t119;
t43 = t118 * t83 + t121 * t82;
t42 = -t118 * t82 + t121 * t83;
t91 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t171;
t93 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t170;
t138 = t119 * t93 - t122 * t91;
t137 = t87 * qJ(3) + t92 * qJD(3);
t129 = t138 * qJD(4);
t126 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t22 * t191 + t32 * t195 - t75 * (mrSges(6,1) * t74 - mrSges(6,2) * t73) - t74 * (-Ifges(6,1) * t73 - t190) / 0.2e1 - t114 * (-Ifges(6,5) * t73 - Ifges(6,6) * t74) / 0.2e1 - Ifges(6,6) * t37 - Ifges(6,5) * t36 + (-Ifges(6,2) * t74 + t33 - t67) * t209;
t96 = t107 + t113;
t88 = -pkin(2) * t115 + t143;
t86 = t101 + t112;
t76 = t115 * t130;
t56 = qJD(4) * t83 + t149;
t55 = t111 + t200;
t54 = mrSges(6,1) * t114 - t187;
t53 = -mrSges(6,2) * t114 - t191;
t25 = t121 * t50 - t177;
t24 = -t118 * t50 - t174;
t10 = mrSges(6,1) * t37 - mrSges(6,2) * t36;
t9 = -t43 * qJD(5) - t118 * t56 + t121 * t55;
t8 = t42 * qJD(5) + t118 * t55 + t121 * t56;
t1 = [m(6) * (t22 * t9 + t23 * t8 + t3 * t43 + t4 * t42 + t64 * t96 + t75 * t86) + m(4) * ((t154 * qJD(1) + t88) * t158 + t140) + m(5) * (-t106 * t198 + t158 * t197 + t140) + mrSges(4,2) * t150 + t96 * t10 + t107 * t76 + t86 * t39 + t87 * mrSges(4,3) + t8 * t53 + t9 * t54 + t200 * t93 + (t106 * t163 + t149) * t91 + t202 * t101 + (-t115 * t157 - t146) * mrSges(3,2) + (-t150 - t147) * mrSges(3,1) + (t42 * t36 - t43 * t37 - t144) * mrSges(6,3) + t210; (t87 + t166) * mrSges(4,3) + m(5) * (-t124 * t198 + t137) + t208 * t54 + t207 * t53 + ((t160 * mrSges(3,1) - t115 * mrSges(4,2) - t119 * t91 - t122 * t93) * t120 - m(5) * (t120 * t197 + t172) + (t160 * mrSges(3,2) + t201) * t123 + (-pkin(2) * t167 - t88 * t120 - t172) * m(4)) * t179 + t102 * t39 + t108 * t10 + qJ(3) * t76 + qJD(3) * t80 + (t51 * t36 - t52 * t37 - t144) * mrSges(6,3) + m(4) * t137 - t124 * t129 + (t102 * t75 + t108 * t64 + t207 * t23 + t208 * t22 + t3 * t52 + t4 * t51) * m(6) + t210; m(4) * t147 + t47 * t53 - t48 * t54 - t129 + (-t205 + t206) * mrSges(6,3) - m(5) * t198 + m(6) * (t144 - t203) + ((-m(4) - m(5)) * t92 + t201) * t115; (m(6) * (t118 * t3 + t121 * t4 + t23 * t161 - t22 * t162) - t54 * t162 + t53 * t161 + (-t118 * t37 + t121 * t36) * mrSges(6,3)) * pkin(4) + (-t92 * t142 + t119 * t72 / 0.2e1 + (-t133 / 0.2e1 + t134 / 0.2e1) * t115 + (-t181 / 0.2e1 - t180 / 0.2e1) * qJD(4) + (t71 / 0.2e1 + t204 * pkin(4)) * t122) * t115 - m(6) * (t22 * t24 + t23 * t25) + t138 * t79 - t25 * t53 - t24 * t54 + t44 * mrSges(5,1) - t45 * mrSges(5,2) + t23 * t187 + t126; -t22 * t53 + (t54 + t187) * t23 + t126;];
tauc = t1(:);
