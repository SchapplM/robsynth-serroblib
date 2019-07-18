% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:21
% EndTime: 2019-07-18 13:28:28
% DurationCPUTime: 2.53s
% Computational Cost: add. (2005->305), mult. (5755->441), div. (0->0), fcn. (4226->8), ass. (0->153)
t107 = qJD(3) + qJD(4);
t221 = t107 * Ifges(5,6) / 0.2e1;
t117 = cos(qJ(2));
t156 = t117 * qJD(1);
t116 = cos(qJ(3));
t161 = qJD(2) * t116;
t100 = -pkin(2) * t161 - t156;
t177 = t107 * Ifges(5,5);
t220 = t100 * mrSges(5,2) + t177 / 0.2e1;
t111 = sin(qJ(4));
t112 = sin(qJ(3));
t115 = cos(qJ(4));
t130 = t111 * t116 + t112 * t115;
t94 = t130 * qJD(2);
t183 = t94 * Ifges(5,4);
t166 = t111 * t112;
t129 = -t115 * t116 + t166;
t93 = t129 * qJD(2);
t219 = t221 + t183 / 0.2e1 - t93 * Ifges(5,2) / 0.2e1;
t110 = sin(qJ(5));
t114 = cos(qJ(5));
t113 = sin(qJ(2));
t164 = qJD(1) * t113;
t178 = qJD(3) * pkin(2);
t126 = -t112 * t164 + t178;
t152 = t116 * t164;
t76 = t111 * t152 - t115 * t126;
t168 = qJD(4) * t76;
t159 = qJD(3) * t113;
t160 = qJD(2) * t117;
t212 = -t112 * t159 + t116 * t160;
t85 = t212 * qJD(1);
t151 = t112 * t160;
t86 = (-t116 * t159 - t151) * qJD(1);
t29 = t111 * t86 + t115 * t85 - t168;
t77 = t111 * t126 + t115 * t152;
t45 = t100 * t114 - t110 * t77;
t96 = (t112 * t178 + t164) * qJD(2);
t7 = qJD(5) * t45 + t110 * t96 + t114 * t29;
t46 = t100 * t110 + t114 * t77;
t8 = -qJD(5) * t46 - t110 * t29 + t114 * t96;
t216 = -t8 * t110 + t114 * t7;
t134 = t110 * t46 + t114 * t45;
t137 = Ifges(6,5) * t114 - Ifges(6,6) * t110;
t179 = Ifges(6,4) * t114;
t139 = -Ifges(6,2) * t110 + t179;
t180 = Ifges(6,4) * t110;
t141 = Ifges(6,1) * t114 - t180;
t197 = t114 / 0.2e1;
t200 = -t110 / 0.2e1;
t75 = t107 * t110 + t114 * t94;
t204 = t75 / 0.2e1;
t194 = Ifges(6,4) * t75;
t74 = t107 * t114 - t110 * t94;
t90 = qJD(5) + t93;
t26 = Ifges(6,2) * t74 + Ifges(6,6) * t90 + t194;
t71 = Ifges(6,4) * t74;
t27 = Ifges(6,1) * t75 + Ifges(6,5) * t90 + t71;
t215 = -t134 * mrSges(6,3) + t90 * t137 / 0.2e1 + t141 * t204 + t74 * t139 / 0.2e1 + t27 * t197 + t26 * t200;
t214 = -t100 * mrSges(5,1) - t45 * mrSges(6,1) + t46 * mrSges(6,2) + t219;
t213 = -t26 / 0.2e1;
t87 = t130 * t113;
t211 = -t110 * t45 + t114 * t46;
t121 = t107 * t129;
t60 = t121 * qJD(2);
t33 = t74 * qJD(5) - t114 * t60;
t34 = -t75 * qJD(5) + t110 * t60;
t210 = t8 * mrSges(6,1) - t7 * mrSges(6,2) + Ifges(6,5) * t33 + Ifges(6,6) * t34;
t209 = t33 / 0.2e1;
t208 = t34 / 0.2e1;
t122 = t107 * t130;
t61 = t122 * qJD(2);
t207 = t61 / 0.2e1;
t206 = -t74 / 0.2e1;
t205 = -t75 / 0.2e1;
t203 = -t90 / 0.2e1;
t199 = -t112 / 0.2e1;
t198 = t112 / 0.2e1;
t195 = mrSges(5,3) * t93;
t89 = Ifges(5,4) * t93;
t192 = t74 * Ifges(6,6);
t191 = t75 * Ifges(6,5);
t81 = t129 * t164;
t190 = t76 * t81;
t82 = t130 * t156;
t189 = t76 * t82;
t187 = t90 * Ifges(6,3);
t185 = t94 * mrSges(5,3);
t184 = t94 * Ifges(5,1);
t182 = -mrSges(5,1) * t107 - mrSges(6,1) * t74 + mrSges(6,2) * t75 + t185;
t181 = Ifges(4,4) * t112;
t174 = t110 * t93;
t172 = t114 * t93;
t171 = t116 * Ifges(4,2);
t170 = Ifges(4,5) * qJD(3);
t169 = Ifges(4,6) * qJD(3);
t167 = qJD(4) * t77;
t165 = t112 ^ 2 + t116 ^ 2;
t163 = qJD(2) * t112;
t162 = qJD(2) * t113;
t158 = qJD(5) * t110;
t157 = qJD(5) * t114;
t155 = pkin(2) * t163;
t63 = mrSges(5,1) * t93 + mrSges(5,2) * t94;
t148 = m(5) * t100 + t63;
t146 = -t110 * t7 - t114 * t8;
t30 = t111 * t85 - t115 * t86 + t167;
t36 = -qJD(4) * t113 * t166 + (t107 * t116 * t113 + t151) * t115 + t212 * t111;
t145 = t30 * t87 + t36 * t76;
t144 = mrSges(4,1) * t112 + mrSges(4,2) * t116;
t143 = mrSges(6,1) * t114 - mrSges(6,2) * t110;
t142 = mrSges(6,1) * t110 + mrSges(6,2) * t114;
t140 = Ifges(6,1) * t110 + t179;
t138 = Ifges(6,2) * t114 + t180;
t136 = Ifges(6,5) * t110 + Ifges(6,6) * t114;
t42 = -mrSges(6,2) * t90 + mrSges(6,3) * t74;
t43 = mrSges(6,1) * t90 - mrSges(6,3) * t75;
t135 = -t110 * t42 - t114 * t43;
t88 = t129 * t113;
t72 = t110 * t88 - t114 * t117;
t132 = t110 * t117 + t114 * t88;
t103 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t163;
t104 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t161;
t131 = t103 * t116 + t104 * t112;
t128 = mrSges(5,3) + t142;
t78 = -mrSges(5,2) * t107 - t195;
t127 = -t110 * t43 + t114 * t42 + t78;
t25 = t187 + t191 + t192;
t3 = t33 * Ifges(6,4) + t34 * Ifges(6,2) + t61 * Ifges(6,6);
t4 = t33 * Ifges(6,1) + t34 * Ifges(6,4) + t61 * Ifges(6,5);
t52 = t177 - t89 + t184;
t120 = t174 * t213 + t110 * t4 / 0.2e1 - Ifges(5,5) * t60 - Ifges(5,6) * t61 - t29 * mrSges(5,2) + t136 * t207 + t138 * t208 + t140 * t209 + t3 * t197 + t76 * t195 + t27 * t172 / 0.2e1 + (-t137 * t203 - t139 * t206 - t141 * t205 + t220) * t93 + (Ifges(6,5) * t205 + Ifges(6,6) * t206 + Ifges(6,3) * t203 + t214 + t221) * t94 + (-t143 - mrSges(5,1)) * t30 + (-Ifges(5,2) * t94 + t52 - t89) * t93 / 0.2e1 - (-Ifges(5,1) * t93 - t183 + t25) * t94 / 0.2e1 + (-t45 * t172 - t46 * t174 + t216) * mrSges(6,3) + (t142 * t76 + t215) * qJD(5);
t118 = qJD(2) ^ 2;
t106 = Ifges(4,4) * t161;
t99 = (-mrSges(4,1) * t116 + mrSges(4,2) * t112) * qJD(2);
t92 = Ifges(4,1) * t163 + t106 + t170;
t91 = t169 + (t171 + t181) * qJD(2);
t84 = t129 * t156;
t83 = t130 * t164;
t67 = t110 * t164 - t114 * t84;
t66 = t110 * t84 + t114 * t164;
t65 = t110 * t155 - t114 * t83;
t64 = t110 * t83 + t114 * t155;
t57 = Ifges(6,3) * t61;
t53 = t142 * t93;
t35 = -t107 * t87 - t117 * t93;
t24 = mrSges(5,1) * t61 - mrSges(5,2) * t60;
t12 = t132 * qJD(5) - t110 * t35 + t114 * t162;
t11 = t72 * qJD(5) + t110 * t162 + t114 * t35;
t10 = -mrSges(6,2) * t61 + mrSges(6,3) * t34;
t9 = mrSges(6,1) * t61 - mrSges(6,3) * t33;
t6 = -mrSges(6,1) * t34 + mrSges(6,2) * t33;
t1 = [-t132 * t10 + t11 * t42 + t12 * t43 + t35 * t78 + t87 * t6 + t72 * t9 + t182 * t36 + (-t118 * mrSges(3,2) - t24) * t117 + (-t60 * t87 + t61 * t88) * mrSges(5,3) + m(5) * (-t117 * t96 - t29 * t88 + t35 * t77 + t145) + m(6) * (t11 * t46 + t12 * t45 - t132 * t7 + t72 * t8 + t145) + (-t144 * qJD(3) - t112 * t103 + t116 * t104) * t160 + (-t118 * mrSges(3,1) + m(4) * (-t112 * t86 + t116 * t85) - t131 * qJD(3) + (t99 + m(4) * (-0.2e1 + t165) * t156 + t148) * qJD(2)) * t113; -t67 * t42 - t66 * t43 + t84 * t78 - t182 * t82 + (-t63 - t99) * t164 + (t187 / 0.2e1 + t192 / 0.2e1 + t191 / 0.2e1 + t25 / 0.2e1 - t77 * mrSges(5,3) - t214 - t219) * t122 - (-t89 / 0.2e1 + t184 / 0.2e1 + t52 / 0.2e1 + t128 * t76 + t215 + t220) * t121 - m(5) * (t100 * t164 - t77 * t84 + t189) - m(6) * (t45 * t66 + t46 * t67 + t189) - m(4) * (-0.1e1 + t165) * qJD(1) ^ 2 * t113 * t117 + (Ifges(5,4) * t60 + t96 * mrSges(5,1) + t57 / 0.2e1 - t29 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t61 + t210) * t129 + (t137 * t207 + t141 * t209 + t139 * t208 - Ifges(5,1) * t60 - Ifges(5,4) * t61 + t96 * mrSges(5,2) + t4 * t197 + t3 * t200 + t128 * t30 + t146 * mrSges(6,3) + (-mrSges(6,3) * t211 + t114 * t213 + t136 * t203 + t138 * t206 + t140 * t205 + t76 * t143 + t27 * t200) * qJD(5)) * t130 + (-t86 * mrSges(4,3) + (mrSges(4,2) * t162 + t103 * t117) * qJD(1) + (-mrSges(4,1) * t156 - t91 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,4) * t163 - t169 / 0.2e1 + (m(6) * t134 - t135 + t148) * pkin(2)) * qJD(3)) * t112 + (t85 * mrSges(4,3) + (-mrSges(4,1) * t162 - t104 * t117) * qJD(1) + (t92 / 0.2e1 - mrSges(4,2) * t156 + t170 / 0.2e1 + (0.3e1 / 0.2e1 * Ifges(4,4) * t116 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t112) * qJD(2)) * qJD(3) + (-t114 * t9 - t110 * t10 - m(5) * t96 + m(6) * (-t46 * t157 + t45 * t158 + t146) - t24 + t43 * t158 - t42 * t157) * pkin(2)) * t116; (-t148 * t163 + (t60 * mrSges(5,3) - t6 - m(6) * t30 + m(5) * (-t30 + t167) + (m(6) * t211 + t127) * qJD(4)) * t115 + (-t61 * mrSges(5,3) + t114 * t10 - t110 * t9 + t135 * qJD(5) + t182 * qJD(4) + m(6) * (-t45 * t157 - t46 * t158 + t168 + t216) + m(5) * (t29 + t168)) * t111) * pkin(2) - m(5) * (-t77 * t83 - t190) - t85 * mrSges(4,2) + t86 * mrSges(4,1) + t83 * t78 + t76 * t53 - t64 * t43 - t65 * t42 + t77 * t185 + t120 - m(6) * (t45 * t64 + t46 * t65 - t190) + (t144 * t156 + t91 * t198 + (t171 * t198 + (Ifges(4,1) * t116 - t181) * t199) * qJD(2) + (Ifges(4,5) * t116 / 0.2e1 + Ifges(4,6) * t199) * qJD(3) - (t106 + t92) * t116 / 0.2e1) * qJD(2) + t182 * t81 + t131 * t164; (-t182 + t185) * t77 + (t53 - m(6) * (-t211 + t77) + t127) * t76 + t120; t57 - t76 * (mrSges(6,1) * t75 + mrSges(6,2) * t74) + (Ifges(6,1) * t74 - t194) * t205 + t26 * t204 + (Ifges(6,5) * t74 - Ifges(6,6) * t75) * t203 - t45 * t42 + t46 * t43 + (t45 * t74 + t46 * t75) * mrSges(6,3) + (-Ifges(6,2) * t75 + t27 + t71) * t206 + t210;];
tauc  = t1(:);
