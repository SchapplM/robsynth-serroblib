% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR15_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR15_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR15_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:38
% EndTime: 2019-12-31 18:36:45
% DurationCPUTime: 2.80s
% Computational Cost: add. (2434->358), mult. (5638->534), div. (0->0), fcn. (3454->6), ass. (0->172)
t223 = -Ifges(5,3) / 0.2e1;
t222 = qJD(1) / 0.2e1;
t128 = sin(qJ(5));
t130 = cos(qJ(5));
t126 = sin(pkin(8));
t127 = cos(pkin(8));
t143 = t126 * t128 - t127 * t130;
t131 = cos(qJ(3));
t129 = sin(qJ(3));
t168 = t127 * t129;
t160 = pkin(7) * t168;
t138 = (pkin(4) * t131 + t160) * qJD(1);
t147 = pkin(3) * t131 + qJ(4) * t129;
t112 = t147 * qJD(1);
t132 = -pkin(1) - pkin(6);
t122 = qJD(1) * t132 + qJD(2);
t169 = t126 * t131;
t59 = t127 * t112 - t122 * t169;
t36 = t138 + t59;
t165 = qJD(1) * t129;
t158 = t126 * t165;
t167 = t127 * t131;
t60 = t126 * t112 + t122 * t167;
t42 = pkin(7) * t158 + t60;
t189 = pkin(7) + qJ(4);
t116 = t189 * t126;
t117 = t189 * t127;
t65 = -t116 * t130 - t117 * t128;
t221 = -qJD(4) * t143 + qJD(5) * t65 - t128 * t36 - t130 * t42;
t111 = t126 * t130 + t127 * t128;
t66 = -t116 * t128 + t117 * t130;
t220 = -qJD(4) * t111 - qJD(5) * t66 + t128 * t42 - t130 * t36;
t213 = qJD(5) * t143;
t84 = t111 * t131;
t219 = t143 * qJD(1) - qJD(3) * t84 + t129 * t213;
t86 = t143 * t131;
t95 = t111 * qJD(1);
t97 = t111 * qJD(5);
t218 = -qJD(3) * t86 - t129 * t97 - t95;
t164 = qJD(1) * t131;
t108 = t127 * qJD(3) - t126 * t164;
t162 = t126 * qJD(3);
t109 = t127 * t164 + t162;
t153 = t130 * t108 - t109 * t128;
t217 = (qJ(2) * (m(3) + m(4)));
t216 = t143 * t129;
t87 = qJD(3) * t147 - qJD(4) * t131 + qJD(2);
t67 = t87 * qJD(1);
t170 = t122 * t131;
t90 = (qJD(4) + t170) * qJD(3);
t29 = -t126 * t90 + t127 * t67;
t20 = qJD(3) * t138 + t29;
t152 = pkin(7) * t129 * t162;
t30 = t126 * t67 + t127 * t90;
t27 = qJD(1) * t152 + t30;
t115 = pkin(3) * t129 - qJ(4) * t131 + qJ(2);
t100 = t115 * qJD(1);
t114 = t129 * t122;
t103 = qJD(3) * qJ(4) + t114;
t43 = t127 * t100 - t103 * t126;
t26 = pkin(4) * t165 - pkin(7) * t109 + t43;
t44 = t126 * t100 + t127 * t103;
t28 = pkin(7) * t108 + t44;
t8 = -t128 * t28 + t130 * t26;
t1 = qJD(5) * t8 + t128 * t20 + t130 * t27;
t9 = t128 * t26 + t130 * t28;
t2 = -qJD(5) * t9 - t128 * t27 + t130 * t20;
t215 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t214 = qJD(3) * t131;
t145 = t44 * t126 + t43 * t127;
t182 = Ifges(5,2) * t126;
t184 = Ifges(5,4) * t127;
t148 = t182 - t184;
t185 = Ifges(5,4) * t126;
t149 = -Ifges(5,1) * t127 + t185;
t150 = mrSges(5,1) * t126 + mrSges(5,2) * t127;
t174 = Ifges(4,5) * qJD(3);
t194 = t127 / 0.2e1;
t196 = -t126 / 0.2e1;
t92 = -qJD(3) * pkin(3) + qJD(4) - t170;
t212 = -t145 * mrSges(5,3) - t109 * t149 / 0.2e1 - t108 * t148 / 0.2e1 + (t109 * Ifges(5,1) + t108 * Ifges(5,4) + Ifges(5,5) * t165) * t194 + (t109 * Ifges(5,4) + t108 * Ifges(5,2) + Ifges(5,6) * t165) * t196 + t174 / 0.2e1 + (Ifges(4,1) * t131 - Ifges(4,4) * t129) * t222 + t92 * t150;
t136 = qJD(3) * t216;
t24 = qJD(1) * t136 + qJD(5) * t153;
t211 = t24 / 0.2e1;
t163 = qJD(3) * t129;
t137 = t111 * t163;
t53 = t108 * t128 + t109 * t130;
t25 = qJD(1) * t137 - qJD(5) * t53;
t210 = t25 / 0.2e1;
t209 = -t153 / 0.2e1;
t208 = t153 / 0.2e1;
t207 = -t53 / 0.2e1;
t206 = t53 / 0.2e1;
t205 = -t84 / 0.2e1;
t204 = -t86 / 0.2e1;
t200 = -t143 / 0.2e1;
t199 = t111 / 0.2e1;
t124 = qJD(5) + t165;
t198 = -t124 / 0.2e1;
t197 = t124 / 0.2e1;
t195 = t126 / 0.2e1;
t193 = Ifges(6,4) * t53;
t192 = pkin(4) * t126;
t75 = t129 * t95;
t188 = t75 + t97;
t76 = qJD(1) * t216;
t187 = t76 + t213;
t186 = Ifges(4,4) * t131;
t183 = Ifges(5,5) * t127;
t181 = Ifges(5,6) * t126;
t180 = qJ(2) * mrSges(4,1);
t179 = qJ(2) * mrSges(4,2);
t157 = t132 * t214;
t55 = t126 * t87 + t127 * t157;
t175 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t108 + mrSges(5,2) * t109 + mrSges(4,3) * t164;
t173 = Ifges(4,6) * qJD(3);
t172 = qJD(3) * mrSges(4,2);
t166 = t129 * t132;
t69 = t126 * t115 + t127 * t166;
t161 = qJD(1) * qJD(3);
t156 = t131 * t161;
t159 = Ifges(6,5) * t24 + Ifges(6,6) * t25 + Ifges(6,3) * t156;
t7 = -t25 * mrSges(6,1) + t24 * mrSges(6,2);
t155 = -t126 * t132 + pkin(4);
t154 = -t132 + t192;
t151 = mrSges(4,1) * t129 + mrSges(4,2) * t131;
t146 = -t126 * t29 + t127 * t30;
t144 = -t126 * t43 + t127 * t44;
t106 = t127 * t115;
t49 = -pkin(7) * t167 + t129 * t155 + t106;
t58 = -pkin(7) * t169 + t69;
t15 = -t128 * t58 + t130 * t49;
t16 = t128 * t49 + t130 * t58;
t142 = t183 / 0.2e1 - t181 / 0.2e1;
t141 = t150 * qJD(1);
t135 = -t109 * Ifges(5,5) - t108 * Ifges(5,6) - t43 * mrSges(5,1) + t44 * mrSges(5,2) - t124 * Ifges(6,3) + t9 * mrSges(6,2) - t153 * Ifges(6,6) - t53 * Ifges(6,5) - t8 * mrSges(6,1) + t173 / 0.2e1 + (-Ifges(4,2) * t129 + t186) * t222 + t165 * t223;
t125 = -pkin(4) * t127 - pkin(3);
t119 = -mrSges(4,3) * t165 - t172;
t113 = t151 * qJD(1);
t107 = t154 * t131;
t93 = t154 * t163;
t89 = (mrSges(5,1) * t131 + mrSges(5,3) * t168) * t161;
t88 = (mrSges(5,3) * t126 * t129 - mrSges(5,2) * t131) * t161;
t83 = t111 * t129;
t78 = -pkin(4) * t158 + t114;
t77 = t141 * t163;
t74 = mrSges(5,1) * t165 - mrSges(5,3) * t109;
t73 = -mrSges(5,2) * t165 + mrSges(5,3) * t108;
t72 = t127 * t87;
t70 = (-qJD(1) * t192 + t122) * t163;
t68 = -t126 * t166 + t106;
t64 = (Ifges(5,5) * t131 + t129 * t149) * t161;
t63 = (Ifges(5,6) * t131 + t129 * t148) * t161;
t56 = -pkin(4) * t108 + t92;
t54 = -t126 * t157 + t72;
t48 = Ifges(6,4) * t153;
t41 = t152 + t55;
t40 = t131 * t213 + t137;
t38 = -t131 * t97 + t136;
t35 = mrSges(6,1) * t124 - mrSges(6,3) * t53;
t34 = -mrSges(6,2) * t124 + mrSges(6,3) * t153;
t31 = t72 + (t131 * t155 + t160) * qJD(3);
t19 = -mrSges(6,2) * t156 + mrSges(6,3) * t25;
t18 = mrSges(6,1) * t156 - mrSges(6,3) * t24;
t17 = -mrSges(6,1) * t153 + mrSges(6,2) * t53;
t14 = Ifges(6,1) * t53 + Ifges(6,5) * t124 + t48;
t13 = Ifges(6,2) * t153 + Ifges(6,6) * t124 + t193;
t6 = t24 * Ifges(6,1) + t25 * Ifges(6,4) + Ifges(6,5) * t156;
t5 = t24 * Ifges(6,4) + t25 * Ifges(6,2) + Ifges(6,6) * t156;
t4 = -qJD(5) * t16 - t128 * t41 + t130 * t31;
t3 = qJD(5) * t15 + t128 * t31 + t130 * t41;
t10 = [(Ifges(6,4) * t38 + Ifges(6,2) * t40) * t208 + (Ifges(6,5) * t38 + Ifges(6,6) * t40) * t197 + t6 * t204 + t5 * t205 + (Ifges(6,1) * t38 + Ifges(6,4) * t40) * t206 + (t132 * t77 + t63 * t196 + t64 * t194 + (-t126 * t30 - t127 * t29) * mrSges(5,3)) * t131 + qJD(2) * t113 + t107 * t7 + t69 * t88 + t68 * t89 - t93 * t17 + t55 * t73 + t54 * t74 + t56 * (-mrSges(6,1) * t40 + mrSges(6,2) * t38) + t3 * t34 + t4 * t35 + t38 * t14 / 0.2e1 + t40 * t13 / 0.2e1 + t15 * t18 + t16 * t19 + (-t135 + t132 * t119 - t173 / 0.2e1) * t214 + (-t30 * mrSges(5,2) + t29 * mrSges(5,1) + Ifges(6,6) * t210 + Ifges(6,5) * t211 + t159 / 0.2e1 + (-t174 / 0.2e1 + t150 * t170 + (m(5) * (t92 - t170) + t175) * t132 + (-0.2e1 * t179 + (-0.3e1 / 0.2e1 * t183 + 0.3e1 / 0.2e1 * t181 + 0.3e1 / 0.2e1 * Ifges(4,4)) * t129 + (0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2) + Ifges(6,3) / 0.2e1 - Ifges(5,1) * t127 ^ 2 / 0.2e1 + (t184 - t182 / 0.2e1) * t126) * t131) * qJD(1) - t212) * qJD(3) + t215) * t129 + (((2 * mrSges(3,3)) + t151 + (2 * t217)) * qJD(2) + (0.2e1 * t180 + Ifges(6,5) * t204 + Ifges(6,6) * t205 + (-0.3e1 / 0.2e1 * Ifges(4,4) + t142) * t131) * t214) * qJD(1) + t70 * (mrSges(6,1) * t84 - mrSges(6,2) * t86) + (-Ifges(6,4) * t86 - Ifges(6,2) * t84) * t210 + (-Ifges(6,1) * t86 - Ifges(6,4) * t84) * t211 + (-t1 * t84 + t2 * t86 - t38 * t8 + t40 * t9) * mrSges(6,3) + m(6) * (t1 * t16 + t107 * t70 + t15 * t2 + t3 * t9 + t4 * t8 - t56 * t93) + m(5) * (t29 * t68 + t30 * t69 + t43 * t54 + t44 * t55); -t83 * t18 - t216 * t19 + t219 * t35 + t218 * t34 + (-t7 + t77) * t131 + (-t126 * t89 + t127 * t88) * t129 + ((-t126 * t74 + t127 * t73 + t119) * t131 + (t17 + t175) * t129) * qJD(3) + m(5) * (t146 * t129 + (t129 * t92 + (t144 - t114) * t131) * qJD(3)) + (-t1 * t216 - t131 * t70 + t56 * t163 - t2 * t83 + t218 * t9 + t219 * t8) * m(6) + (-m(5) * t145 - t126 * t73 - t127 * t74 - t113 + (-mrSges(3,3) - t217) * qJD(1)) * qJD(1); (Ifges(6,4) * t76 + Ifges(6,2) * t75) * t209 + (Ifges(6,5) * t76 + Ifges(6,6) * t75) * t198 + t6 * t199 + t5 * t200 + (Ifges(6,1) * t76 + Ifges(6,4) * t75) * t207 + (-t213 / 0.2e1 - t76 / 0.2e1) * t14 + (-Ifges(6,4) * t213 - Ifges(6,2) * t97) * t208 + (-Ifges(6,5) * t213 - Ifges(6,6) * t97) * t197 + (-Ifges(6,1) * t213 - Ifges(6,4) * t97) * t206 + (Ifges(6,4) * t111 - Ifges(6,2) * t143) * t210 + (Ifges(6,1) * t111 - Ifges(6,4) * t143) * t211 + (-t1 * t143 - t111 * t2 + t187 * t8 - t188 * t9) * mrSges(6,3) + t70 * (mrSges(6,1) * t143 + mrSges(6,2) * t111) + (mrSges(6,1) * t188 - mrSges(6,2) * t187) * t56 + ((-t119 - t172) * t131 + ((-mrSges(5,1) * t127 + mrSges(5,2) * t126 - mrSges(4,1)) * qJD(3) - t175) * t129) * t122 + t125 * t7 + pkin(3) * t77 - t78 * t17 + t65 * t18 + t66 * t19 - t60 * t73 - t59 * t74 + (qJD(4) * t73 + qJ(4) * t88 + t30 * mrSges(5,3) + t63 / 0.2e1) * t127 + (-qJD(4) * t74 - qJ(4) * t89 - t29 * mrSges(5,3) + t64 / 0.2e1) * t126 + (-pkin(3) * t122 * t163 + qJ(4) * t146 + qJD(4) * t144 - t114 * t92 - t43 * t59 - t44 * t60) * m(5) + (-t97 / 0.2e1 - t75 / 0.2e1) * t13 + t220 * t35 + t221 * t34 + (t1 * t66 + t125 * t70 + t2 * t65 + t220 * t8 + t221 * t9 - t56 * t78) * m(6) + ((t135 + (-t180 + t186 / 0.2e1) * qJD(1)) * t131 + ((t179 + (-Ifges(4,4) / 0.2e1 + t142) * t129 + (t223 + Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1) * t131) * qJD(1) + t212) * t129 + ((-Ifges(4,6) / 0.2e1 + Ifges(5,5) * t195 + Ifges(5,6) * t194 + Ifges(6,5) * t199 + Ifges(6,6) * t200) * t131 + (-Ifges(4,5) / 0.2e1 - t127 * (Ifges(5,1) * t126 + t184) / 0.2e1 + (Ifges(5,2) * t127 + t185) * t195) * t129) * qJD(3)) * qJD(1); -t108 * t73 + t109 * t74 - t153 * t34 + t53 * t35 + (m(5) * t122 - t141) * t163 - m(5) * (t108 * t44 - t109 * t43) + t7 + (-t153 * t9 + t53 * t8 + t70) * m(6); -t56 * (mrSges(6,1) * t53 + mrSges(6,2) * t153) + (Ifges(6,1) * t153 - t193) * t207 + t13 * t206 + (Ifges(6,5) * t153 - Ifges(6,6) * t53) * t198 - t8 * t34 + t9 * t35 + (t153 * t8 + t53 * t9) * mrSges(6,3) + t159 + (-Ifges(6,2) * t53 + t14 + t48) * t209 + t215;];
tauc = t10(:);
