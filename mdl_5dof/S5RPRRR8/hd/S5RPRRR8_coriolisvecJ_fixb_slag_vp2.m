% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:43
% EndTime: 2019-12-31 19:05:49
% DurationCPUTime: 2.68s
% Computational Cost: add. (2894->261), mult. (4708->384), div. (0->0), fcn. (2457->6), ass. (0->150)
t199 = -qJD(3) + qJD(1);
t114 = cos(qJ(3));
t198 = t114 * t199;
t105 = qJD(4) + qJD(5);
t111 = sin(qJ(3));
t112 = cos(qJ(5));
t113 = cos(qJ(4));
t109 = sin(qJ(5));
t110 = sin(qJ(4));
t157 = t109 * t110;
t78 = -t112 * t113 + t157;
t75 = t78 * t111;
t79 = t109 * t113 + t110 * t112;
t197 = t105 * t75 + t79 * t198;
t195 = t105 * t79;
t196 = -t111 * t195 + t78 * t198;
t154 = qJ(2) * qJD(1);
t115 = -pkin(1) - pkin(2);
t95 = qJD(1) * t115 + qJD(2);
t73 = t111 * t95 + t114 * t154;
t57 = -pkin(7) * t199 + t73;
t140 = (t110 ^ 2 + t113 ^ 2) * t57;
t72 = -t111 * t154 + t114 * t95;
t56 = pkin(3) * t199 - t72;
t194 = t111 * t56 + t114 * t140;
t193 = (m(3) * qJ(2) + mrSges(3,3)) * qJD(1);
t66 = t79 * t199;
t192 = t66 / 0.2e1;
t191 = -t105 / 0.2e1;
t190 = t199 / 0.2e1;
t179 = -pkin(8) - pkin(7);
t93 = t179 * t110;
t94 = t179 * t113;
t54 = t109 * t93 - t112 * t94;
t144 = qJD(4) * t179;
t86 = t110 * t144;
t87 = t113 * t144;
t189 = -qJD(5) * t54 - t109 * t86 + t112 * t87 + t79 * t72;
t53 = t109 * t94 + t112 * t93;
t188 = qJD(5) * t53 + t109 * t87 + t112 * t86 + t78 * t72;
t169 = Ifges(5,4) * t110;
t91 = t113 * Ifges(5,2) + t169;
t187 = t199 * t91;
t65 = t78 * t199;
t159 = t199 * t110;
t88 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t159;
t158 = t199 * t113;
t89 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t158;
t185 = (-mrSges(4,2) * t199 + t110 * t88 - t113 * t89) * t114;
t139 = -t111 * qJ(2) + t114 * t115;
t155 = t114 * qJ(2) + t111 * t115;
t172 = t66 * Ifges(6,4);
t27 = t65 * Ifges(6,2) + t105 * Ifges(6,6) - t172;
t58 = Ifges(6,4) * t65;
t28 = -t66 * Ifges(6,1) + t105 * Ifges(6,5) + t58;
t37 = t105 * t65;
t38 = t195 * t199;
t149 = qJD(4) * t110;
t145 = pkin(4) * t149;
t150 = qJD(3) * t114;
t151 = qJD(3) * t111;
t153 = qJD(2) * t111;
t52 = t95 * t151 + (qJ(2) * t150 + t153) * qJD(1);
t42 = -t145 * t199 + t52;
t174 = pkin(4) * t113;
t101 = -pkin(3) - t174;
t45 = -t101 * t199 - t72;
t146 = qJD(5) * t112;
t148 = qJD(4) * t113;
t47 = t105 * t157 - t112 * t148 - t113 * t146;
t152 = qJD(2) * t114;
t51 = t95 * t150 + (-qJ(2) * t151 + t152) * qJD(1);
t80 = (mrSges(5,1) * t110 + mrSges(5,2) * t113) * qJD(4);
t90 = -mrSges(5,1) * t113 + mrSges(5,2) * t110;
t184 = t56 * t80 - t47 * t28 / 0.2e1 - t195 * t27 / 0.2e1 - t51 * mrSges(4,2) - t38 * Ifges(6,2) * t78 + t37 * Ifges(6,1) * t79 + (t90 - mrSges(4,1)) * t52 + t42 * (mrSges(6,1) * t78 + mrSges(6,2) * t79) + t45 * (mrSges(6,1) * t195 - mrSges(6,2) * t47) + (Ifges(6,5) * t47 + Ifges(6,6) * t195) * t191 + (-t37 * t78 + t38 * t79) * Ifges(6,4);
t160 = Ifges(5,6) * qJD(4);
t63 = t160 - t187;
t183 = t63 / 0.2e1;
t161 = Ifges(5,5) * qJD(4);
t96 = Ifges(5,4) * t158;
t64 = -Ifges(5,1) * t159 + t161 - t96;
t182 = -t64 / 0.2e1;
t181 = -t65 / 0.2e1;
t180 = -t66 / 0.2e1;
t178 = m(6) * t45;
t177 = -t110 / 0.2e1;
t85 = -pkin(7) + t155;
t176 = pkin(8) - t85;
t175 = mrSges(6,3) * t65;
t173 = t66 * mrSges(6,3);
t40 = -mrSges(6,1) * t65 - mrSges(6,2) * t66;
t76 = t90 * t199;
t171 = -t40 + t76;
t168 = Ifges(5,4) * t113;
t142 = -pkin(8) * t199 + t57;
t44 = t142 * t113;
t167 = t109 * t44;
t166 = t110 * t51;
t164 = t112 * t44;
t163 = t114 * t52;
t147 = qJD(5) * t109;
t141 = qJD(4) * t176;
t46 = t113 * t51;
t24 = -t149 * t57 + t46;
t82 = (-Ifges(5,2) * t110 + t168) * qJD(4);
t138 = -t24 * mrSges(5,3) + t190 * t82;
t25 = -t148 * t57 - t166;
t83 = (Ifges(5,1) * t113 - t169) * qJD(4);
t137 = t25 * mrSges(5,3) + t190 * t83;
t136 = -mrSges(4,1) * t199 + t171;
t84 = pkin(3) - t139;
t133 = Ifges(6,1) * t47 + Ifges(6,4) * t195;
t131 = Ifges(6,4) * t47 + Ifges(6,2) * t195;
t43 = t142 * t110;
t128 = qJD(4) * t142;
t41 = qJD(4) * pkin(4) - t43;
t10 = t112 * t41 - t167;
t11 = t109 * t41 + t164;
t61 = t176 * t110;
t62 = t176 * t113;
t31 = t109 * t62 + t112 * t61;
t32 = t109 * t61 - t112 * t62;
t127 = -t25 * t110 + t24 * t113;
t126 = t110 * t89 + t113 * t88;
t125 = -t111 * t72 + t114 * t73;
t124 = t83 * t177 - t113 * t82 / 0.2e1;
t16 = -t110 * t128 + t46;
t17 = -t113 * t128 - t166;
t2 = qJD(5) * t10 + t109 * t17 + t112 * t16;
t3 = -qJD(5) * t11 - t109 * t16 + t112 * t17;
t120 = -t10 * t47 + t11 * t195 + t2 * t78 + t3 * t79;
t71 = qJD(3) * t155 + t153;
t119 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + t10 * t175 + t27 * t180 - t45 * (-mrSges(6,1) * t66 + mrSges(6,2) * t65) + (Ifges(6,1) * t65 + t172) * t192 + (Ifges(6,5) * t65 + Ifges(6,6) * t66) * t191 + Ifges(6,6) * t38 + Ifges(6,5) * t37 + (Ifges(6,2) * t66 + t28 + t58) * t181;
t92 = t110 * Ifges(5,1) + t168;
t81 = qJD(4) * (Ifges(5,5) * t113 - Ifges(5,6) * t110);
t77 = t84 + t174;
t74 = t79 * t111;
t70 = qJD(3) * t139 + t152;
t67 = t199 * t80;
t55 = t71 - t145;
t50 = mrSges(6,1) * t105 + t173;
t49 = -mrSges(6,2) * t105 + t175;
t36 = -t110 * t70 + t113 * t141;
t35 = t110 * t141 + t113 * t70;
t13 = -t112 * t43 - t167;
t12 = t109 * t43 - t164;
t7 = -mrSges(6,1) * t38 + mrSges(6,2) * t37;
t5 = -qJD(5) * t32 - t109 * t35 + t112 * t36;
t4 = qJD(5) * t31 + t109 * t36 + t112 * t35;
t1 = [-t184 + t133 * t180 + 0.2e1 * qJD(2) * t193 + m(6) * (t10 * t5 + t11 * t4 + t2 * t32 + t3 * t31 + t42 * t77 + t45 * t55) - t84 * t67 - t71 * t76 + t77 * t7 + t5 * t50 + t55 * t40 + t4 * t49 + (-t81 / 0.2e1 + (t190 * t92 - t85 * t88 + t182) * t113 + (t183 - t187 / 0.2e1 - t85 * t89) * t110) * qJD(4) + m(4) * (-t139 * t52 + t155 * t51 + t73 * t70 - t72 * t71) + (-t70 * t88 + t137) * t110 + (t70 * t89 + t138) * t113 + m(5) * (t127 * t85 + t140 * t70 + t52 * t84 + t56 * t71) + t65 * t131 / 0.2e1 - (-t71 * mrSges(4,1) - t70 * mrSges(4,2) + t124) * t199 + (-t31 * t37 + t32 * t38 + t120) * mrSges(6,3); t197 * t50 + t196 * t49 + (t67 - t7) * t114 - t126 * t111 * qJD(4) + (t37 * t74 - t38 * t75) * mrSges(6,3) + (-t136 * t111 - t185) * qJD(3) + m(5) * (qJD(3) * t194 + t127 * t111 - t163) + m(4) * (qJD(3) * t125 + t111 * t51 - t163) + (t185 - m(5) * t194 - m(4) * t125 + (t136 - t178) * t111 - t193) * qJD(1) + (t10 * t197 + t196 * t11 - t114 * t42 + t45 * t151 - t2 * t75 - t3 * t74) * m(6); t184 + t131 * t181 + t101 * t7 + pkin(3) * t67 + t189 * t50 + t188 * t49 + (t101 * t42 + t2 * t54 + t3 * t53 + (t145 - t73) * t45 + t188 * t11 + t189 * t10) * m(6) - (t72 * mrSges(4,2) + t73 * mrSges(4,1) + (t113 * t92 / 0.2e1 + t91 * t177) * qJD(4) - t124) * t199 + t171 * t73 + (t72 * t88 - t137) * t110 + (-t72 * t89 - t138) * t113 + t133 * t192 + (-t37 * t53 + t38 * t54 - t120) * mrSges(6,3) + (-pkin(3) * t52 + pkin(7) * t127 - t140 * t72 - t56 * t73) * m(5) + (t81 / 0.2e1 + (t64 / 0.2e1 - pkin(7) * t88) * t113 + (-t63 / 0.2e1 - pkin(7) * t89 + pkin(4) * t40) * t110) * qJD(4); -m(6) * (t10 * t12 + t11 * t13) + t119 + t126 * t57 - t13 * t49 - t12 * t50 - t24 * mrSges(5,2) + t25 * mrSges(5,1) - t11 * t173 + (m(6) * (-t10 * t147 + t109 * t2 + t11 * t146 + t112 * t3) + t49 * t146 - t50 * t147 + (t109 * t38 - t112 * t37) * mrSges(6,3)) * pkin(4) - ((-t56 * mrSges(5,2) + t161 / 0.2e1 + t96 / 0.2e1 + t182) * t113 + (-t56 * mrSges(5,1) - t160 / 0.2e1 + t183 - (t169 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t113) * t199 + (-t40 - t178) * pkin(4)) * t110) * t199; t119 + (t50 - t173) * t11 - t10 * t49;];
tauc = t1(:);
