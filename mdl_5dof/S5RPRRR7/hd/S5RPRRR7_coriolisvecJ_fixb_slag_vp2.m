% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:12
% EndTime: 2019-12-31 19:03:23
% DurationCPUTime: 4.81s
% Computational Cost: add. (3813->402), mult. (9174->567), div. (0->0), fcn. (5691->8), ass. (0->194)
t140 = sin(qJ(4));
t226 = -pkin(8) - pkin(7);
t170 = qJD(4) * t226;
t144 = cos(qJ(3));
t183 = qJD(1) * t144;
t141 = sin(qJ(3));
t165 = pkin(3) * t141 - pkin(7) * t144;
t120 = t165 * qJD(1);
t143 = cos(qJ(4));
t134 = sin(pkin(9)) * pkin(1) + pkin(6);
t124 = t134 * qJD(1);
t94 = qJD(2) * t144 - t141 * t124;
t57 = t140 * t120 + t143 * t94;
t258 = -t57 + (pkin(8) * t183 + t170) * t140;
t186 = t143 * t144;
t152 = pkin(4) * t141 - pkin(8) * t186;
t56 = t143 * t120 - t140 * t94;
t257 = -qJD(1) * t152 + t143 * t170 - t56;
t256 = -Ifges(4,1) / 0.2e1;
t136 = Ifges(4,4) * t183;
t255 = -t136 / 0.2e1;
t254 = qJD(3) / 0.2e1;
t139 = sin(qJ(5));
t142 = cos(qJ(5));
t179 = qJD(3) * t143;
t184 = qJD(1) * t141;
t115 = -t140 * t184 + t179;
t175 = t141 * qJD(2);
t95 = t124 * t144 + t175;
t87 = qJD(3) * pkin(7) + t95;
t171 = -cos(pkin(9)) * pkin(1) - pkin(2);
t111 = -pkin(3) * t144 - t141 * pkin(7) + t171;
t90 = t111 * qJD(1);
t47 = t140 * t90 + t143 * t87;
t35 = pkin(8) * t115 + t47;
t193 = t142 * t35;
t133 = qJD(4) - t183;
t181 = qJD(3) * t140;
t116 = t143 * t184 + t181;
t46 = -t140 * t87 + t143 * t90;
t34 = -pkin(8) * t116 + t46;
t32 = pkin(4) * t133 + t34;
t10 = t139 * t32 + t193;
t180 = qJD(3) * t141;
t167 = qJD(1) * t180;
t131 = Ifges(6,3) * t167;
t166 = t142 * t115 - t116 * t139;
t66 = t115 * t139 + t116 * t142;
t215 = Ifges(6,4) * t66;
t130 = qJD(5) + t133;
t220 = -t130 / 0.2e1;
t232 = -t66 / 0.2e1;
t234 = -t166 / 0.2e1;
t60 = Ifges(6,4) * t166;
t30 = Ifges(6,1) * t66 + Ifges(6,5) * t130 + t60;
t86 = -qJD(3) * pkin(3) - t94;
t59 = -t115 * pkin(4) + t86;
t195 = t139 * t35;
t9 = t142 * t32 - t195;
t253 = t131 + (Ifges(6,5) * t166 - Ifges(6,6) * t66) * t220 + (t10 * t66 + t166 * t9) * mrSges(6,3) + (-Ifges(6,2) * t66 + t30 + t60) * t234 - t59 * (mrSges(6,1) * t66 + mrSges(6,2) * t166) + (Ifges(6,1) * t166 - t215) * t232;
t169 = Ifges(4,5) * t254;
t128 = t226 * t140;
t129 = t226 * t143;
t76 = t128 * t142 + t129 * t139;
t252 = qJD(5) * t76 + t257 * t139 + t258 * t142;
t77 = t128 * t139 - t129 * t142;
t251 = -qJD(5) * t77 - t258 * t139 + t257 * t142;
t156 = t140 * t47 + t143 * t46;
t203 = Ifges(5,4) * t143;
t160 = -Ifges(5,2) * t140 + t203;
t204 = Ifges(5,4) * t140;
t162 = Ifges(5,1) * t143 - t204;
t163 = mrSges(5,1) * t140 + mrSges(5,2) * t143;
t201 = Ifges(5,6) * t140;
t202 = Ifges(5,5) * t143;
t216 = t143 / 0.2e1;
t217 = -t140 / 0.2e1;
t221 = t116 / 0.2e1;
t205 = Ifges(5,4) * t116;
t54 = Ifges(5,2) * t115 + Ifges(5,6) * t133 + t205;
t112 = Ifges(5,4) * t115;
t55 = Ifges(5,1) * t116 + Ifges(5,5) * t133 + t112;
t145 = -t156 * mrSges(5,3) + t133 * (-t201 + t202) / 0.2e1 + t86 * t163 + t115 * t160 / 0.2e1 + t162 * t221 + t54 * t217 + t55 * t216;
t250 = t94 * mrSges(4,3) + t184 * t256 - t145 - t169 + t255;
t123 = t165 * qJD(3);
t110 = qJD(1) * t123;
t88 = t94 * qJD(3);
t20 = -qJD(4) * t47 + t143 * t110 - t140 * t88;
t174 = qJD(3) * qJD(4);
t177 = qJD(4) * t140;
t178 = qJD(3) * t144;
t80 = t143 * t174 + (-t141 * t177 + t143 * t178) * qJD(1);
t13 = pkin(4) * t167 - pkin(8) * t80 + t20;
t176 = qJD(4) * t143;
t19 = t140 * t110 + t143 * t88 + t90 * t176 - t177 * t87;
t148 = t140 * t178 + t141 * t176;
t81 = -qJD(1) * t148 - t140 * t174;
t14 = pkin(8) * t81 + t19;
t2 = qJD(5) * t9 + t13 * t139 + t14 * t142;
t23 = qJD(5) * t166 + t139 * t81 + t142 * t80;
t24 = -qJD(5) * t66 - t139 * t80 + t142 * t81;
t3 = -qJD(5) * t10 + t13 * t142 - t139 * t14;
t249 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t23 + Ifges(6,6) * t24;
t29 = Ifges(6,2) * t166 + Ifges(6,6) * t130 + t215;
t247 = t29 / 0.2e1;
t168 = -Ifges(4,6) * qJD(3) / 0.2e1;
t33 = -mrSges(6,1) * t166 + mrSges(6,2) * t66;
t242 = m(6) * t59 + t33;
t153 = t139 * t140 - t142 * t143;
t99 = t153 * t141;
t119 = t134 * t186;
t68 = t140 * t111 + t119;
t157 = -t140 * t20 + t143 * t19;
t241 = qJD(4) + qJD(5);
t238 = -t20 * mrSges(5,1) + t19 * mrSges(5,2) - Ifges(5,5) * t80 - Ifges(5,6) * t81 - t249;
t172 = mrSges(4,3) * t184;
t191 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t115 + mrSges(5,2) * t116 + t172;
t237 = -m(4) * t94 + m(5) * t86 + t191;
t236 = t23 / 0.2e1;
t235 = t24 / 0.2e1;
t233 = t166 / 0.2e1;
t231 = t66 / 0.2e1;
t230 = t80 / 0.2e1;
t229 = t81 / 0.2e1;
t118 = t139 * t143 + t140 * t142;
t98 = t118 * t141;
t228 = -t98 / 0.2e1;
t227 = -t99 / 0.2e1;
t223 = -t115 / 0.2e1;
t222 = -t116 / 0.2e1;
t219 = t130 / 0.2e1;
t218 = -t133 / 0.2e1;
t214 = pkin(4) * t140;
t213 = pkin(8) * t141;
t69 = t241 * t153;
t150 = t153 * t144;
t92 = qJD(1) * t150;
t208 = -t69 + t92;
t70 = t241 * t118;
t151 = t118 * t144;
t91 = qJD(1) * t151;
t207 = -t70 + t91;
t206 = Ifges(4,4) * t141;
t126 = t171 * qJD(1);
t198 = t126 * mrSges(4,2);
t89 = qJD(3) * t95;
t187 = t134 * t140;
t185 = t143 * t123 + t180 * t187;
t173 = mrSges(4,3) * t183;
t164 = mrSges(5,1) * t143 - mrSges(5,2) * t140;
t161 = Ifges(5,1) * t140 + t203;
t159 = Ifges(5,2) * t143 + t204;
t158 = Ifges(5,5) * t140 + Ifges(5,6) * t143;
t101 = t143 * t111;
t52 = -t143 * t213 + t101 + (-pkin(4) - t187) * t144;
t58 = -t140 * t213 + t68;
t25 = -t139 * t58 + t142 * t52;
t26 = t139 * t52 + t142 * t58;
t61 = mrSges(5,1) * t167 - mrSges(5,3) * t80;
t62 = -mrSges(5,2) * t167 + mrSges(5,3) * t81;
t155 = -t140 * t61 + t143 * t62;
t83 = -mrSges(5,2) * t133 + mrSges(5,3) * t115;
t84 = mrSges(5,1) * t133 - mrSges(5,3) * t116;
t154 = -t140 * t83 - t143 * t84;
t36 = t140 * t123 + t111 * t176 + (-t141 * t179 - t144 * t177) * t134;
t146 = t126 * mrSges(4,1) + t46 * mrSges(5,1) + t9 * mrSges(6,1) + t168 - (Ifges(4,2) * t144 + t206) * qJD(1) / 0.2e1 + t130 * Ifges(6,3) + t66 * Ifges(6,5) + t166 * Ifges(6,6) + t133 * Ifges(5,3) + t116 * Ifges(5,5) + t115 * Ifges(5,6) - t10 * mrSges(6,2) - t47 * mrSges(5,2) - t95 * mrSges(4,3);
t135 = -pkin(4) * t143 - pkin(3);
t132 = Ifges(5,3) * t167;
t127 = -qJD(3) * mrSges(4,2) + t173;
t102 = (t134 + t214) * t141;
t73 = t175 + (qJD(1) * t214 + t124) * t144;
t72 = pkin(4) * t148 + t134 * t178;
t67 = -t144 * t187 + t101;
t51 = mrSges(6,1) * t130 - mrSges(6,3) * t66;
t50 = -mrSges(6,2) * t130 + mrSges(6,3) * t166;
t48 = -t81 * pkin(4) + t89;
t44 = -mrSges(5,1) * t81 + mrSges(5,2) * t80;
t41 = t80 * Ifges(5,1) + t81 * Ifges(5,4) + Ifges(5,5) * t167;
t40 = t80 * Ifges(5,4) + t81 * Ifges(5,2) + Ifges(5,6) * t167;
t39 = -qJD(3) * t151 + t241 * t99;
t38 = -qJD(3) * t150 - t141 * t70;
t37 = -qJD(4) * t68 + t185;
t31 = -pkin(8) * t148 + t36;
t27 = t152 * qJD(3) + (-t119 + (-t111 + t213) * t140) * qJD(4) + t185;
t18 = -mrSges(6,2) * t167 + mrSges(6,3) * t24;
t17 = mrSges(6,1) * t167 - mrSges(6,3) * t23;
t12 = t142 * t34 - t195;
t11 = -t139 * t34 - t193;
t8 = -mrSges(6,1) * t24 + mrSges(6,2) * t23;
t7 = t23 * Ifges(6,1) + t24 * Ifges(6,4) + Ifges(6,5) * t167;
t6 = t23 * Ifges(6,4) + t24 * Ifges(6,2) + Ifges(6,6) * t167;
t5 = -qJD(5) * t26 - t139 * t31 + t142 * t27;
t4 = qJD(5) * t25 + t139 * t27 + t142 * t31;
t1 = [t59 * (-mrSges(6,1) * t39 + mrSges(6,2) * t38) + t4 * t50 + t5 * t51 + t38 * t30 / 0.2e1 + t25 * t17 + t26 * t18 + t39 * t247 + m(6) * (t10 * t4 + t102 * t48 + t2 * t26 + t25 * t3 + t5 * t9 + t59 * t72) + m(5) * (t19 * t68 + t20 * t67 + t47 * t36 + t46 * t37) + (t10 * t39 - t2 * t98 + t3 * t99 - t38 * t9) * mrSges(6,3) + t48 * (mrSges(6,1) * t98 - mrSges(6,2) * t99) + (-Ifges(6,4) * t99 - Ifges(6,2) * t98) * t235 + (-Ifges(6,1) * t99 - Ifges(6,4) * t98) * t236 + (-t131 / 0.2e1 - t132 / 0.2e1 + (m(4) * t134 + mrSges(4,3)) * t88 + (0.3e1 / 0.2e1 * t136 + t169 + 0.2e1 * t198 + t237 * t134 - t250) * qJD(3) + t238) * t144 + (t40 * t217 + t41 * t216 + t162 * t230 + t160 * t229 + (-t140 * t19 - t143 * t20) * mrSges(5,3) + (mrSges(4,3) + t163) * t89 + (t86 * t164 + t159 * t223 + t161 * t222 + t158 * t218 - t143 * t54 / 0.2e1 + t55 * t217 + (t140 * t46 - t143 * t47) * mrSges(5,3)) * qJD(4) + (t146 + (Ifges(6,5) * t227 + Ifges(6,6) * t228 + (-0.3e1 / 0.2e1 * Ifges(4,4) + t202 / 0.2e1 - t201 / 0.2e1) * t141 + t171 * mrSges(4,1) + (-Ifges(6,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,2) - Ifges(5,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(4,1)) * t144) * qJD(1) + t168) * qJD(3) + (t44 + (m(4) + m(5)) * t89 + (-m(4) * t95 - t127) * qJD(3)) * t134) * t141 + t67 * t61 + t68 * t62 + t72 * t33 + t36 * t83 + t37 * t84 + t102 * t8 + (Ifges(6,5) * t38 + Ifges(6,6) * t39) * t219 + t7 * t227 + t6 * t228 + (Ifges(6,1) * t38 + Ifges(6,4) * t39) * t231 + (Ifges(6,4) * t38 + Ifges(6,2) * t39) * t233; t39 * t51 + t38 * t50 + m(6) * (t10 * t38 - t2 * t99 - t3 * t98 + t9 * t39) - t98 * t17 - t99 * t18 + (-t44 - t8 + (-t140 * t84 + t143 * t83 + t127 - t173) * qJD(3) + m(5) * (t179 * t47 - t181 * t46 - t89) - m(6) * t48) * t144 + (t154 * qJD(4) + m(4) * t88 + m(5) * (-t176 * t46 - t177 * t47 + t157) + t155 + (-t172 + t237 + t242) * qJD(3)) * t141; (t10 * t207 - t118 * t3 - t153 * t2 - t208 * t9) * mrSges(6,3) + t48 * (mrSges(6,1) * t153 + mrSges(6,2) * t118) + (Ifges(6,4) * t118 - Ifges(6,2) * t153) * t235 + (Ifges(6,1) * t118 - Ifges(6,4) * t153) * t236 - t153 * t6 / 0.2e1 + (-Ifges(6,5) * t69 - Ifges(6,6) * t70) * t219 + (-Ifges(6,1) * t69 - Ifges(6,4) * t70) * t231 + (-Ifges(6,4) * t69 - Ifges(6,2) * t70) * t233 + (-t70 / 0.2e1 + t91 / 0.2e1) * t29 + (-t69 / 0.2e1 + t92 / 0.2e1) * t30 + (-Ifges(6,5) * t92 - Ifges(6,6) * t91) * t220 + (-Ifges(6,1) * t92 - Ifges(6,4) * t91) * t232 + (-Ifges(6,4) * t92 - Ifges(6,2) * t91) * t234 + ((t169 + t255 - t198 + t250) * t144 + ((t206 / 0.2e1 + (t256 + Ifges(4,2) / 0.2e1) * t144) * qJD(1) - t146 + t168 + (Ifges(6,5) * t118 - Ifges(6,6) * t153 + t158) * t254) * t141) * qJD(1) + (-mrSges(4,1) - t164) * t89 + t155 * pkin(7) + t157 * mrSges(5,3) - t191 * t95 - m(5) * (t46 * t56 + t47 * t57 + t86 * t95) - pkin(3) * t44 + (t242 * t214 + (-m(5) * t156 + t154) * pkin(7) + t145) * qJD(4) + m(5) * (-pkin(3) * t89 + pkin(7) * t157) + t251 * t51 + t252 * t50 + (t252 * t10 + t135 * t48 + t2 * t77 + t251 * t9 + t3 * t76 - t59 * t73) * m(6) + (-mrSges(6,1) * t207 + mrSges(6,2) * t208) * t59 - t73 * t33 + t76 * t17 + t77 * t18 - t57 * t83 - t56 * t84 - t88 * mrSges(4,2) + t118 * t7 / 0.2e1 - t94 * t127 + t135 * t8 + t140 * t41 / 0.2e1 + t40 * t216 + t159 * t229 + t161 * t230; t66 * t247 - t12 * t50 - t11 * t51 - t238 - m(6) * (t10 * t12 + t11 * t9) + (-Ifges(5,2) * t116 + t112 + t55) * t223 + (t139 * t18 + t142 * t17 + m(6) * (t139 * t2 + t142 * t3) - t242 * t116 + (-t139 * t51 + t142 * t50 + m(6) * (t10 * t142 - t139 * t9)) * qJD(5)) * pkin(4) + (t115 * t46 + t116 * t47) * mrSges(5,3) - t46 * t83 + t47 * t84 - t86 * (mrSges(5,1) * t116 + mrSges(5,2) * t115) + (Ifges(5,5) * t115 - Ifges(5,6) * t116) * t218 + t54 * t221 + (Ifges(5,1) * t115 - t205) * t222 + t132 + t253; t10 * t51 + t29 * t231 - t9 * t50 + t249 + t253;];
tauc = t1(:);
