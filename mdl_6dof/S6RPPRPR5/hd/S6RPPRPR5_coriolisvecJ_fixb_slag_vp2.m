% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:41
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:41:41
% EndTime: 2018-11-23 15:41:44
% DurationCPUTime: 3.25s
% Computational Cost: add. (3087->411), mult. (6657->590), div. (0->0), fcn. (3928->6), ass. (0->190)
t249 = -Ifges(5,5) / 0.2e1;
t142 = qJ(2) - pkin(7);
t147 = cos(qJ(4));
t182 = qJD(4) * t147;
t145 = sin(qJ(4));
t187 = qJD(2) * t145;
t248 = t142 * t182 + t187;
t140 = sin(pkin(9));
t144 = sin(qJ(6));
t141 = cos(pkin(9));
t146 = cos(qJ(6));
t191 = t146 * t141;
t112 = t140 * t144 - t191;
t247 = qJD(6) * t112;
t246 = -Ifges(6,3) / 0.2e1;
t136 = qJD(1) * qJ(2) + qJD(3);
t126 = -qJD(1) * pkin(7) + t136;
t95 = -qJD(4) * pkin(4) - t126 * t147 + qJD(5);
t245 = m(6) * t95;
t244 = qJD(4) * t249;
t194 = t141 * t145;
t178 = pkin(8) * t194;
t153 = (pkin(5) * t147 + t178) * qJD(1);
t162 = pkin(4) * t147 + qJ(5) * t145;
t115 = t162 * qJD(1);
t195 = t140 * t147;
t58 = t141 * t115 - t126 * t195;
t35 = t153 + t58;
t190 = qJD(1) * t145;
t175 = t140 * t190;
t193 = t141 * t147;
t59 = t140 * t115 + t126 * t193;
t46 = pkin(8) * t175 + t59;
t215 = pkin(8) + qJ(5);
t121 = t215 * t140;
t122 = t215 * t141;
t66 = -t121 * t146 - t122 * t144;
t243 = -qJD(5) * t112 + qJD(6) * t66 - t144 * t35 - t146 * t46;
t113 = t140 * t146 + t141 * t144;
t67 = -t121 * t144 + t122 * t146;
t242 = -qJD(5) * t113 - qJD(6) * t67 + t144 * t46 - t146 * t35;
t89 = t113 * t147;
t241 = t112 * qJD(1) - qJD(4) * t89 + t145 * t247;
t239 = t113 * qJD(6);
t91 = t112 * t147;
t98 = t113 * qJD(1);
t240 = -qJD(4) * t91 - t145 * t239 - t98;
t189 = qJD(1) * t147;
t110 = t141 * qJD(4) - t140 * t189;
t181 = t140 * qJD(4);
t111 = t141 * t189 + t181;
t170 = t146 * t110 - t111 * t144;
t92 = qJD(4) * t162 - qJD(5) * t147 + qJD(3);
t71 = t92 * qJD(1);
t180 = qJD(1) * qJD(2);
t86 = t126 * t182 + t145 * t180;
t74 = qJD(4) * qJD(5) + t86;
t29 = -t140 * t74 + t141 * t71;
t17 = qJD(4) * t153 + t29;
t169 = pkin(8) * t145 * t181;
t30 = t140 * t71 + t141 * t74;
t21 = qJD(1) * t169 + t30;
t117 = t145 * t126;
t103 = qJD(4) * qJ(5) + t117;
t143 = pkin(1) + qJ(3);
t118 = pkin(4) * t145 - qJ(5) * t147 + t143;
t84 = qJD(1) * t118 - qJD(2);
t38 = -t103 * t140 + t141 * t84;
t22 = pkin(5) * t190 - pkin(8) * t111 + t38;
t39 = t141 * t103 + t140 * t84;
t28 = pkin(8) * t110 + t39;
t7 = -t144 * t28 + t146 * t22;
t1 = qJD(6) * t7 + t144 * t17 + t146 * t21;
t8 = t144 * t22 + t146 * t28;
t2 = -qJD(6) * t8 - t144 * t21 + t146 * t17;
t183 = qJD(4) * t145;
t151 = t112 * t183;
t26 = qJD(1) * t151 + qJD(6) * t170;
t152 = t113 * t183;
t55 = t110 * t144 + t111 * t146;
t27 = qJD(1) * t152 - qJD(6) * t55;
t238 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t26 + Ifges(7,6) * t27;
t127 = qJD(1) * t143 - qJD(2);
t159 = t39 * t140 + t38 * t141;
t210 = Ifges(6,2) * t140;
t212 = Ifges(6,4) * t141;
t163 = t210 - t212;
t213 = Ifges(6,4) * t140;
t164 = -Ifges(6,1) * t141 + t213;
t165 = -mrSges(6,1) * t140 - mrSges(6,2) * t141;
t222 = -t141 / 0.2e1;
t223 = t140 / 0.2e1;
t237 = t159 * mrSges(6,3) + t95 * t165 + t244 - (t147 * Ifges(5,1) - Ifges(5,4) * t145) * qJD(1) / 0.2e1 + t110 * t163 / 0.2e1 + t111 * t164 / 0.2e1 - t127 * mrSges(5,2) + (t111 * Ifges(6,4) + t110 * Ifges(6,2) + Ifges(6,6) * t190) * t223 + (t111 * Ifges(6,1) + t110 * Ifges(6,4) + Ifges(6,5) * t190) * t222;
t236 = (t145 ^ 2 + t147 ^ 2) * t126;
t235 = t26 / 0.2e1;
t234 = t27 / 0.2e1;
t233 = -t170 / 0.2e1;
t232 = t170 / 0.2e1;
t231 = -t55 / 0.2e1;
t230 = t55 / 0.2e1;
t229 = -t89 / 0.2e1;
t228 = -t91 / 0.2e1;
t227 = -t112 / 0.2e1;
t226 = t113 / 0.2e1;
t131 = qJD(6) + t190;
t225 = -t131 / 0.2e1;
t224 = t131 / 0.2e1;
t221 = t141 / 0.2e1;
t179 = qJD(1) * qJD(4);
t82 = t165 * t145 * t179;
t9 = -t27 * mrSges(7,1) + t26 * mrSges(7,2);
t220 = -t82 - t9;
t219 = Ifges(7,4) * t55;
t218 = pkin(5) * t140;
t81 = -t144 * t175 + t190 * t191;
t214 = t247 - t81;
t211 = Ifges(6,5) * t141;
t209 = Ifges(6,6) * t140;
t135 = t147 * t180;
t85 = t126 * t183 - t135;
t205 = t85 * t147;
t80 = t145 * t98;
t204 = t239 + t80;
t200 = qJD(4) * mrSges(5,1);
t203 = -mrSges(6,1) * t110 + mrSges(6,2) * t111 + mrSges(5,3) * t189 - t200;
t201 = Ifges(5,6) * qJD(4);
t199 = qJD(4) * mrSges(5,2);
t196 = t140 * t145;
t192 = t142 * t145;
t70 = t140 * t118 + t141 * t192;
t186 = qJD(2) * t147;
t185 = qJD(3) * t127;
t184 = qJD(3) * t143;
t18 = -mrSges(7,1) * t170 + mrSges(7,2) * t55;
t177 = t18 + t203;
t176 = t142 * t205;
t45 = t140 * t92 + t248 * t141;
t173 = t147 * t179;
t172 = -t140 * t142 + pkin(5);
t171 = -t142 + t218;
t168 = m(3) * qJ(2) + mrSges(4,2) + mrSges(3,3);
t166 = mrSges(5,1) * t145 + mrSges(5,2) * t147;
t161 = -t30 * t140 - t29 * t141;
t160 = -t140 * t29 + t141 * t30;
t158 = -t140 * t38 + t141 * t39;
t107 = t141 * t118;
t50 = -pkin(8) * t193 + t145 * t172 + t107;
t56 = -pkin(8) * t195 + t70;
t15 = -t144 * t56 + t146 * t50;
t16 = t144 * t50 + t146 * t56;
t123 = -mrSges(5,3) * t190 - t199;
t77 = -mrSges(6,2) * t190 + mrSges(6,3) * t110;
t78 = mrSges(6,1) * t190 - mrSges(6,3) * t111;
t157 = -t140 * t78 + t141 * t77 + t123;
t156 = t211 / 0.2e1 - t209 / 0.2e1;
t150 = t39 * mrSges(6,2) + t8 * mrSges(7,2) + t201 / 0.2e1 + (Ifges(5,4) * t147 - t145 * Ifges(5,2)) * qJD(1) / 0.2e1 - t131 * Ifges(7,3) - t55 * Ifges(7,5) - t170 * Ifges(7,6) + t190 * t246 - t111 * Ifges(6,5) - t110 * Ifges(6,6) - t127 * mrSges(5,1) - t38 * mrSges(6,1) - t7 * mrSges(7,1);
t148 = qJD(1) ^ 2;
t133 = -pkin(5) * t141 - pkin(4);
t130 = Ifges(7,3) * t173;
t116 = qJD(1) * t166;
t109 = t171 * t147;
t94 = (mrSges(6,1) * t147 + mrSges(6,3) * t194) * t179;
t93 = (-mrSges(6,2) * t147 + mrSges(6,3) * t196) * t179;
t90 = t112 * t145;
t88 = t113 * t145;
t83 = -pkin(5) * t175 + t117;
t79 = -t171 * t183 - t186;
t76 = t141 * t92;
t69 = -t140 * t192 + t107;
t65 = (Ifges(6,5) * t147 + t145 * t164) * t179;
t64 = (Ifges(6,6) * t147 + t145 * t163) * t179;
t61 = -t135 + (-qJD(1) * t218 + t126) * t183;
t57 = -pkin(5) * t110 + t95;
t51 = Ifges(7,4) * t170;
t44 = -t248 * t140 + t76;
t43 = t147 * t247 + t152;
t41 = -t147 * t239 + t151;
t37 = mrSges(7,1) * t131 - mrSges(7,3) * t55;
t36 = -mrSges(7,2) * t131 + mrSges(7,3) * t170;
t32 = t169 + t45;
t31 = -t140 * t187 + t76 + (t147 * t172 + t178) * qJD(4);
t20 = -mrSges(7,2) * t173 + mrSges(7,3) * t27;
t19 = mrSges(7,1) * t173 - mrSges(7,3) * t26;
t14 = Ifges(7,1) * t55 + Ifges(7,5) * t131 + t51;
t13 = Ifges(7,2) * t170 + Ifges(7,6) * t131 + t219;
t6 = t26 * Ifges(7,1) + t27 * Ifges(7,4) + Ifges(7,5) * t173;
t5 = t26 * Ifges(7,4) + t27 * Ifges(7,2) + Ifges(7,6) * t173;
t4 = -qJD(6) * t16 - t144 * t32 + t146 * t31;
t3 = qJD(6) * t15 + t144 * t31 + t146 * t32;
t10 = [m(7) * (t1 * t16 + t109 * t61 + t15 * t2 + t3 * t8 + t4 * t7 + t57 * t79) + ((t142 * t123 - t150 - t201 / 0.2e1) * t147 + (t244 + (t203 + t245) * t142 + t237) * t145) * qJD(4) + (t29 * mrSges(6,1) - t30 * mrSges(6,2) + t130 / 0.2e1 + qJD(2) * t123 - t86 * mrSges(5,3) + t238) * t145 + m(5) * (qJD(2) * t236 + t192 * t86 - t176 + t185) + qJD(3) * t116 + t109 * t9 + t70 * t93 + t69 * t94 + t45 * t77 + t44 * t78 + t79 * t18 + t57 * (-mrSges(7,1) * t43 + mrSges(7,2) * t41) + t41 * t14 / 0.2e1 + t43 * t13 / 0.2e1 + t3 * t36 + t4 * t37 + t15 * t19 + t16 * t20 + (((2 * mrSges(4,3)) + t166) * qJD(3) + 0.2e1 * t168 * qJD(2) + m(5) * t184 + m(4) * (qJ(2) * qJD(2) + t184) + ((t143 * mrSges(5,1) + Ifges(7,5) * t228 + Ifges(7,6) * t229 + (-0.3e1 / 0.2e1 * Ifges(5,4) + t156) * t147) * t147 + (-t143 * mrSges(5,2) + (-0.3e1 / 0.2e1 * t211 + 0.3e1 / 0.2e1 * t209 + 0.3e1 / 0.2e1 * Ifges(5,4)) * t145 + (0.3e1 / 0.2e1 * Ifges(6,3) - 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,2) + Ifges(7,3) / 0.2e1 - Ifges(6,1) * t141 ^ 2 / 0.2e1 + (t212 - t210 / 0.2e1) * t140) * t147) * t145) * qJD(4)) * qJD(1) + t61 * (mrSges(7,1) * t89 - mrSges(7,2) * t91) + (-t1 * t89 + t2 * t91 - t7 * t41 + t8 * t43) * mrSges(7,3) + (-Ifges(7,4) * t91 - Ifges(7,2) * t89) * t234 + (-Ifges(7,1) * t91 - Ifges(7,4) * t89) * t235 + (t65 * t221 - t142 * t82 - t140 * t64 / 0.2e1 + (mrSges(5,3) - t165) * t85 - t203 * qJD(2) + t161 * mrSges(6,3)) * t147 + m(4) * (qJD(2) * t136 + t185) + m(6) * (-t186 * t95 + t29 * t69 + t30 * t70 + t38 * t44 + t39 * t45 - t176) + (Ifges(7,5) * t41 + Ifges(7,6) * t43) * t224 + t6 * t228 + t5 * t229 + (Ifges(7,1) * t41 + Ifges(7,4) * t43) * t230 + (Ifges(7,4) * t41 + Ifges(7,2) * t43) * t232; t112 * t19 - t113 * t20 - t140 * t93 - t141 * t94 + t204 * t37 + t214 * t36 - t168 * t148 - m(7) * (-t7 * t80 + t8 * t81) + m(6) * t161 + m(7) * (-t1 * t113 + t112 * t2 + t239 * t7 + t247 * t8) + ((-t136 - qJD(3)) * m(4) + (-t157 + t199) * t145 - m(6) * (t194 * t39 - t196 * t38) + (m(7) * t57 + t177 - t200 + t245) * t147 + (-qJD(3) - t236) * m(5)) * qJD(1); -t148 * mrSges(4,3) - t88 * t19 - t90 * t20 + t241 * t37 + t240 * t36 + t220 * t147 + (-t140 * t94 + t141 * t93) * t145 + (t177 * t145 + t157 * t147) * qJD(4) + m(5) * (t86 * t145 - t205) + m(6) * (-t205 + t160 * t145 + (t145 * t95 + t147 * t158) * qJD(4)) + (-t116 - m(6) * t159 - t141 * t78 - t140 * t77 - m(5) * t127 + (-t127 + qJD(2)) * m(4)) * qJD(1) + (-t1 * t90 - t147 * t61 + t57 * t183 - t2 * t88 + t240 * t8 + t241 * t7) * m(7); (-pkin(4) * t85 + qJ(5) * t160 + qJD(5) * t158 - t117 * t95 - t38 * t58 - t39 * t59) * m(6) + t242 * t37 + (t1 * t67 + t133 * t61 + t2 * t66 + t242 * t7 + t243 * t8 - t57 * t83) * m(7) + t243 * t36 + t133 * t9 + t61 * (mrSges(7,1) * t112 + mrSges(7,2) * t113) - pkin(4) * t82 - t83 * t18 - t85 * mrSges(5,1) - t86 * mrSges(5,2) - t59 * t77 - t58 * t78 + t66 * t19 + t67 * t20 + (-t239 / 0.2e1 - t80 / 0.2e1) * t13 + (-t85 * mrSges(6,1) + t64 / 0.2e1 + qJ(5) * t93 + qJD(5) * t77 + t30 * mrSges(6,3)) * t141 + (t85 * mrSges(6,2) + t65 / 0.2e1 - qJ(5) * t94 - qJD(5) * t78 - t29 * mrSges(6,3)) * t140 + (-t1 * t112 - t113 * t2 - t204 * t8 + t214 * t7) * mrSges(7,3) + (mrSges(7,1) * t204 - mrSges(7,2) * t214) * t57 + (-t147 * t123 - t145 * t203) * t126 + ((Ifges(5,4) * t189 / 0.2e1 + t150) * t147 + ((-Ifges(5,4) / 0.2e1 + t156) * t190 + (t246 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t189 - t237) * t145 + ((Ifges(6,5) * t223 + Ifges(6,6) * t221 + Ifges(7,5) * t226 + Ifges(7,6) * t227 - Ifges(5,6) / 0.2e1) * t147 + (t249 + (Ifges(6,2) * t141 + t213) * t223 + (Ifges(6,1) * t140 + t212) * t222) * t145) * qJD(4)) * qJD(1) + (-t247 / 0.2e1 + t81 / 0.2e1) * t14 + (-Ifges(7,5) * t247 - Ifges(7,6) * t239) * t224 + (-Ifges(7,1) * t247 - Ifges(7,4) * t239) * t230 + (-Ifges(7,4) * t247 - Ifges(7,2) * t239) * t232 + (-Ifges(7,5) * t81 + Ifges(7,6) * t80) * t225 + t6 * t226 + t5 * t227 + (-Ifges(7,1) * t81 + Ifges(7,4) * t80) * t231 + (-Ifges(7,4) * t81 + Ifges(7,2) * t80) * t233 + (Ifges(7,4) * t113 - Ifges(7,2) * t112) * t234 + (Ifges(7,1) * t113 - Ifges(7,4) * t112) * t235; -t110 * t77 + t111 * t78 - t170 * t36 + t55 * t37 - t220 + (-t170 * t8 + t55 * t7 + t61) * m(7) + (-t39 * t110 + t38 * t111 + t85) * m(6); t130 - t57 * (mrSges(7,1) * t55 + mrSges(7,2) * t170) + (Ifges(7,1) * t170 - t219) * t231 + t13 * t230 + (Ifges(7,5) * t170 - Ifges(7,6) * t55) * t225 - t7 * t36 + t8 * t37 + (t170 * t7 + t55 * t8) * mrSges(7,3) + (-Ifges(7,2) * t55 + t14 + t51) * t233 + t238;];
tauc  = t10(:);
