% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:43:00
% EndTime: 2019-03-08 19:43:10
% DurationCPUTime: 5.94s
% Computational Cost: add. (3888->430), mult. (10181->577), div. (0->0), fcn. (7553->10), ass. (0->203)
t132 = sin(pkin(11));
t134 = cos(pkin(11));
t137 = sin(qJ(4));
t228 = cos(qJ(4));
t112 = t132 * t228 + t137 * t134;
t133 = sin(pkin(6));
t140 = cos(qJ(2));
t195 = t133 * t140;
t146 = t112 * t195;
t214 = pkin(8) + qJ(3);
t118 = t214 * t132;
t119 = t214 * t134;
t78 = -t137 * t118 + t119 * t228;
t248 = -qJD(1) * t146 + qJD(3) * t112 + qJD(4) * t78;
t179 = qJD(4) * t228;
t197 = t132 * t137;
t180 = qJD(4) * t197;
t105 = -t134 * t179 + t180;
t256 = -t105 * pkin(5) + t248;
t106 = t112 * qJD(4);
t155 = qJ(5) * t105 - qJD(5) * t112;
t138 = sin(qJ(2));
t196 = t133 * t138;
t183 = qJD(1) * t196;
t234 = pkin(4) + pkin(9);
t255 = -t106 * t234 - t155 + t183;
t216 = mrSges(6,2) - mrSges(5,1);
t254 = mrSges(5,3) + mrSges(6,1);
t136 = sin(qJ(6));
t139 = cos(qJ(6));
t185 = t228 * t134;
t111 = -t185 + t197;
t128 = -pkin(3) * t134 - pkin(2);
t152 = -qJ(5) * t112 + t128;
t50 = t111 * t234 + t152;
t77 = t228 * t118 + t119 * t137;
t57 = pkin(5) * t112 + t77;
t19 = -t136 * t50 + t139 * t57;
t253 = qJD(6) * t19 + t136 * t256 - t255 * t139;
t20 = t136 * t57 + t139 * t50;
t252 = -qJD(6) * t20 + t255 * t136 + t139 * t256;
t251 = Ifges(6,4) - Ifges(5,5);
t250 = Ifges(6,5) - Ifges(5,6);
t55 = t137 * (qJD(3) * t132 + qJD(4) * t119) - qJD(3) * t185 + t118 * t179;
t156 = t185 * t195;
t190 = qJD(1) * t140;
t182 = t133 * t190;
t80 = qJD(1) * t156 - t182 * t197;
t249 = t55 + t80;
t192 = t132 ^ 2 + t134 ^ 2;
t247 = mrSges(4,3) * t192;
t175 = qJD(2) * t185;
t189 = qJD(2) * t132;
t103 = t137 * t189 - t175;
t82 = -qJD(4) * t136 + t103 * t139;
t83 = qJD(4) * t139 + t103 * t136;
t44 = -mrSges(7,1) * t82 + mrSges(7,2) * t83;
t88 = mrSges(6,1) * t103 - qJD(4) * mrSges(6,3);
t210 = -t88 + t44;
t86 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t103;
t187 = -t86 - t210;
t104 = t112 * qJD(2);
t202 = pkin(8) * qJD(2);
t117 = qJD(2) * qJ(3) + t183;
t135 = cos(pkin(6));
t191 = qJD(1) * t135;
t85 = t134 * t117 + t132 * t191;
t76 = t134 * t202 + t85;
t201 = t137 * t76;
t123 = t134 * t191;
t75 = t123 + (-t117 - t202) * t132;
t38 = -t228 * t75 + t201;
t151 = pkin(5) * t104 + t38;
t246 = qJD(5) + t151;
t71 = -mrSges(6,2) * t103 - mrSges(6,3) * t104;
t245 = -mrSges(5,1) * t103 - mrSges(5,2) * t104 - t71;
t211 = -qJD(4) * t216 - t254 * t104;
t244 = -qJD(5) - t38;
t110 = (qJD(3) + t182) * qJD(2);
t39 = t137 * t75 + t228 * t76;
t17 = qJD(4) * t39 + t110 * t112;
t96 = qJD(2) * t180 - qJD(4) * t175;
t12 = -t96 * pkin(5) + t17;
t188 = qJD(2) * t138;
t181 = t133 * t188;
t120 = qJD(1) * t181;
t148 = qJ(5) * t96 - qJD(5) * t104 + t120;
t97 = qJD(2) * t106;
t22 = t234 * t97 + t148;
t21 = -qJD(4) * t234 + t246;
t161 = qJD(3) - t182;
t98 = qJD(2) * t128 + t161;
t143 = -qJ(5) * t104 + t98;
t31 = t103 * t234 + t143;
t5 = -t136 * t31 + t139 * t21;
t1 = qJD(6) * t5 + t12 * t136 + t139 * t22;
t6 = t136 * t21 + t139 * t31;
t2 = -qJD(6) * t6 + t12 * t139 - t136 * t22;
t243 = t1 * t136 + t139 * t2;
t101 = qJD(6) + t104;
t51 = qJD(6) * t82 + t136 * t97;
t52 = -qJD(6) * t83 + t139 * t97;
t242 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t51 + Ifges(7,6) * t52;
t213 = t110 * t185 + t75 * t179;
t15 = qJD(4) * (-qJD(5) + t201) + t110 * t197 - t213;
t209 = Ifges(7,4) * t136;
t163 = Ifges(7,2) * t139 + t209;
t208 = Ifges(7,4) * t139;
t165 = Ifges(7,1) * t136 + t208;
t168 = mrSges(7,1) * t139 - mrSges(7,2) * t136;
t170 = t136 * t5 - t139 * t6;
t204 = Ifges(7,6) * t139;
t207 = Ifges(7,5) * t136;
t223 = t103 * pkin(5);
t33 = -qJD(4) * qJ(5) - t39;
t23 = -t33 - t223;
t230 = -t136 / 0.2e1;
t233 = -t101 / 0.2e1;
t236 = -t83 / 0.2e1;
t238 = -t82 / 0.2e1;
t227 = Ifges(7,4) * t83;
t27 = Ifges(7,2) * t82 + Ifges(7,6) * t101 + t227;
t81 = Ifges(7,4) * t82;
t28 = Ifges(7,1) * t83 + Ifges(7,5) * t101 + t81;
t241 = (t204 + t207) * t233 + t163 * t238 + t165 * t236 + t23 * t168 + t170 * mrSges(7,3) + t28 * t230 - t139 * t27 / 0.2e1;
t240 = t51 / 0.2e1;
t239 = t52 / 0.2e1;
t237 = t82 / 0.2e1;
t235 = t83 / 0.2e1;
t232 = t101 / 0.2e1;
t229 = t139 / 0.2e1;
t226 = Ifges(7,5) * t83;
t225 = Ifges(7,6) * t82;
t102 = t132 * t135 + t134 * t196;
t149 = t132 * t196 - t134 * t135;
t147 = t228 * t149;
t63 = t102 * t137 + t147;
t221 = t17 * t63;
t220 = t17 * t77;
t64 = t102 * t228 - t137 * t149;
t219 = t64 * t97;
t218 = t96 * mrSges(6,1);
t215 = -Ifges(5,4) - Ifges(6,6);
t212 = t86 - t88;
t206 = Ifges(7,5) * t139;
t205 = Ifges(7,6) * t136;
t203 = Ifges(7,3) * t101;
t200 = qJ(5) * t103;
t199 = t106 * t136;
t198 = t106 * t139;
t184 = t133 ^ 2 * t190;
t178 = t192 * t110;
t173 = t1 * t139 - t2 * t136;
t171 = t136 * t6 + t139 * t5;
t169 = -t77 * t96 - t78 * t97;
t167 = mrSges(7,1) * t136 + mrSges(7,2) * t139;
t166 = Ifges(7,1) * t139 - t209;
t164 = -Ifges(7,2) * t136 + t208;
t160 = -t132 * (-t117 * t132 + t123) + t134 * t85;
t29 = -mrSges(7,1) * t96 - mrSges(7,3) * t51;
t30 = mrSges(7,2) * t96 + mrSges(7,3) * t52;
t159 = t136 * t30 + t139 * t29;
t53 = -mrSges(7,2) * t101 + mrSges(7,3) * t82;
t54 = mrSges(7,1) * t101 - mrSges(7,3) * t83;
t158 = -t136 * t54 + t139 * t53;
t157 = -t136 * t53 - t139 * t54;
t150 = -t136 * t63 + t139 * t195;
t45 = t136 * t195 + t139 * t63;
t145 = ((-qJD(2) * pkin(2) + t161) * t138 + t140 * t160) * t133;
t142 = qJD(2) ^ 2;
t100 = Ifges(5,4) * t103;
t99 = Ifges(6,6) * t103;
t93 = Ifges(7,3) * t96;
t92 = t96 * mrSges(5,2);
t91 = t96 * mrSges(6,3);
t72 = pkin(4) * t111 + t152;
t69 = pkin(4) * t104 + t200;
t68 = Ifges(5,1) * t104 + Ifges(5,5) * qJD(4) - t100;
t67 = t104 * Ifges(5,4) - Ifges(5,2) * t103 + Ifges(5,6) * qJD(4);
t66 = Ifges(6,4) * qJD(4) - Ifges(6,2) * t104 + t99;
t65 = Ifges(6,5) * qJD(4) - t104 * Ifges(6,6) + Ifges(6,3) * t103;
t60 = t97 * mrSges(5,1) - t92;
t59 = -t97 * mrSges(6,2) + t91;
t58 = -t111 * pkin(5) + t78;
t47 = pkin(4) * t106 + t155;
t43 = t104 * t234 + t200;
t42 = pkin(4) * t103 + t143;
t40 = -pkin(5) * t106 - t55;
t37 = pkin(4) * t97 + t148;
t35 = qJD(2) * t146 + qJD(4) * t64;
t32 = -qJD(4) * pkin(4) - t244;
t26 = t203 + t225 + t226;
t25 = t39 - t223;
t18 = -mrSges(7,1) * t52 + mrSges(7,2) * t51;
t16 = (-qJD(4) * t76 - t110 * t132) * t137 + t213;
t14 = Ifges(7,1) * t51 + Ifges(7,4) * t52 - t96 * Ifges(7,5);
t13 = Ifges(7,4) * t51 + Ifges(7,2) * t52 - Ifges(7,6) * t96;
t11 = -pkin(5) * t97 - t15;
t10 = qJD(6) * t45 + t136 * t35 + t139 * t181;
t9 = qJD(6) * t150 - t136 * t181 + t139 * t35;
t8 = t136 * t25 + t139 * t43;
t7 = -t136 * t43 + t139 * t25;
t3 = [m(7) * (-t1 * t150 + t10 * t6 + t11 * t64 + t2 * t45 + t5 * t9) + m(4) * ((t134 * t102 + t132 * t149) * t110 + (-t138 * t184 + t145) * qJD(2)) - t63 * t218 - mrSges(6,1) * t219 + m(5) * (t16 * t64 + t221 + (t133 * t98 - t184) * t188) + m(6) * (-t15 * t64 + t221 + (-t140 * t37 + t188 * t42) * t133) + t45 * t29 - t150 * t30 + t10 * t53 + t9 * t54 + t64 * t18 + (m(5) * t38 + m(6) * t32 - t211) * t35 + (-t60 - t59) * t195 + (-t63 * t96 - t219) * mrSges(5,3) + (-m(5) * t39 + m(6) * t33 - m(7) * t23 + t187) * (-qJD(2) * t156 + qJD(4) * t147 + (qJD(4) * t102 + t189 * t195) * t137) + (qJD(2) * (-mrSges(4,1) * t134 + mrSges(4,2) * t132) - t245) * t181 + (-t138 * mrSges(3,1) + (-mrSges(3,2) + t247) * t140) * t133 * t142; (qJD(2) * t161 * t192 + t178) * mrSges(4,3) + (-t105 * t38 - t106 * t39 + t169) * mrSges(5,3) + (-t105 * t32 + t106 * t33 + t169) * mrSges(6,1) + (-pkin(2) * t120 + qJ(3) * t178 - qJD(1) * t145 + qJD(3) * t160) * m(4) + t245 * t183 + t27 * t198 / 0.2e1 + t6 * (mrSges(7,2) * t105 + mrSges(7,3) * t198) + t187 * t80 + t252 * t54 + t253 * t53 + (t1 * t20 + t11 * t58 + t19 * t2 + t253 * t6 + t252 * t5 + (t40 - t80) * t23) * m(7) + (mrSges(5,2) * t120 - t93 / 0.2e1 - t37 * mrSges(6,3) + t215 * t97 + t254 * t17 + (-Ifges(5,1) - Ifges(6,2) - Ifges(7,3) / 0.2e1) * t96 + t242) * t112 + (Ifges(7,5) * t199 + Ifges(7,6) * t198 - Ifges(7,3) * t105) * t232 + (Ifges(7,1) * t199 + Ifges(7,4) * t198 - Ifges(7,5) * t105) * t235 + (Ifges(7,4) * t199 + Ifges(7,2) * t198 - Ifges(7,6) * t105) * t237 + t28 * t199 / 0.2e1 + t5 * (-mrSges(7,1) * t105 - mrSges(7,3) * t199) + t23 * (-mrSges(7,1) * t198 + mrSges(7,2) * t199) + t105 * t66 / 0.2e1 - t104 * (Ifges(6,2) * t105 + Ifges(6,6) * t106) / 0.2e1 + t103 * (Ifges(6,6) * t105 + Ifges(6,3) * t106) / 0.2e1 + t42 * (-mrSges(6,2) * t106 + mrSges(6,3) * t105) + t106 * t65 / 0.2e1 + t104 * (-Ifges(5,1) * t105 - Ifges(5,4) * t106) / 0.2e1 + t98 * (mrSges(5,1) * t106 - mrSges(5,2) * t105) - t103 * (-Ifges(5,4) * t105 - Ifges(5,2) * t106) / 0.2e1 - t106 * t67 / 0.2e1 + t19 * t29 + t20 * t30 - t212 * t55 + t128 * t60 + (-t11 * t168 + t165 * t240 + t163 * t239 + t136 * t14 / 0.2e1 + t13 * t229 - t37 * mrSges(6,2) + t15 * mrSges(6,1) - t16 * mrSges(5,3) + mrSges(5,1) * t120 + (Ifges(6,3) + Ifges(5,2)) * t97 + (-t207 / 0.2e1 - t204 / 0.2e1 - t215) * t96 + t173 * mrSges(7,3) + ((-t205 + t206) * t232 + t164 * t237 + t166 * t235 + t23 * t167 + t28 * t229 + t27 * t230 - t171 * mrSges(7,3)) * qJD(6)) * t111 + t40 * t44 + t58 * t18 + t47 * t71 + t72 * t59 - t248 * t211 + (-t15 * t78 + t37 * t72 + t220 + (-t183 + t47) * t42 + t249 * t33 + t248 * t32) * m(6) + (t120 * t128 + t16 * t78 - t183 * t98 + t248 * t38 - t249 * t39 + t220) * m(5) + (t105 * t251 + t106 * t250) * qJD(4) / 0.2e1 - (t68 + t26) * t105 / 0.2e1; -t136 * t29 + t139 * t30 + t91 - t92 - t216 * t97 + t157 * qJD(6) + (m(4) + m(5)) * t120 - t142 * t247 - t187 * t103 + (t157 + t211) * t104 - m(5) * (-t103 * t39 + t104 * t38) - m(4) * t160 * qJD(2) + (-t101 * t171 + t103 * t23 + t173) * m(7) + (-t33 * t103 - t32 * t104 + t37) * m(6); -(t159 + m(7) * t243 + (-m(7) * t170 + t158) * qJD(6)) * t234 + (-pkin(4) * t17 - qJ(5) * t15 + t244 * t33 - t32 * t39 - t42 * t69) * m(6) + (t42 * mrSges(6,2) - t65 / 0.2e1 - t98 * mrSges(5,1) + t67 / 0.2e1 - t33 * mrSges(6,1) + t39 * mrSges(5,3) + (Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t104 + (-Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1) * qJD(4) + (-Ifges(6,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t103 + t241) * t104 + t241 * qJD(6) + t11 * t167 - t243 * mrSges(7,3) + t14 * t229 + t13 * t230 + t164 * t239 + t166 * t240 + (qJ(5) * t11 + t246 * t23 - t5 * t7 - t6 * t8) * m(7) + t210 * qJD(5) + t211 * t39 + t212 * t38 + t216 * t17 + (-t99 / 0.2e1 - t100 / 0.2e1 + t5 * mrSges(7,1) - t6 * mrSges(7,2) - t42 * mrSges(6,3) + t98 * mrSges(5,2) - t66 / 0.2e1 + t26 / 0.2e1 + t68 / 0.2e1 + t32 * mrSges(6,1) + t38 * mrSges(5,3) + t203 / 0.2e1 + t225 / 0.2e1 + t226 / 0.2e1 + (Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1) * qJD(4)) * t103 - t8 * t53 - t7 * t54 - t69 * t71 + (-mrSges(6,1) * qJ(5) + t250) * t97 + (pkin(4) * mrSges(6,1) - t206 / 0.2e1 + t205 / 0.2e1 + t251) * t96 - t16 * mrSges(5,2) + qJ(5) * t18 + t151 * t44 - t15 * mrSges(6,3); -t218 + t158 * qJD(6) - t210 * qJD(4) + (t158 + t71) * t104 + t159 + (-qJD(4) * t23 - t101 * t170 + t243) * m(7) + (qJD(4) * t33 + t104 * t42 + t17) * m(6); -t93 - t23 * (mrSges(7,1) * t83 + mrSges(7,2) * t82) + (Ifges(7,1) * t82 - t227) * t236 + t27 * t235 + (Ifges(7,5) * t82 - Ifges(7,6) * t83) * t233 - t5 * t53 + t6 * t54 + (t5 * t82 + t6 * t83) * mrSges(7,3) + (-Ifges(7,2) * t83 + t28 + t81) * t238 + t242;];
tauc  = t3(:);
