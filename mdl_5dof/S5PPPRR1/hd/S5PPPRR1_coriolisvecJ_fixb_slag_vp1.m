% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPPRR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPPRR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:50
% EndTime: 2019-12-05 14:58:04
% DurationCPUTime: 6.42s
% Computational Cost: add. (13098->546), mult. (19819->860), div. (0->0), fcn. (21606->8), ass. (0->218)
t173 = pkin(9) + qJ(4);
t171 = sin(t173);
t175 = sin(pkin(7));
t176 = cos(pkin(8));
t225 = t175 * t176;
t172 = cos(t173);
t177 = cos(pkin(7));
t227 = t172 * t177;
t150 = t171 * t225 + t227;
t174 = sin(pkin(8));
t216 = qJD(4) * t174;
t206 = t175 * t216;
t134 = qJD(5) * t150 + t206;
t251 = t134 / 0.2e1;
t224 = t177 * t171;
t152 = -t175 * t172 + t176 * t224;
t214 = qJD(4) * t177;
t205 = t174 * t214;
t135 = qJD(5) * t152 + t205;
t249 = t135 / 0.2e1;
t210 = qJD(5) * t174;
t215 = qJD(4) * t176;
t163 = t171 * t210 - t215;
t247 = t163 / 0.2e1;
t255 = qJD(5) / 0.2e1;
t154 = (-Icges(5,5) * t171 - Icges(5,6) * t172) * t174;
t143 = qJD(4) * t154;
t178 = sin(qJ(5));
t179 = cos(qJ(5));
t228 = t172 * t174;
t161 = -t176 * t179 - t178 * t228;
t162 = -t176 * t178 + t179 * t228;
t229 = t171 * t174;
t102 = Icges(6,5) * t162 + Icges(6,6) * t161 + Icges(6,3) * t229;
t230 = Icges(6,4) * t162;
t103 = Icges(6,2) * t161 + Icges(6,6) * t229 + t230;
t157 = Icges(6,4) * t161;
t104 = Icges(6,1) * t162 + Icges(6,5) * t229 + t157;
t209 = t172 * t225;
t151 = t209 - t224;
t226 = t175 * t174;
t128 = -t151 * t178 + t179 * t226;
t129 = t151 * t179 + t178 * t226;
t140 = qJD(4) * t209 - t171 * t214;
t208 = t171 * t216;
t132 = -qJD(5) * t162 + t178 * t208;
t133 = qJD(5) * t161 - t179 * t208;
t207 = t172 * t216;
t62 = Icges(6,5) * t133 + Icges(6,6) * t132 + Icges(6,3) * t207;
t63 = Icges(6,4) * t133 + Icges(6,2) * t132 + Icges(6,6) * t207;
t64 = Icges(6,1) * t133 + Icges(6,4) * t132 + Icges(6,5) * t207;
t139 = t150 * qJD(4);
t82 = -qJD(5) * t129 + t139 * t178;
t83 = qJD(5) * t128 - t139 * t179;
t14 = t102 * t140 + t103 * t82 + t104 * t83 + t128 * t63 + t129 * t64 + t150 * t62;
t153 = t171 * t175 + t176 * t227;
t142 = t153 * qJD(4);
t52 = Icges(6,5) * t129 + Icges(6,6) * t128 + Icges(6,3) * t150;
t232 = Icges(6,4) * t129;
t54 = Icges(6,2) * t128 + Icges(6,6) * t150 + t232;
t122 = Icges(6,4) * t128;
t56 = Icges(6,1) * t129 + Icges(6,5) * t150 + t122;
t19 = t128 * t54 + t129 * t56 + t150 * t52;
t223 = t177 * t174;
t130 = -t153 * t178 + t179 * t223;
t131 = t153 * t179 + t178 * t223;
t53 = Icges(6,5) * t131 + Icges(6,6) * t130 + Icges(6,3) * t152;
t231 = Icges(6,4) * t131;
t55 = Icges(6,2) * t130 + Icges(6,6) * t152 + t231;
t123 = Icges(6,4) * t130;
t57 = Icges(6,1) * t131 + Icges(6,5) * t152 + t123;
t20 = t128 * t55 + t129 * t57 + t150 * t53;
t197 = t140 * t19 + t142 * t20;
t34 = t102 * t150 + t103 * t128 + t104 * t129;
t40 = Icges(6,5) * t83 + Icges(6,6) * t82 + Icges(6,3) * t140;
t42 = Icges(6,4) * t83 + Icges(6,2) * t82 + Icges(6,6) * t140;
t44 = Icges(6,1) * t83 + Icges(6,4) * t82 + Icges(6,5) * t140;
t6 = t128 * t42 + t129 * t44 + t140 * t52 + t150 * t40 + t54 * t82 + t56 * t83;
t141 = t152 * qJD(4);
t84 = -qJD(5) * t131 + t141 * t178;
t85 = qJD(5) * t130 - t141 * t179;
t41 = Icges(6,5) * t85 + Icges(6,6) * t84 + Icges(6,3) * t142;
t43 = Icges(6,4) * t85 + Icges(6,2) * t84 + Icges(6,6) * t142;
t45 = Icges(6,1) * t85 + Icges(6,4) * t84 + Icges(6,5) * t142;
t7 = t128 * t43 + t129 * t45 + t140 * t53 + t150 * t41 + t55 * t82 + t57 * t83;
t254 = t6 * t251 + t7 * t249 + t14 * t247 + (t207 * t34 + t197) * t255;
t15 = t102 * t142 + t103 * t84 + t104 * t85 + t130 * t63 + t131 * t64 + t152 * t62;
t21 = t130 * t54 + t131 * t56 + t152 * t52;
t22 = t130 * t55 + t131 * t57 + t152 * t53;
t196 = t140 * t21 + t142 * t22;
t35 = t102 * t152 + t103 * t130 + t104 * t131;
t8 = t130 * t42 + t131 * t44 + t142 * t52 + t152 * t40 + t54 * t84 + t56 * t85;
t9 = t130 * t43 + t131 * t45 + t142 * t53 + t152 * t41 + t55 * t84 + t57 * t85;
t253 = t8 * t251 + t9 * t249 + t15 * t247 + (t207 * t35 + t196) * t255;
t252 = -t134 / 0.2e1;
t250 = -t135 / 0.2e1;
t248 = -t163 / 0.2e1;
t100 = -pkin(4) * t139 + pkin(6) * t140;
t46 = rSges(6,1) * t83 + rSges(6,2) * t82 + rSges(6,3) * t140;
t245 = t100 + t46;
t101 = -pkin(4) * t141 + pkin(6) * t142;
t47 = rSges(6,1) * t85 + rSges(6,2) * t84 + rSges(6,3) * t142;
t244 = -t101 - t47;
t147 = Icges(5,4) * t150;
t88 = Icges(5,1) * t151 + Icges(5,5) * t226 - t147;
t243 = Icges(5,2) * t151 + t147 - t88;
t148 = Icges(5,4) * t152;
t89 = Icges(5,1) * t153 + Icges(5,5) * t223 - t148;
t242 = Icges(5,2) * t153 + t148 - t89;
t236 = Icges(5,4) * t151;
t86 = -Icges(5,2) * t150 + Icges(5,6) * t226 + t236;
t241 = -Icges(5,1) * t150 - t236 - t86;
t235 = Icges(5,4) * t153;
t87 = -Icges(5,2) * t152 + Icges(5,6) * t223 + t235;
t240 = -Icges(5,1) * t152 - t235 - t87;
t115 = pkin(4) * t151 + pkin(6) * t150;
t58 = rSges(6,1) * t129 + rSges(6,2) * t128 + rSges(6,3) * t150;
t239 = t115 + t58;
t117 = pkin(4) * t153 + pkin(6) * t152;
t59 = rSges(6,1) * t131 + rSges(6,2) * t130 + rSges(6,3) * t152;
t238 = -t117 - t59;
t185 = (-pkin(4) * t171 + pkin(6) * t172) * t174;
t149 = qJD(4) * t185;
t65 = rSges(6,1) * t133 + rSges(6,2) * t132 + rSges(6,3) * t207;
t237 = t149 + t65;
t234 = Icges(5,4) * t171;
t233 = Icges(5,4) * t172;
t105 = rSges(6,1) * t162 + rSges(6,2) * t161 + rSges(6,3) * t229;
t160 = (pkin(4) * t172 + pkin(6) * t171) * t174;
t222 = t105 + t160;
t136 = -Icges(5,6) * t176 + (-Icges(5,2) * t171 + t233) * t174;
t156 = (-Icges(5,1) * t171 - t233) * t174;
t221 = t136 - t156;
t137 = -Icges(5,5) * t176 + (Icges(5,1) * t172 - t234) * t174;
t155 = (-Icges(5,2) * t172 - t234) * t174;
t220 = t137 + t155;
t218 = qJD(3) * t174;
t219 = qJD(2) * t175 + t177 * t218;
t217 = qJD(4) * t172;
t213 = qJD(5) * t140;
t212 = qJD(5) * t142;
t211 = qJD(5) * t172;
t204 = t213 / 0.2e1;
t203 = t212 / 0.2e1;
t201 = -qJD(2) * t177 + t175 * t218;
t200 = -qJD(3) * t176 + qJD(1);
t199 = t207 / 0.2e1;
t198 = -rSges(6,1) * t179 + rSges(6,2) * t178;
t28 = t161 * t54 + t162 * t56 + t229 * t52;
t29 = t161 * t55 + t162 * t57 + t229 * t53;
t195 = t140 * t28 + t142 * t29;
t194 = -t140 * t59 + t142 * t58;
t138 = -rSges(5,3) * t176 + (rSges(5,1) * t172 - rSges(5,2) * t171) * t174;
t90 = rSges(5,1) * t151 - rSges(5,2) * t150 + rSges(5,3) * t226;
t50 = (t138 * t226 + t176 * t90) * qJD(4) + t219;
t91 = rSges(5,1) * t153 - rSges(5,2) * t152 + rSges(5,3) * t223;
t51 = (-t138 * t223 - t176 * t91) * qJD(4) + t201;
t193 = t175 * t50 - t177 * t51;
t184 = (-rSges(5,1) * t171 - rSges(5,2) * t172) * t174;
t146 = qJD(4) * t184;
t98 = -rSges(5,1) * t139 - rSges(5,2) * t140;
t60 = (t146 * t226 + t176 * t98) * qJD(4);
t99 = -rSges(5,1) * t141 - rSges(5,2) * t142;
t61 = (-t146 * t223 - t176 * t99) * qJD(4);
t192 = t175 * t60 - t177 * t61;
t191 = -t175 * t91 + t177 * t90;
t190 = -t175 * t99 + t177 * t98;
t189 = -Icges(6,1) * t179 + Icges(6,4) * t178;
t188 = -Icges(6,4) * t179 + Icges(6,2) * t178;
t187 = -Icges(6,5) * t179 + Icges(6,6) * t178;
t186 = qJD(5) * t199;
t183 = (Icges(6,5) * t161 - Icges(6,6) * t162) * t163 + t134 * (Icges(6,5) * t128 - Icges(6,6) * t129) + t135 * (Icges(6,5) * t130 - Icges(6,6) * t131);
t182 = (Icges(6,1) * t130 - t231 - t55) * t135 + (Icges(6,1) * t128 - t232 - t54) * t134 + (Icges(6,1) * t161 - t103 - t230) * t163;
t181 = (-Icges(6,2) * t131 + t123 + t57) * t135 + (-Icges(6,2) * t129 + t122 + t56) * t134 + (-Icges(6,2) * t162 + t104 + t157) * t163;
t180 = (Icges(6,3) * t153 + t152 * t187 + t178 * t55 - t179 * t57) * t135 + (Icges(6,3) * t151 + t150 * t187 + t178 * t54 - t179 * t56) * t134 + (t103 * t178 - t104 * t179 + (Icges(6,3) * t172 + t171 * t187) * t174) * t163;
t145 = qJD(4) * t156;
t144 = qJD(4) * t155;
t127 = (rSges(6,3) * t172 + t171 * t198) * t174;
t126 = (Icges(6,5) * t172 + t171 * t189) * t174;
t125 = (Icges(6,6) * t172 + t171 * t188) * t174;
t121 = rSges(6,1) * t161 - rSges(6,2) * t162;
t116 = -pkin(4) * t152 + pkin(6) * t153;
t114 = -pkin(4) * t150 + pkin(6) * t151;
t113 = -rSges(5,1) * t152 - rSges(5,2) * t153;
t112 = -rSges(5,1) * t150 - rSges(5,2) * t151;
t107 = -Icges(5,5) * t152 - Icges(5,6) * t153;
t106 = -Icges(5,5) * t150 - Icges(5,6) * t151;
t97 = -Icges(5,1) * t141 - Icges(5,4) * t142;
t96 = -Icges(5,1) * t139 - Icges(5,4) * t140;
t95 = -Icges(5,4) * t141 - Icges(5,2) * t142;
t94 = -Icges(5,4) * t139 - Icges(5,2) * t140;
t93 = -Icges(5,5) * t141 - Icges(5,6) * t142;
t92 = -Icges(5,5) * t139 - Icges(5,6) * t140;
t81 = rSges(6,1) * t130 - rSges(6,2) * t131;
t80 = rSges(6,1) * t128 - rSges(6,2) * t129;
t73 = t153 * rSges(6,3) + t152 * t198;
t72 = t151 * rSges(6,3) + t150 * t198;
t71 = Icges(6,5) * t153 + t152 * t189;
t70 = Icges(6,5) * t151 + t150 * t189;
t69 = Icges(6,6) * t153 + t152 * t188;
t68 = Icges(6,6) * t151 + t150 * t188;
t49 = t190 * t216;
t48 = t191 * t216 + t200;
t38 = t102 * t229 + t103 * t161 + t104 * t162;
t33 = -t105 * t135 + t163 * t59 + (-t117 * t176 - t160 * t223) * qJD(4) + t201;
t32 = t105 * t134 - t163 * t58 + (t115 * t176 + t160 * t226) * qJD(4) + t219;
t23 = -t134 * t59 + t135 * t58 + (t115 * t177 - t117 * t175) * t216 + t200;
t18 = t103 * t132 + t104 * t133 + t161 * t63 + t162 * t64 + (t102 * t217 + t171 * t62) * t174;
t17 = -t105 * t212 - t135 * t65 + t163 * t47 + (-t101 * t176 + (-t149 * t177 + t211 * t59) * t174) * qJD(4);
t16 = t105 * t213 + t134 * t65 - t163 * t46 + (t100 * t176 + (t149 * t175 - t211 * t58) * t174) * qJD(4);
t13 = -t134 * t47 + t135 * t46 + t194 * qJD(5) + (t100 * t177 - t101 * t175) * t216;
t12 = t132 * t55 + t133 * t57 + t161 * t43 + t162 * t45 + (t171 * t41 + t217 * t53) * t174;
t11 = t132 * t54 + t133 * t56 + t161 * t42 + t162 * t44 + (t171 * t40 + t217 * t52) * t174;
t10 = t134 * t28 + t135 * t29 + t163 * t38;
t5 = t134 * t21 + t135 * t22 + t163 * t35;
t4 = t134 * t19 + t135 * t20 + t163 * t34;
t3 = t11 * t134 + t12 * t135 + t163 * t18 + (t207 * t38 + t195) * qJD(5);
t1 = [m(5) * t49 + m(6) * t13; m(5) * t192 + m(6) * (t16 * t175 - t17 * t177); m(5) * (-t176 * t49 + (t175 * t61 + t177 * t60) * t174) + m(6) * (-t13 * t176 + (t16 * t177 + t17 * t175) * t174); (-t176 * t34 + (t175 * t19 + t177 * t20) * t174) * t204 + (-t176 * t35 + (t175 * t21 + t177 * t22) * t174) * t203 + (-t176 * t18 + (t11 * t175 + t12 * t177) * t174) * t247 + ((t161 * t69 + t162 * t71) * t135 + (t161 * t68 + t162 * t70) * t134 + (t161 * t125 + t162 * t126) * t163 + (t151 * t28 + t153 * t29) * qJD(5) + ((qJD(5) * t38 + t102 * t163 + t134 * t52 + t135 * t53) * t172 + t180 * t171) * t174) * t248 + (-t15 * t176 + (t175 * t8 + t177 * t9) * t174) * t249 + ((t130 * t69 + t131 * t71 + t153 * t53) * t135 + (t130 * t68 + t131 * t70 + t153 * t52) * t134 + (t153 * t102 + t130 * t125 + t131 * t126) * t163 + (t151 * t21 + t153 * t22 + t228 * t35) * qJD(5) + t180 * t152) * t250 + (-t14 * t176 + (t175 * t6 + t177 * t7) * t174) * t251 + ((t128 * t69 + t129 * t71 + t151 * t53) * t135 + (t128 * t68 + t129 * t70 + t151 * t52) * t134 + (t151 * t102 + t128 * t125 + t129 * t126) * t163 + (t151 * t19 + t153 * t20 + t228 * t34) * qJD(5) + t180 * t150) * t252 + t223 * t253 + t226 * t254 + (t176 ^ 2 * t143 + (((t175 * t241 + t177 * t240) * t172 + (t175 * t243 + t177 * t242) * t171) * t174 + (-t106 * t175 - t107 * t177 + t171 * t220 + t172 * t221) * t176) * t216) * t215 / 0.2e1 + (-t176 * (-t136 * t140 - t137 * t139 + t143 * t226 - t144 * t150 + t145 * t151) + (t175 * (-t139 * t88 - t140 * t86 - t150 * t94 + t151 * t96 + t226 * t92) + t177 * (-t139 * t89 - t140 * t87 - t150 * t95 + t151 * t97 + t226 * t93)) * t174) * t206 + (-t176 * (-t136 * t142 - t137 * t141 + t143 * t223 - t144 * t152 + t145 * t153) + (t175 * (-t141 * t88 - t142 * t86 - t152 * t94 + t153 * t96 + t223 * t92) + t177 * (-t141 * t89 - t142 * t87 - t152 * t95 + t153 * t97 + t223 * t93)) * t174) * t205 - (-t176 * (-t143 * t176 + (-t144 * t171 + t145 * t172 + (-t136 * t172 - t137 * t171) * qJD(4)) * t174) + (t175 * (-t176 * t92 + (-t171 * t94 + t172 * t96 + (-t171 * t88 - t172 * t86) * qJD(4)) * t174) + t177 * (-t176 * t93 + (-t171 * t95 + t172 * t97 + (-t171 * t89 - t172 * t87) * qJD(4)) * t174)) * t174) * t215 + (-t176 * t38 + (t175 * t28 + t177 * t29) * t174) * t186 - t172 * t10 * t210 / 0.2e1 - t176 * t3 / 0.2e1 - (t151 * t4 + t153 * t5) * qJD(5) / 0.2e1 - (t177 * ((t107 * t223 + t152 * t242 + t153 * t240) * t223 + (t106 * t223 + t152 * t243 + t153 * t241) * t226 - (-t152 * t220 - t153 * t221 + t154 * t223) * t176) + t175 * ((t107 * t226 + t150 * t242 + t151 * t240) * t223 + (t106 * t226 + t150 * t243 + t151 * t241) * t226 - (-t150 * t220 - t151 * t221 + t154 * t226) * t176)) * t174 * qJD(4) ^ 2 / 0.2e1 + ((t16 * t239 + t17 * t238 + t244 * t33 + t245 * t32) * t176 + ((t13 * t239 - t17 * t222 + t23 * t245 - t237 * t33) * t177 + (t13 * t238 + t16 * t222 + t23 * t244 + t237 * t32) * t175) * t174 - t32 * (t127 * t134 - t163 * t72) - t33 * (-t127 * t135 + t163 * t73) - t23 * (-t134 * t73 + t135 * t72) - (t32 * (t105 * t151 - t228 * t58) + t33 * (-t105 * t153 + t228 * t59) + t23 * (-t151 * t59 + t153 * t58)) * qJD(5) - ((t32 * t114 - t33 * t116) * t176 + (t23 * (t114 * t177 - t116 * t175) + (t175 * t32 - t177 * t33) * t185) * t174) * qJD(4)) * m(6) + ((t50 * t98 - t51 * t99 + t60 * t90 - t61 * t91) * t176 + (t138 * t192 + t146 * t193 + t190 * t48 + t191 * t49) * t174 - ((t50 * t112 - t51 * t113) * t176 + (t48 * (t112 * t177 - t113 * t175) + t193 * t184) * t174) * qJD(4)) * m(5); t142 * t5 / 0.2e1 + t152 * t253 + (t150 * t21 + t152 * t22 + t229 * t35) * t203 + (t150 * t8 + t152 * t9 + (t15 * t171 + t217 * t35) * t174 + t196) * t249 + t140 * t4 / 0.2e1 + t150 * t254 + (t150 * t19 + t152 * t20 + t229 * t34) * t204 + (t150 * t6 + t152 * t7 + (t14 * t171 + t217 * t34) * t174 + t197) * t251 + t10 * t199 + t3 * t229 / 0.2e1 + (t150 * t28 + t152 * t29 + t229 * t38) * t186 + (t11 * t150 + t12 * t152 + (t171 * t18 + t217 * t38) * t174 + t195) * t247 + (t130 * t181 + t131 * t182 + t152 * t183) * t250 + (t128 * t181 + t129 * t182 + t150 * t183) * t252 + (t161 * t181 + t162 * t182 + t183 * t229) * t248 + (t13 * (-t150 * t59 + t152 * t58) + (t150 * t32 - t152 * t33) * t65 + (t140 * t32 - t142 * t33 + t150 * t16 - t152 * t17) * t105 + ((-t32 * t58 + t33 * t59) * t217 + (-t16 * t58 + t17 * t59 - t32 * t46 + t33 * t47) * t171) * t174 - t32 * (t121 * t134 - t163 * t80) - t33 * (-t121 * t135 + t163 * t81) + (t134 * t81 - t135 * t80 - t150 * t47 + t152 * t46 + t194) * t23) * m(6);];
tauc = t1(:);
