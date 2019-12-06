% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPPRR2
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPPRR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPPRR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:26
% EndTime: 2019-12-05 14:59:35
% DurationCPUTime: 6.00s
% Computational Cost: add. (11982->536), mult. (33175->789), div. (0->0), fcn. (39660->10), ass. (0->210)
t180 = sin(pkin(9));
t182 = cos(pkin(9));
t234 = sin(pkin(7));
t235 = cos(pkin(8));
t212 = t235 * t234;
t236 = cos(pkin(7));
t168 = -t180 * t236 + t182 * t212;
t181 = sin(pkin(8));
t216 = t181 * t234;
t245 = sin(qJ(4));
t246 = cos(qJ(4));
t190 = -t168 * t245 + t216 * t246;
t167 = t180 * t212 + t182 * t236;
t224 = qJD(4) * t167;
t134 = -qJD(5) * t190 + t224;
t251 = t134 / 0.2e1;
t213 = t236 * t235;
t170 = t180 * t234 + t182 * t213;
t217 = t181 * t236;
t191 = -t170 * t245 + t217 * t246;
t169 = t180 * t213 - t182 * t234;
t223 = qJD(4) * t169;
t135 = -qJD(5) * t191 + t223;
t249 = t135 / 0.2e1;
t226 = t181 * t182;
t171 = t226 * t245 + t235 * t246;
t227 = t181 * t180;
t221 = qJD(4) * t227;
t163 = qJD(5) * t171 + t221;
t247 = t163 / 0.2e1;
t243 = qJD(5) / 0.2e1;
t172 = t246 * t226 - t235 * t245;
t158 = t168 * t246 + t216 * t245;
t183 = sin(qJ(5));
t184 = cos(qJ(5));
t130 = -t158 * t183 + t167 * t184;
t131 = t158 * t184 + t167 * t183;
t151 = t158 * qJD(4);
t162 = t172 * t184 + t183 * t227;
t164 = t171 * qJD(4);
t128 = -qJD(5) * t162 + t164 * t183;
t161 = -t172 * t183 + t184 * t227;
t129 = qJD(5) * t161 - t164 * t184;
t165 = t172 * qJD(4);
t62 = Icges(6,5) * t129 + Icges(6,6) * t128 + Icges(6,3) * t165;
t63 = Icges(6,4) * t129 + Icges(6,2) * t128 + Icges(6,6) * t165;
t64 = Icges(6,1) * t129 + Icges(6,4) * t128 + Icges(6,5) * t165;
t150 = t190 * qJD(4);
t82 = -qJD(5) * t131 - t150 * t183;
t83 = qJD(5) * t130 + t150 * t184;
t92 = Icges(6,5) * t162 + Icges(6,6) * t161 + Icges(6,3) * t171;
t228 = Icges(6,4) * t162;
t93 = Icges(6,2) * t161 + Icges(6,6) * t171 + t228;
t156 = Icges(6,4) * t161;
t94 = Icges(6,1) * t162 + Icges(6,5) * t171 + t156;
t14 = t130 * t63 + t131 * t64 + t151 * t92 - t190 * t62 + t82 * t93 + t83 * t94;
t160 = t170 * t246 + t217 * t245;
t153 = t160 * qJD(4);
t52 = Icges(6,5) * t131 + Icges(6,6) * t130 - Icges(6,3) * t190;
t230 = Icges(6,4) * t131;
t54 = Icges(6,2) * t130 - Icges(6,6) * t190 + t230;
t126 = Icges(6,4) * t130;
t56 = Icges(6,1) * t131 - Icges(6,5) * t190 + t126;
t19 = t130 * t54 + t131 * t56 - t190 * t52;
t132 = -t160 * t183 + t169 * t184;
t133 = t160 * t184 + t169 * t183;
t53 = Icges(6,5) * t133 + Icges(6,6) * t132 - Icges(6,3) * t191;
t229 = Icges(6,4) * t133;
t55 = Icges(6,2) * t132 - Icges(6,6) * t191 + t229;
t127 = Icges(6,4) * t132;
t57 = Icges(6,1) * t133 - Icges(6,5) * t191 + t127;
t20 = t130 * t55 + t131 * t57 - t190 * t53;
t34 = t130 * t93 + t131 * t94 - t190 * t92;
t196 = t19 * t151 + t153 * t20 + t165 * t34;
t40 = Icges(6,5) * t83 + Icges(6,6) * t82 + Icges(6,3) * t151;
t42 = Icges(6,4) * t83 + Icges(6,2) * t82 + Icges(6,6) * t151;
t44 = Icges(6,1) * t83 + Icges(6,4) * t82 + Icges(6,5) * t151;
t7 = t130 * t42 + t131 * t44 + t151 * t52 - t190 * t40 + t54 * t82 + t56 * t83;
t152 = t191 * qJD(4);
t84 = -qJD(5) * t133 - t152 * t183;
t85 = qJD(5) * t132 + t152 * t184;
t41 = Icges(6,5) * t85 + Icges(6,6) * t84 + Icges(6,3) * t153;
t43 = Icges(6,4) * t85 + Icges(6,2) * t84 + Icges(6,6) * t153;
t45 = Icges(6,1) * t85 + Icges(6,4) * t84 + Icges(6,5) * t153;
t8 = t130 * t43 + t131 * t45 + t151 * t53 - t190 * t41 + t55 * t82 + t57 * t83;
t255 = t14 * t247 + t196 * t243 + t8 * t249 + t7 * t251;
t10 = t132 * t43 + t133 * t45 + t153 * t53 - t191 * t41 + t55 * t84 + t57 * t85;
t15 = t132 * t63 + t133 * t64 + t153 * t92 - t191 * t62 + t84 * t93 + t85 * t94;
t21 = t132 * t54 + t133 * t56 - t191 * t52;
t22 = t132 * t55 + t133 * t57 - t191 * t53;
t35 = t132 * t93 + t133 * t94 - t191 * t92;
t195 = t21 * t151 + t22 * t153 + t35 * t165;
t9 = t132 * t42 + t133 * t44 + t153 * t52 - t191 * t40 + t54 * t84 + t56 * t85;
t254 = t10 * t249 + t15 * t247 + t195 * t243 + t9 * t251;
t11 = t128 * t54 + t129 * t56 + t161 * t42 + t162 * t44 + t165 * t52 + t171 * t40;
t12 = t128 * t55 + t129 * t57 + t161 * t43 + t162 * t45 + t165 * t53 + t171 * t41;
t18 = t128 * t93 + t129 * t94 + t161 * t63 + t162 * t64 + t165 * t92 + t171 * t62;
t28 = t161 * t54 + t162 * t56 + t171 * t52;
t29 = t161 * t55 + t162 * t57 + t171 * t53;
t38 = t161 * t93 + t162 * t94 + t171 * t92;
t194 = t28 * t151 + t29 * t153 + t38 * t165;
t253 = t11 * t251 + t12 * t249 + t18 * t247 + t194 * t243;
t252 = -t134 / 0.2e1;
t250 = -t135 / 0.2e1;
t248 = -t163 / 0.2e1;
t104 = pkin(4) * t150 + pkin(6) * t151;
t46 = rSges(6,1) * t83 + rSges(6,2) * t82 + rSges(6,3) * t151;
t242 = t104 + t46;
t105 = pkin(4) * t152 + pkin(6) * t153;
t47 = rSges(6,1) * t85 + rSges(6,2) * t84 + rSges(6,3) * t153;
t241 = t105 + t47;
t115 = pkin(4) * t158 - pkin(6) * t190;
t58 = rSges(6,1) * t131 + rSges(6,2) * t130 - rSges(6,3) * t190;
t240 = t115 + t58;
t117 = pkin(4) * t160 - pkin(6) * t191;
t59 = rSges(6,1) * t133 + rSges(6,2) * t132 - rSges(6,3) * t191;
t239 = t117 + t59;
t143 = -pkin(4) * t164 + pkin(6) * t165;
t65 = rSges(6,1) * t129 + rSges(6,2) * t128 + rSges(6,3) * t165;
t238 = t143 + t65;
t149 = pkin(4) * t172 + pkin(6) * t171;
t95 = rSges(6,1) * t162 + rSges(6,2) * t161 + rSges(6,3) * t171;
t237 = t149 + t95;
t233 = Icges(5,4) * t158;
t232 = Icges(5,4) * t160;
t231 = Icges(5,4) * t172;
t225 = qJD(2) * t234 + qJD(3) * t217;
t220 = t151 * t243;
t219 = t153 * t243;
t218 = t165 * t243;
t211 = -rSges(6,1) * t184 + rSges(6,2) * t183;
t210 = -t151 * t59 + t153 * t58;
t209 = t151 * t95 - t165 * t58;
t208 = -t153 * t95 + t165 * t59;
t90 = rSges(5,1) * t158 + rSges(5,2) * t190 + rSges(5,3) * t167;
t91 = rSges(5,1) * t160 + rSges(5,2) * t191 + rSges(5,3) * t169;
t207 = -t167 * t91 + t169 * t90;
t206 = -qJD(2) * t236 + qJD(3) * t216;
t205 = -Icges(6,1) * t184 + Icges(6,4) * t183;
t204 = -Icges(6,4) * t184 + Icges(6,2) * t183;
t203 = -Icges(6,5) * t184 + Icges(6,6) * t183;
t102 = rSges(5,1) * t150 - rSges(5,2) * t151;
t103 = rSges(5,1) * t152 - rSges(5,2) * t153;
t202 = t102 * t169 - t103 * t167;
t201 = -qJD(3) * t235 + qJD(1);
t138 = rSges(5,1) * t172 - rSges(5,2) * t171 + rSges(5,3) * t227;
t200 = t138 * t167 - t227 * t90;
t199 = -t138 * t169 + t227 * t91;
t142 = -rSges(5,1) * t164 - rSges(5,2) * t165;
t198 = -t102 * t227 + t142 * t167;
t197 = t103 * t227 - t142 * t169;
t193 = (Icges(6,5) * t161 - Icges(6,6) * t162) * t163 + t134 * (Icges(6,5) * t130 - Icges(6,6) * t131) + t135 * (Icges(6,5) * t132 - Icges(6,6) * t133);
t192 = (Icges(5,5) * t190 - Icges(5,6) * t158) * t167 + (Icges(5,5) * t191 - Icges(5,6) * t160) * t169 + (-Icges(5,5) * t171 - Icges(5,6) * t172) * t227;
t189 = (Icges(6,1) * t132 - t229 - t55) * t135 + (Icges(6,1) * t130 - t230 - t54) * t134 + (Icges(6,1) * t161 - t228 - t93) * t163;
t188 = (-Icges(6,2) * t133 + t127 + t57) * t135 + (-Icges(6,2) * t131 + t126 + t56) * t134 + (-Icges(6,2) * t162 + t156 + t94) * t163;
t166 = Icges(5,4) * t171;
t137 = Icges(5,1) * t172 + Icges(5,5) * t227 - t166;
t154 = Icges(5,4) * t190;
t155 = Icges(5,4) * t191;
t88 = Icges(5,1) * t158 + Icges(5,5) * t167 + t154;
t89 = Icges(5,1) * t160 + Icges(5,5) * t169 + t155;
t187 = (Icges(5,2) * t160 - t155 - t89) * t169 + (Icges(5,2) * t158 - t154 - t88) * t167 + (Icges(5,2) * t172 - t137 + t166) * t227;
t136 = -Icges(5,2) * t171 + Icges(5,6) * t227 + t231;
t86 = Icges(5,2) * t190 + Icges(5,6) * t167 + t233;
t87 = Icges(5,2) * t191 + Icges(5,6) * t169 + t232;
t186 = (Icges(5,1) * t191 - t232 - t87) * t169 + (Icges(5,1) * t190 - t233 - t86) * t167 + (-Icges(5,1) * t171 - t136 - t231) * t227;
t185 = (Icges(6,3) * t160 + t183 * t55 - t184 * t57 - t191 * t203) * t135 + (Icges(6,3) * t158 + t183 * t54 - t184 * t56 - t190 * t203) * t134 + (Icges(6,3) * t172 + t171 * t203 + t183 * t93 - t184 * t94) * t163;
t148 = -pkin(4) * t171 + pkin(6) * t172;
t147 = -rSges(5,1) * t171 - rSges(5,2) * t172;
t141 = -Icges(5,1) * t164 - Icges(5,4) * t165;
t140 = -Icges(5,4) * t164 - Icges(5,2) * t165;
t139 = -Icges(5,5) * t164 - Icges(5,6) * t165;
t125 = t172 * rSges(6,3) + t171 * t211;
t124 = Icges(6,5) * t172 + t171 * t205;
t123 = Icges(6,6) * t172 + t171 * t204;
t121 = rSges(6,1) * t161 - rSges(6,2) * t162;
t116 = pkin(4) * t191 + pkin(6) * t160;
t114 = pkin(4) * t190 + pkin(6) * t158;
t113 = rSges(5,1) * t191 - rSges(5,2) * t160;
t112 = rSges(5,1) * t190 - rSges(5,2) * t158;
t101 = Icges(5,1) * t152 - Icges(5,4) * t153;
t100 = Icges(5,1) * t150 - Icges(5,4) * t151;
t99 = Icges(5,4) * t152 - Icges(5,2) * t153;
t98 = Icges(5,4) * t150 - Icges(5,2) * t151;
t97 = Icges(5,5) * t152 - Icges(5,6) * t153;
t96 = Icges(5,5) * t150 - Icges(5,6) * t151;
t81 = t160 * rSges(6,3) - t191 * t211;
t80 = t158 * rSges(6,3) - t190 * t211;
t79 = Icges(6,5) * t160 - t191 * t205;
t78 = Icges(6,5) * t158 - t190 * t205;
t77 = Icges(6,6) * t160 - t191 * t204;
t76 = Icges(6,6) * t158 - t190 * t204;
t73 = rSges(6,1) * t132 - rSges(6,2) * t133;
t72 = rSges(6,1) * t130 - rSges(6,2) * t131;
t61 = t197 * qJD(4);
t60 = t198 * qJD(4);
t51 = qJD(4) * t199 + t206;
t50 = qJD(4) * t200 + t225;
t49 = t202 * qJD(4);
t48 = qJD(4) * t207 + t201;
t33 = -t135 * t95 + t163 * t59 + (t117 * t227 - t149 * t169) * qJD(4) + t206;
t32 = t134 * t95 - t163 * t58 + (-t115 * t227 + t149 * t167) * qJD(4) + t225;
t23 = -t134 * t59 + t135 * t58 + (t115 * t169 - t117 * t167) * qJD(4) + t201;
t17 = -t135 * t65 + t163 * t47 + t208 * qJD(5) + (t105 * t227 - t143 * t169) * qJD(4);
t16 = t134 * t65 - t163 * t46 + t209 * qJD(5) + (-t104 * t227 + t143 * t167) * qJD(4);
t13 = -t134 * t47 + t135 * t46 + t210 * qJD(5) + (t104 * t169 - t105 * t167) * qJD(4);
t6 = t134 * t28 + t135 * t29 + t38 * t163;
t5 = t134 * t21 + t135 * t22 + t163 * t35;
t4 = t19 * t134 + t135 * t20 + t163 * t34;
t1 = [m(5) * t49 + m(6) * t13; m(5) * (t234 * t60 - t236 * t61) + m(6) * (t16 * t234 - t17 * t236); m(5) * (-t49 * t235 + (t234 * t61 + t236 * t60) * t181) + m(6) * (-t13 * t235 + (t16 * t236 + t17 * t234) * t181); (t167 * t19 + t169 * t20 + t227 * t34) * t220 + ((-t136 * t153 + t137 * t152 + t139 * t169 + t140 * t191 + t141 * t160) * t227 + t167 * (t100 * t160 + t152 * t88 - t153 * t86 + t169 * t96 + t191 * t98) + t169 * (t101 * t160 + t152 * t89 - t153 * t87 + t169 * t97 + t191 * t99)) * t223 + ((-t136 * t151 + t137 * t150 + t139 * t167 + t140 * t190 + t141 * t158) * t227 + t167 * (t100 * t158 + t150 * t88 - t151 * t86 + t167 * t96 + t190 * t98) + t169 * (t101 * t158 + t150 * t89 - t151 * t87 + t167 * t97 + t190 * t99)) * t224 + ((-t136 * t165 - t137 * t164 + t139 * t227 - t140 * t171 + t141 * t172) * t227 + t167 * (t100 * t172 - t164 * t88 - t165 * t86 - t171 * t98 + t227 * t96) + t169 * (t101 * t172 - t164 * t89 - t165 * t87 - t171 * t99 + t227 * t97)) * t221 + (t167 * t28 + t169 * t29 + t227 * t38) * t218 + (t167 * t21 + t169 * t22 + t227 * t35) * t219 + (t11 * t167 + t12 * t169 + t18 * t227) * t247 + ((t161 * t77 + t162 * t79 + t172 * t53) * t135 + (t161 * t76 + t162 * t78 + t172 * t52) * t134 + (t161 * t123 + t162 * t124 + t172 * t92) * t163 + (t158 * t28 + t160 * t29 + t172 * t38) * qJD(5) + t185 * t171) * t248 + (t10 * t169 + t15 * t227 + t167 * t9) * t249 + ((t132 * t77 + t133 * t79 + t160 * t53) * t135 + (t132 * t76 + t133 * t78 + t160 * t52) * t134 + (t132 * t123 + t133 * t124 + t160 * t92) * t163 + (t158 * t21 + t160 * t22 + t172 * t35) * qJD(5) - t185 * t191) * t250 + (t14 * t227 + t167 * t7 + t169 * t8) * t251 + ((t130 * t77 + t131 * t79 + t158 * t53) * t135 + (t130 * t76 + t131 * t78 + t158 * t52) * t134 + (t130 * t123 + t131 * t124 + t158 * t92) * t163 + (t158 * t19 + t160 * t20 + t172 * t34) * qJD(5) - t185 * t190) * t252 + t227 * t253 + t169 * t254 + t167 * t255 + (-t32 * (t125 * t134 - t163 * t80) - t33 * (-t125 * t135 + t163 * t81) - t23 * (-t134 * t81 + t135 * t80) - (t32 * (t158 * t95 - t172 * t58) + t33 * (-t160 * t95 + t172 * t59) + t23 * (-t158 * t59 + t160 * t58)) * qJD(5) - (t32 * (-t114 * t227 + t148 * t167) + t33 * (t116 * t227 - t148 * t169) + t23 * (t114 * t169 - t116 * t167)) * qJD(4) + (-t16 * t240 + t17 * t239 + t241 * t33 - t242 * t32) * t227 + (t13 * t240 - t17 * t237 + t23 * t242 - t238 * t33) * t169 + (-t13 * t239 + t16 * t237 - t23 * t241 + t238 * t32) * t167) * m(6) + (-(t50 * (-t112 * t227 + t147 * t167) + t51 * (t113 * t227 - t147 * t169) + t48 * (t112 * t169 - t113 * t167)) * qJD(4) + t197 * t51 + t198 * t50 + t199 * t61 + t200 * t60 + t202 * t48 + t207 * t49) * m(5) - (t158 * t4 + t160 * t5 + t172 * t6) * qJD(5) / 0.2e1 - ((t171 * t187 + t172 * t186 + t192 * t227) * t227 + t169 * (t160 * t186 + t169 * t192 - t187 * t191) + t167 * (t158 * t186 + t167 * t192 - t187 * t190)) * qJD(4) ^ 2 / 0.2e1; t153 * t5 / 0.2e1 - t191 * t254 + (t171 * t35 - t190 * t21 - t191 * t22) * t219 + (-t10 * t191 + t15 * t171 - t190 * t9 + t195) * t249 + t151 * t4 / 0.2e1 - t190 * t255 + (t171 * t34 - t19 * t190 - t191 * t20) * t220 + (t14 * t171 - t190 * t7 - t191 * t8 + t196) * t251 + t165 * t6 / 0.2e1 + t171 * t253 + (t171 * t38 - t190 * t28 - t191 * t29) * t218 + (-t11 * t190 - t12 * t191 + t18 * t171 + t194) * t247 + (t132 * t188 + t133 * t189 - t191 * t193) * t250 + (t130 * t188 + t131 * t189 - t190 * t193) * t252 + (t161 * t188 + t162 * t189 + t171 * t193) * t248 + (t16 * (-t171 * t58 - t190 * t95) + t17 * (t171 * t59 + t191 * t95) + t13 * (t190 * t59 - t191 * t58) + (t121 * t135 - t163 * t73 + t171 * t47 + t191 * t65 + t208) * t33 + (-t121 * t134 + t163 * t72 - t171 * t46 - t190 * t65 + t209) * t32 + (t134 * t73 - t135 * t72 + t190 * t47 - t191 * t46 + t210) * t23) * m(6);];
tauc = t1(:);
