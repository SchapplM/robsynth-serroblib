% Calculate time derivative of joint inertia matrix for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR6_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR6_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:40
% EndTime: 2019-12-31 16:34:48
% DurationCPUTime: 4.66s
% Computational Cost: add. (11207->452), mult. (18147->717), div. (0->0), fcn. (17490->8), ass. (0->239)
t186 = sin(pkin(7));
t187 = cos(pkin(7));
t290 = t187 ^ 2;
t297 = t186 ^ 2 + t290;
t189 = sin(qJ(2));
t191 = cos(qJ(2));
t296 = qJD(2) * (rSges(3,1) * t189 + rSges(3,2) * t191);
t190 = cos(qJ(3));
t285 = pkin(3) * t190;
t146 = -pkin(6) * t191 + t285 * t189;
t293 = 2 * m(4);
t292 = 2 * m(5);
t289 = t186 / 0.2e1;
t288 = -t187 / 0.2e1;
t287 = t187 / 0.2e1;
t286 = -t191 / 0.2e1;
t188 = sin(qJ(3));
t247 = qJD(3) * t188;
t244 = pkin(3) * t247;
t193 = -t146 * qJD(2) - t191 * t244;
t245 = qJD(3) * t190;
t243 = pkin(3) * t245;
t185 = qJ(3) + qJ(4);
t182 = sin(t185);
t183 = cos(t185);
t184 = qJD(3) + qJD(4);
t249 = qJD(2) * t189;
t202 = -t184 * t186 + t187 * t249;
t263 = t187 * t191;
t240 = t184 * t263;
t135 = t202 * t182 - t183 * t240;
t136 = -t182 * t240 - t202 * t183;
t248 = qJD(2) * t191;
t236 = t187 * t248;
t80 = rSges(5,1) * t136 + rSges(5,2) * t135 + rSges(5,3) * t236;
t283 = t186 * t243 + t193 * t187 + t80;
t266 = t186 * t191;
t157 = -t182 * t266 - t183 * t187;
t158 = -t182 * t187 + t183 * t266;
t267 = t186 * t189;
t110 = Icges(5,5) * t158 + Icges(5,6) * t157 + Icges(5,3) * t267;
t112 = Icges(5,4) * t158 + Icges(5,2) * t157 + Icges(5,6) * t267;
t114 = Icges(5,1) * t158 + Icges(5,4) * t157 + Icges(5,5) * t267;
t159 = -t182 * t263 + t183 * t186;
t160 = t182 * t186 + t183 * t263;
t264 = t187 * t189;
t47 = t110 * t264 + t112 * t159 + t114 * t160;
t282 = t186 * t47;
t262 = t188 * t191;
t173 = -t186 * t262 - t187 * t190;
t261 = t190 * t191;
t265 = t187 * t188;
t174 = t186 * t261 - t265;
t119 = Icges(4,5) * t174 + Icges(4,6) * t173 + Icges(4,3) * t267;
t121 = Icges(4,4) * t174 + Icges(4,2) * t173 + Icges(4,6) * t267;
t123 = Icges(4,1) * t174 + Icges(4,4) * t173 + Icges(4,5) * t267;
t175 = t186 * t190 - t187 * t262;
t268 = t186 * t188;
t176 = t187 * t261 + t268;
t56 = t119 * t264 + t121 * t175 + t123 * t176;
t281 = t186 * t56;
t111 = Icges(5,5) * t160 + Icges(5,6) * t159 + Icges(5,3) * t264;
t113 = Icges(5,4) * t160 + Icges(5,2) * t159 + Icges(5,6) * t264;
t115 = Icges(5,1) * t160 + Icges(5,4) * t159 + Icges(5,5) * t264;
t46 = t111 * t267 + t113 * t157 + t115 * t158;
t280 = t187 * t46;
t120 = Icges(4,5) * t176 + Icges(4,6) * t175 + Icges(4,3) * t264;
t122 = Icges(4,4) * t176 + Icges(4,2) * t175 + Icges(4,6) * t264;
t124 = Icges(4,1) * t176 + Icges(4,4) * t175 + Icges(4,5) * t264;
t55 = t120 * t267 + t122 * t173 + t124 * t174;
t279 = t187 * t55;
t116 = rSges(5,1) * t158 + rSges(5,2) * t157 + rSges(5,3) * t267;
t201 = t184 * t187 + t186 * t249;
t241 = t184 * t266;
t133 = t201 * t182 - t183 * t241;
t134 = -t182 * t241 - t201 * t183;
t237 = t186 * t248;
t79 = rSges(5,1) * t134 + rSges(5,2) * t133 + rSges(5,3) * t237;
t278 = t116 * t236 + t79 * t264;
t275 = Icges(4,4) * t188;
t274 = Icges(4,4) * t190;
t273 = Icges(5,4) * t182;
t272 = Icges(5,4) * t183;
t218 = Icges(5,5) * t183 - Icges(5,6) * t182;
t147 = -Icges(5,3) * t191 + t218 * t189;
t271 = t147 * t191;
t269 = t184 * t189;
t260 = t191 * ((-Icges(5,5) * t182 - Icges(5,6) * t183) * t269 + (Icges(5,3) * t189 + t218 * t191) * qJD(2));
t219 = Icges(4,5) * t190 - Icges(4,6) * t188;
t246 = qJD(3) * t189;
t259 = t191 * ((-Icges(4,5) * t188 - Icges(4,6) * t190) * t246 + (Icges(4,3) * t189 + t219 * t191) * qJD(2));
t161 = -Icges(4,3) * t191 + t219 * t189;
t258 = t191 * t161;
t228 = rSges(5,1) * t183 - rSges(5,2) * t182;
t109 = (-rSges(5,1) * t182 - rSges(5,2) * t183) * t269 + (rSges(5,3) * t189 + t228 * t191) * qJD(2);
t194 = pkin(6) * t189 + t285 * t191;
t137 = t194 * qJD(2) - t189 * t244;
t257 = -t109 - t137;
t127 = -pkin(3) * t265 + t194 * t186;
t256 = t116 + t127;
t117 = rSges(5,1) * t160 + rSges(5,2) * t159 + rSges(5,3) * t264;
t128 = pkin(3) * t268 + t194 * t187;
t255 = t117 + t128;
t229 = rSges(4,1) * t190 - rSges(4,2) * t188;
t132 = (-rSges(4,1) * t188 - rSges(4,2) * t190) * t246 + (rSges(4,3) * t189 + t229 * t191) * qJD(2);
t231 = pkin(2) * t191 + pkin(5) * t189;
t178 = t231 * qJD(2);
t254 = -t132 - t178;
t150 = -rSges(5,3) * t191 + t228 * t189;
t82 = t191 * t116 + t150 * t267;
t253 = -t146 - t150;
t180 = t189 * pkin(2) - t191 * pkin(5);
t252 = t297 * qJD(2) * t180;
t164 = -rSges(4,3) * t191 + t229 * t189;
t251 = -t164 - t180;
t250 = t297 * t231;
t242 = t109 * t267 + t150 * t237 + t191 * t79;
t239 = -t178 + t257;
t238 = -t180 + t253;
t235 = t188 * t249;
t234 = t190 * t249;
t232 = t255 * t191;
t217 = -t112 * t182 + t114 * t183;
t52 = -t110 * t191 + t217 * t189;
t216 = -t113 * t182 + t115 * t183;
t53 = -t111 * t191 + t216 * t189;
t227 = t52 * t186 + t53 * t187;
t215 = -t121 * t188 + t123 * t190;
t60 = -t119 * t191 + t215 * t189;
t214 = -t122 * t188 + t124 * t190;
t61 = -t120 * t191 + t214 * t189;
t226 = t60 * t186 + t61 * t187;
t224 = Icges(4,1) * t190 - t275;
t223 = Icges(5,1) * t183 - t273;
t221 = -Icges(4,2) * t188 + t274;
t220 = -Icges(5,2) * t182 + t272;
t125 = rSges(4,1) * t174 + rSges(4,2) * t173 + rSges(4,3) * t267;
t126 = rSges(4,1) * t176 + rSges(4,2) * t175 + rSges(4,3) * t264;
t213 = t125 * t187 - t126 * t186;
t148 = -Icges(5,6) * t191 + t220 * t189;
t149 = -Icges(5,5) * t191 + t223 * t189;
t212 = t148 * t182 - t149 * t183;
t162 = -Icges(4,6) * t191 + t221 * t189;
t163 = -Icges(4,5) * t191 + t224 * t189;
t210 = t162 * t188 - t163 * t190;
t73 = Icges(5,5) * t134 + Icges(5,6) * t133 + Icges(5,3) * t237;
t207 = t110 * t248 + t189 * t73;
t74 = Icges(5,5) * t136 + Icges(5,6) * t135 + Icges(5,3) * t236;
t206 = t111 * t248 + t189 * t74;
t142 = -t174 * qJD(3) + t186 * t235;
t143 = t173 * qJD(3) - t186 * t234;
t87 = Icges(4,5) * t143 + Icges(4,6) * t142 + Icges(4,3) * t237;
t205 = t119 * t248 + t189 * t87;
t144 = -t176 * qJD(3) + t187 * t235;
t145 = t175 * qJD(3) - t187 * t234;
t88 = Icges(4,5) * t145 + Icges(4,6) * t144 + Icges(4,3) * t236;
t204 = t120 * t248 + t189 * t88;
t107 = (-Icges(5,2) * t183 - t273) * t269 + (Icges(5,6) * t189 + t220 * t191) * qJD(2);
t108 = (-Icges(5,1) * t182 - t272) * t269 + (Icges(5,5) * t189 + t223 * t191) * qJD(2);
t75 = Icges(5,4) * t134 + Icges(5,2) * t133 + Icges(5,6) * t237;
t77 = Icges(5,1) * t134 + Icges(5,4) * t133 + Icges(5,5) * t237;
t17 = (t217 * qJD(2) - t73) * t191 + (qJD(2) * t110 + (-t112 * t184 + t77) * t183 + (-t114 * t184 - t75) * t182) * t189;
t76 = Icges(5,4) * t136 + Icges(5,2) * t135 + Icges(5,6) * t236;
t78 = Icges(5,1) * t136 + Icges(5,4) * t135 + Icges(5,5) * t236;
t18 = (t216 * qJD(2) - t74) * t191 + (qJD(2) * t111 + (-t113 * t184 + t78) * t183 + (-t115 * t184 - t76) * t182) * t189;
t45 = t110 * t267 + t112 * t157 + t114 * t158;
t48 = t111 * t264 + t113 * t159 + t115 * t160;
t19 = t112 * t133 + t114 * t134 + t157 * t75 + t158 * t77 + t207 * t186;
t20 = t113 * t133 + t115 * t134 + t157 * t76 + t158 * t78 + t206 * t186;
t63 = t147 * t267 + t148 * t157 + t149 * t158;
t5 = -(t107 * t157 + t108 * t158 + t133 * t148 + t134 * t149) * t191 + (t20 * t187 + (t19 - t260) * t186) * t189 + (t63 * t189 + (t280 + (t45 - t271) * t186) * t191) * qJD(2);
t21 = t112 * t135 + t114 * t136 + t159 * t75 + t160 * t77 + t207 * t187;
t22 = t113 * t135 + t115 * t136 + t159 * t76 + t160 * t78 + t206 * t187;
t64 = t147 * t264 + t148 * t159 + t149 * t160;
t6 = -(t107 * t159 + t108 * t160 + t135 * t148 + t136 * t149) * t191 + (t21 * t186 + (t22 - t260) * t187) * t189 + (t64 * t189 + (t282 + (t48 - t271) * t187) * t191) * qJD(2);
t81 = -t212 * t189 - t271;
t200 = -t191 * ((t260 + (t212 * t191 + t227) * qJD(2)) * t191 + (t18 * t187 + t17 * t186 - (qJD(2) * t147 + (-t148 * t184 + t108) * t183 + (-t149 * t184 - t107) * t182) * t191 + t81 * qJD(2)) * t189) + t6 * t264 + t5 * t267 + (-t191 * t63 + (t186 * t45 + t280) * t189) * t237 + (-t191 * t64 + (t187 * t48 + t282) * t189) * t236 + (t227 * t189 - t191 * t81) * t249;
t13 = t186 * t20 - t187 * t19;
t14 = t186 * t22 - t187 * t21;
t199 = t5 * t288 + t6 * t289 + (-t17 * t187 + t18 * t186) * t286 + t13 * t267 / 0.2e1 + t14 * t264 / 0.2e1 + (t186 * t53 - t187 * t52) * t249 / 0.2e1 + (t186 * (t186 * t46 - t187 * t45) + t187 * (t186 * t48 - t187 * t47)) * t248 / 0.2e1;
t195 = qJD(2) * (-Icges(3,5) * t189 - Icges(3,6) * t191);
t167 = t186 * t195;
t139 = t251 * t187;
t138 = t251 * t186;
t131 = (-Icges(4,1) * t188 - t274) * t246 + (Icges(4,5) * t189 + t224 * t191) * qJD(2);
t130 = (-Icges(4,2) * t190 - t275) * t246 + (Icges(4,6) * t189 + t221 * t191) * qJD(2);
t118 = t297 * t296;
t104 = t117 * t249;
t103 = t116 * t264;
t100 = t254 * t187;
t99 = t254 * t186;
t97 = t193 * t186 - t187 * t243;
t96 = t238 * t187;
t95 = t238 * t186;
t94 = rSges(4,1) * t145 + rSges(4,2) * t144 + rSges(4,3) * t236;
t93 = rSges(4,1) * t143 + rSges(4,2) * t142 + rSges(4,3) * t237;
t92 = Icges(4,1) * t145 + Icges(4,4) * t144 + Icges(4,5) * t236;
t91 = Icges(4,1) * t143 + Icges(4,4) * t142 + Icges(4,5) * t237;
t90 = Icges(4,4) * t145 + Icges(4,2) * t144 + Icges(4,6) * t236;
t89 = Icges(4,4) * t143 + Icges(4,2) * t142 + Icges(4,6) * t237;
t86 = -t126 * t191 - t164 * t264;
t85 = t125 * t191 + t164 * t267;
t84 = -t210 * t189 - t258;
t83 = -t117 * t191 - t150 * t264;
t70 = t161 * t264 + t162 * t175 + t163 * t176;
t69 = t161 * t267 + t162 * t173 + t163 * t174;
t68 = t213 * t189;
t67 = t239 * t187;
t66 = t239 * t186;
t65 = -t117 * t267 + t103;
t62 = t125 * t186 + t126 * t187 + t250;
t59 = t253 * t264 - t232;
t58 = t127 * t191 + t146 * t267 + t82;
t57 = t120 * t264 + t122 * t175 + t124 * t176;
t54 = t119 * t267 + t121 * t173 + t123 * t174;
t51 = t186 * t93 + t187 * t94 - t252;
t50 = -t132 * t264 - t191 * t94 + (t126 * t189 - t164 * t263) * qJD(2);
t49 = t132 * t267 + t191 * t93 + (-t125 * t189 + t164 * t266) * qJD(2);
t44 = t103 + (t127 * t187 - t255 * t186) * t189;
t43 = -t191 * t80 + t104 + (-t109 * t189 - t150 * t248) * t187;
t42 = -t116 * t249 + t242;
t41 = t256 * t186 + t255 * t187 + t250;
t40 = (-t186 * t94 + t187 * t93) * t189 + t213 * t248;
t39 = (-t117 * t248 - t189 * t80) * t186 + t278;
t38 = t283 * t187 + (t79 + t97) * t186 - t252;
t34 = t128 * t249 + t104 - t283 * t191 + (t257 * t189 + t253 * t248) * t187;
t33 = t137 * t267 + t191 * t97 + (t146 * t266 - t256 * t189) * qJD(2) + t242;
t32 = t122 * t144 + t124 * t145 + t175 * t90 + t176 * t92 + t204 * t187;
t31 = t121 * t144 + t123 * t145 + t175 * t89 + t176 * t91 + t205 * t187;
t30 = t122 * t142 + t124 * t143 + t173 * t90 + t174 * t92 + t204 * t186;
t29 = t121 * t142 + t123 * t143 + t173 * t89 + t174 * t91 + t205 * t186;
t27 = (t214 * qJD(2) - t88) * t191 + (qJD(2) * t120 - t188 * t90 + t190 * t92 + (-t122 * t190 - t124 * t188) * qJD(3)) * t189;
t26 = (t215 * qJD(2) - t87) * t191 + (qJD(2) * t119 - t188 * t89 + t190 * t91 + (-t121 * t190 - t123 * t188) * qJD(3)) * t189;
t23 = (t127 * t248 + t189 * t97) * t187 + (-qJD(2) * t232 - t283 * t189) * t186 + t278;
t16 = t186 * t32 - t187 * t31;
t15 = t186 * t30 - t187 * t29;
t9 = -(t130 * t175 + t131 * t176 + t144 * t162 + t145 * t163) * t191 + (t31 * t186 + (t32 - t259) * t187) * t189 + (t70 * t189 + (t281 + (t57 - t258) * t187) * t191) * qJD(2);
t8 = -(t130 * t173 + t131 * t174 + t142 * t162 + t143 * t163) * t191 + (t30 * t187 + (t29 - t259) * t186) * t189 + (t69 * t189 + (t279 + (t54 - t258) * t186) * t191) * qJD(2);
t1 = [0; -m(3) * t118 + m(4) * t51 + m(5) * t38; (t41 * t38 + t66 * t95 + t67 * t96) * t292 + (t100 * t139 + t138 * t99 + t51 * t62) * t293 + 0.2e1 * m(3) * (-t118 + t296) * t297 * (rSges(3,1) * t191 - rSges(3,2) * t189) + (-t290 * t167 - t13 - t15) * t187 + (t14 + t16 + (-t186 * t167 + t297 * t195) * t187) * t186; m(4) * t40 + m(5) * t23; (t186 * t27 - t187 * t26) * t286 + t9 * t289 + t8 * t288 + (t15 * t289 + t16 * t287) * t189 + m(5) * (t23 * t41 + t33 * t96 + t34 * t95 + t44 * t38 + t58 * t67 + t59 * t66) + m(4) * (t100 * t85 + t138 * t50 + t139 * t49 + t40 * t62 + t51 * t68 + t86 * t99) + (t189 * (t186 * t61 - t187 * t60) / 0.2e1 + ((t186 * t57 - t187 * t56) * t287 + (t186 * t55 - t187 * t54) * t289) * t191) * qJD(2) + t199; (-t191 * t70 + (t187 * t57 + t281) * t189) * t236 + t9 * t264 + (-t191 * t69 + (t186 * t54 + t279) * t189) * t237 + t8 * t267 + (t44 * t23 + t33 * t58 + t34 * t59) * t292 + (t226 * t189 - t84 * t191) * t249 - t191 * ((t259 + (t210 * t191 + t226) * qJD(2)) * t191 + (t27 * t187 + t26 * t186 - (qJD(2) * t161 - t130 * t188 + t131 * t190 - t162 * t245 - t163 * t247) * t191 + t84 * qJD(2)) * t189) + (t40 * t68 + t49 * t85 + t50 * t86) * t293 + t200; m(5) * t39; m(5) * (t38 * t65 + t39 * t41 + t42 * t96 + t43 * t95 + t66 * t83 + t67 * t82) + t199; m(5) * (t23 * t65 + t33 * t82 + t34 * t83 + t39 * t44 + t42 * t58 + t43 * t59) + t200; (t39 * t65 + t42 * t82 + t43 * t83) * t292 + t200;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
