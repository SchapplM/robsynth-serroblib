% Calculate time derivative of joint inertia matrix for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR11_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR11_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR11_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:33
% EndTime: 2019-12-31 18:05:40
% DurationCPUTime: 4.37s
% Computational Cost: add. (4834->450), mult. (13001->677), div. (0->0), fcn. (12281->6), ass. (0->220)
t150 = cos(qJ(4));
t147 = sin(qJ(4));
t151 = cos(qJ(1));
t208 = qJD(4) * t151;
t197 = t147 * t208;
t148 = sin(qJ(1));
t213 = qJD(1) * t148;
t155 = t150 * t213 + t197;
t149 = cos(qJ(5));
t146 = sin(qJ(5));
t235 = Icges(6,4) * t146;
t171 = Icges(6,1) * t149 - t235;
t100 = Icges(6,5) * t147 + t150 * t171;
t234 = Icges(6,4) * t149;
t169 = -Icges(6,2) * t146 + t234;
t97 = Icges(6,6) * t147 + t150 * t169;
t274 = t100 * t149 - t146 * t97;
t263 = -t148 / 0.2e1;
t196 = t150 * t208;
t193 = -qJD(5) * t147 - qJD(1);
t165 = t193 * t151;
t192 = qJD(1) * t147 + qJD(5);
t267 = t148 * t192 - t196;
t59 = t146 * t267 + t149 * t165;
t60 = t146 * t165 - t149 * t267;
t37 = t60 * rSges(6,1) + t59 * rSges(6,2) + t155 * rSges(6,3);
t273 = pkin(4) * t196 + t155 * pkin(7) + t37;
t168 = Icges(5,5) * t147 + Icges(5,6) * t150;
t95 = Icges(5,3) * t151 + t148 * t168;
t272 = qJD(1) * t95;
t144 = t148 ^ 2;
t145 = t151 ^ 2;
t271 = t144 + t145;
t270 = -t100 * t146 - t149 * t97;
t142 = t151 * qJ(2);
t188 = rSges(5,1) * t147 + rSges(5,2) * t150;
t252 = -pkin(1) - qJ(3);
t161 = -t188 + t252;
t259 = -rSges(5,3) - pkin(6);
t152 = t148 * t161 + t151 * t259;
t84 = t142 + t152;
t140 = t148 * qJ(2);
t214 = t151 * pkin(1) + t140;
t200 = t151 * qJ(3) + t214;
t220 = t150 * t151;
t225 = t147 * t151;
t217 = rSges(5,1) * t225 + rSges(5,2) * t220;
t85 = t148 * t259 + t200 + t217;
t269 = t148 * t85 + t151 * t84;
t127 = rSges(5,1) * t196;
t212 = qJD(1) * t151;
t216 = qJ(2) * t212 + qJD(2) * t148;
t201 = qJD(3) * t151 + t216;
t52 = -rSges(5,2) * t197 + qJD(1) * t152 + t127 + t201;
t250 = rSges(5,2) * t147;
t124 = rSges(5,1) * t150 - t250;
t139 = qJD(2) * t151;
t215 = pkin(6) * t213 + t139;
t53 = (-qJD(4) * t124 - qJD(3)) * t148 + ((rSges(5,3) - qJ(2)) * t148 + t161 * t151) * qJD(1) + t215;
t268 = t148 * t52 + t151 * t53;
t237 = Icges(5,4) * t147;
t170 = Icges(5,2) * t150 + t237;
t98 = Icges(5,6) * t151 + t148 * t170;
t236 = Icges(5,4) * t150;
t172 = Icges(5,1) * t147 + t236;
t101 = Icges(5,5) * t151 + t148 * t172;
t266 = 2 * m(5);
t265 = 2 * m(6);
t264 = t147 / 0.2e1;
t262 = -t151 / 0.2e1;
t260 = rSges(3,2) - pkin(1);
t258 = rSges(6,3) + pkin(7);
t257 = m(5) * t124;
t256 = pkin(4) * t147;
t255 = pkin(4) * t150;
t254 = pkin(6) * t151;
t253 = -qJD(1) / 0.2e1;
t209 = qJD(4) * t150;
t211 = qJD(4) * t147;
t167 = Icges(6,5) * t149 - Icges(6,6) * t146;
t207 = qJD(5) * t150;
t72 = (-Icges(6,5) * t146 - Icges(6,6) * t149) * t207 + (Icges(6,3) * t150 - t147 * t167) * qJD(4);
t78 = (-Icges(6,1) * t146 - t234) * t207 + (Icges(6,5) * t150 - t147 * t171) * qJD(4);
t94 = Icges(6,3) * t147 + t150 * t167;
t153 = t150 * t149 * t78 + t147 * t72 + t94 * t209 - t274 * t211;
t75 = (-Icges(6,2) * t149 - t235) * t207 + (Icges(6,6) * t150 - t147 * t169) * qJD(4);
t249 = t146 * t75;
t47 = t147 * t94 + t274 * t150;
t251 = ((qJD(5) * t270 - t249) * t150 + t153) * t147 + t47 * t209;
t247 = t147 * t98;
t99 = -Icges(5,6) * t148 + t151 * t170;
t246 = t147 * t99;
t189 = -pkin(7) * t150 + t256;
t221 = t149 * t151;
t224 = t148 * t146;
t110 = -t147 * t224 + t221;
t223 = t148 * t149;
t111 = t146 * t151 + t147 * t223;
t187 = -t111 * rSges(6,1) - t110 * rSges(6,2);
t222 = t148 * t150;
t70 = -rSges(6,3) * t222 - t187;
t240 = -t189 * t148 - t70;
t134 = pkin(4) * t225;
t205 = pkin(7) * t220;
t112 = -t146 * t225 - t223;
t113 = t147 * t221 - t224;
t218 = t113 * rSges(6,1) + t112 * rSges(6,2);
t71 = -rSges(6,3) * t220 + t218;
t239 = -t134 + t205 - t71;
t186 = rSges(6,1) * t149 - rSges(6,2) * t146;
t81 = (-rSges(6,1) * t146 - rSges(6,2) * t149) * t207 + (rSges(6,3) * t150 - t147 * t186) * qJD(4);
t238 = -t189 * qJD(4) + t81;
t96 = -Icges(5,3) * t148 + t151 * t168;
t230 = qJD(1) * t96;
t227 = t101 * t150;
t102 = -Icges(5,5) * t148 + t151 * t172;
t226 = t102 * t150;
t104 = rSges(6,3) * t147 + t150 * t186;
t125 = pkin(7) * t147 + t255;
t219 = t104 + t125;
t210 = qJD(4) * t148;
t206 = -rSges(4,3) + t252;
t66 = Icges(6,4) * t111 + Icges(6,2) * t110 - Icges(6,6) * t222;
t68 = Icges(6,1) * t111 + Icges(6,4) * t110 - Icges(6,5) * t222;
t185 = t146 * t66 - t149 * t68;
t64 = Icges(6,5) * t111 + Icges(6,6) * t110 - Icges(6,3) * t222;
t28 = t147 * t64 - t150 * t185;
t39 = t100 * t111 + t110 * t97 - t222 * t94;
t204 = t28 / 0.2e1 + t39 / 0.2e1;
t67 = Icges(6,4) * t113 + Icges(6,2) * t112 - Icges(6,6) * t220;
t69 = Icges(6,1) * t113 + Icges(6,4) * t112 - Icges(6,5) * t220;
t184 = t146 * t67 - t149 * t69;
t65 = Icges(6,5) * t113 + Icges(6,6) * t112 - Icges(6,3) * t220;
t29 = t147 * t65 - t150 * t184;
t40 = t113 * t100 + t112 * t97 - t220 * t94;
t203 = t40 / 0.2e1 + t29 / 0.2e1;
t198 = t150 * t212;
t195 = t211 / 0.2e1;
t87 = t219 * t151;
t119 = t188 * qJD(4);
t194 = t271 * t119;
t191 = t252 - t256;
t61 = t193 * t223 + (-t148 * t209 - t151 * t192) * t146;
t62 = t192 * t221 + (t146 * t193 + t149 * t209) * t148;
t190 = t62 * rSges(6,1) + t61 * rSges(6,2);
t24 = t110 * t66 + t111 * t68 - t222 * t64;
t25 = t110 * t67 + t111 * t69 - t222 * t65;
t181 = t25 * t148 - t151 * t24;
t180 = t148 * t24 + t151 * t25;
t26 = t112 * t66 + t113 * t68 - t220 * t64;
t27 = t112 * t67 + t113 * t69 - t220 * t65;
t179 = t27 * t148 - t151 * t26;
t178 = t148 * t26 + t151 * t27;
t177 = t29 * t148 - t151 * t28;
t176 = t148 * t28 + t151 * t29;
t175 = t148 * t71 - t151 * t70;
t174 = t101 * t147 + t150 * t98;
t173 = t102 * t147 + t150 * t99;
t164 = rSges(3,3) * t151 + t148 * t260;
t163 = t174 * t151;
t162 = t173 * t148;
t160 = qJD(4) * (Icges(5,1) * t150 - t237);
t159 = qJD(4) * (-Icges(5,2) * t147 + t236);
t158 = qJD(4) * (Icges(5,5) * t150 - Icges(5,6) * t147);
t157 = rSges(4,2) * t151 + t148 * t206;
t156 = t147 * t210 - t198;
t154 = t150 * t258 + t191;
t107 = -rSges(3,2) * t151 + t148 * rSges(3,3) + t214;
t106 = t142 + t164;
t105 = -rSges(5,3) * t148 + t217;
t103 = t151 * rSges(5,3) + t148 * t188;
t93 = t148 * rSges(4,2) + rSges(4,3) * t151 + t200;
t92 = t142 + t157;
t89 = t139 + (t260 * t151 + (-rSges(3,3) - qJ(2)) * t148) * qJD(1);
t88 = qJD(1) * t164 + t216;
t86 = t219 * t148;
t83 = -qJD(3) * t148 + t139 + ((-rSges(4,2) - qJ(2)) * t148 + t206 * t151) * qJD(1);
t82 = qJD(1) * t157 + t201;
t74 = t148 * t158 + t230;
t73 = t151 * t158 - t272;
t55 = t104 * t220 + t147 * t71;
t54 = -t104 * t222 - t147 * t70;
t51 = -t148 * pkin(6) - t220 * t258 + t134 + t200 + t218;
t50 = t148 * t154 + t142 + t187 - t254;
t49 = -t148 * t96 + t173 * t151;
t48 = -t148 * t95 + t163;
t46 = t151 * t96 + t162;
t45 = t148 * t174 + t151 * t95;
t43 = t175 * t150;
t42 = qJD(1) * t87 + t148 * t238;
t41 = t151 * t238 - t213 * t219;
t38 = rSges(6,3) * t156 + t190;
t36 = Icges(6,1) * t62 + Icges(6,4) * t61 + Icges(6,5) * t156;
t35 = Icges(6,1) * t60 + Icges(6,4) * t59 + Icges(6,5) * t155;
t34 = Icges(6,4) * t62 + Icges(6,2) * t61 + Icges(6,6) * t156;
t33 = Icges(6,4) * t60 + Icges(6,2) * t59 + Icges(6,6) * t155;
t32 = Icges(6,5) * t62 + Icges(6,6) * t61 + Icges(6,3) * t156;
t31 = Icges(6,5) * t60 + Icges(6,6) * t59 + Icges(6,3) * t155;
t30 = t148 * t240 + t151 * t239;
t23 = (-qJD(3) + (-t147 * t258 - t255) * qJD(4)) * t148 + (t151 * t154 - t140) * qJD(1) - t190 + t215;
t22 = (t148 * t191 - t254) * qJD(1) + t201 + t273;
t21 = (t104 * t210 - t38) * t147 + (-qJD(4) * t70 - t104 * t212 - t148 * t81) * t150;
t20 = (-t104 * t208 + t37) * t147 + (qJD(4) * t71 - t104 * t213 + t151 * t81) * t150;
t16 = -t94 * t198 + t62 * t100 + t110 * t75 + t111 * t78 + t61 * t97 + (-t150 * t72 + t211 * t94) * t148;
t15 = t94 * t197 + t60 * t100 + t112 * t75 + t113 * t78 + t59 * t97 + (-t151 * t72 + t213 * t94) * t150;
t14 = -t175 * t211 + (t148 * t37 - t151 * t38 + (t148 * t70 + t151 * t71) * qJD(1)) * t150;
t13 = (qJD(1) * t240 - t273) * t151 + (-t38 - t125 * t210 + (t205 - t239) * qJD(1)) * t148;
t12 = t40 * t147 - t150 * t178;
t11 = t39 * t147 - t150 * t180;
t10 = (qJD(4) * t184 + t31) * t147 + (qJD(4) * t65 - t146 * t33 + t149 * t35 + (-t146 * t69 - t149 * t67) * qJD(5)) * t150;
t9 = (qJD(4) * t185 + t32) * t147 + (qJD(4) * t64 - t146 * t34 + t149 * t36 + (-t146 * t68 - t149 * t66) * qJD(5)) * t150;
t8 = -t65 * t198 + t110 * t33 + t111 * t35 + t61 * t67 + t62 * t69 + (-t150 * t31 + t211 * t65) * t148;
t7 = -t64 * t198 + t110 * t34 + t111 * t36 + t61 * t66 + t62 * t68 + (-t150 * t32 + t211 * t64) * t148;
t6 = t65 * t197 + t112 * t33 + t113 * t35 + t59 * t67 + t60 * t69 + (-t151 * t31 + t213 * t65) * t150;
t5 = t64 * t197 + t112 * t34 + t113 * t36 + t59 * t66 + t60 * t68 + (-t151 * t32 + t213 * t64) * t150;
t4 = -qJD(1) * t180 - t8 * t148 + t151 * t7;
t3 = -qJD(1) * t178 - t6 * t148 + t151 * t5;
t2 = (qJD(4) * t180 + t16) * t147 + (qJD(1) * t181 + qJD(4) * t39 - t148 * t7 - t151 * t8) * t150;
t1 = (qJD(4) * t178 + t15) * t147 + (qJD(1) * t179 + qJD(4) * t40 - t148 * t5 - t151 * t6) * t150;
t17 = [(t22 * t51 + t23 * t50) * t265 - t147 * t160 - t172 * t209 + t170 * t211 + (t52 * t85 + t53 * t84) * t266 + 0.2e1 * m(3) * (t106 * t89 + t107 * t88) + 0.2e1 * m(4) * (t82 * t93 + t83 * t92) + t153 + t270 * t207 + (-t249 - t159) * t150; m(6) * (t148 * t23 - t151 * t22 + (t148 * t51 + t151 * t50) * qJD(1)) + m(5) * (qJD(1) * t269 + t148 * t53 - t151 * t52) + m(3) * (t148 * t89 - t151 * t88 + (t106 * t151 + t107 * t148) * qJD(1)) + m(4) * (t148 * t83 - t151 * t82 + (t148 * t93 + t151 * t92) * qJD(1)); 0; m(6) * (t148 * t22 + t151 * t23 + (-t148 * t50 + t151 * t51) * qJD(1)) + m(5) * ((-t148 * t84 + t151 * t85) * qJD(1) + t268) + m(4) * (t148 * t82 + t151 * t83 + (-t148 * t92 + t151 * t93) * qJD(1)); 0; 0; m(6) * (t22 * t86 + t23 * t87 + t41 * t50 + t42 * t51) + m(5) * (-t119 * t269 + t124 * t268) - (t144 / 0.2e1 + t145 / 0.2e1) * t168 * qJD(4) + ((-t226 / 0.2e1 + t246 / 0.2e1 + t85 * t257 - t203) * t151 + (-t227 / 0.2e1 + t247 / 0.2e1 - t84 * t257 - t204) * t148) * qJD(1) + (-qJD(4) * t173 - t147 * (-qJD(1) * t98 + t151 * t159) + t150 * (-qJD(1) * t101 + t151 * t160) + t10 + t15) * t263 + (-qJD(4) * t174 - t147 * (qJD(1) * t99 + t148 * t159) + t150 * (qJD(1) * t102 + t148 * t160) + t16 + t9) * t151 / 0.2e1; m(6) * (t41 * t148 - t151 * t42 + (t148 * t86 + t151 * t87) * qJD(1)); m(6) * (t42 * t148 + t41 * t151 + (-t148 * t87 + t151 * t86) * qJD(1)) - m(5) * t194; (-t124 * t194 + (-t151 * t127 + (-t124 * t144 + t145 * t250) * qJD(4) + (rSges(5,3) * t271 - t151 * t103 + t148 * t105) * qJD(1)) * (-t148 * t103 - t105 * t151)) * t266 - t148 * ((t148 * t73 + (-t48 + t162) * qJD(1)) * t148 + (-t49 * qJD(1) + (t101 * t209 - t211 * t98 - t272) * t151 + (-t74 + (-t226 + t246) * qJD(4) + (-t174 + t96) * qJD(1)) * t148) * t151) + t151 * ((t151 * t74 + (-t46 + t163) * qJD(1)) * t151 + (-t45 * qJD(1) + (-t102 * t209 + t211 * t99 + t230) * t148 + (-t73 + (t227 - t247) * qJD(4) + (-t173 - t95) * qJD(1)) * t151) * t148) + (t30 * t13 + t41 * t87 + t42 * t86) * t265 - t148 * t3 + t151 * t4 + (t46 * t148 - t151 * t45 + t181) * t213 + (t49 * t148 - t151 * t48 + t179) * t212; m(6) * (t20 * t51 + t21 * t50 + t22 * t55 + t23 * t54) + (t148 * t204 + t151 * t203) * t211 + ((-t15 / 0.2e1 - t10 / 0.2e1) * t151 + (-t9 / 0.2e1 - t16 / 0.2e1) * t148 + (t148 * t203 - t151 * t204) * qJD(1)) * t150 + t251; m(6) * (t21 * t148 - t151 * t20 + (t148 * t55 + t151 * t54) * qJD(1)); m(6) * (t20 * t148 + t151 * t21 + (-t148 * t54 + t151 * t55) * qJD(1)); m(6) * (t43 * t13 + t14 * t30 + t20 * t86 + t21 * t87 + t54 * t41 + t55 * t42) + (t12 * t253 - t179 * t195 + t2 / 0.2e1 + (-qJD(1) * t29 + t9) * t264) * t151 + (-t1 / 0.2e1 + t11 * t253 - t181 * t195 + (-qJD(1) * t28 - t10) * t264) * t148 + (t3 * t262 + t4 * t263 - qJD(4) * t177 / 0.2e1 + (t179 * t263 - t181 * t262) * qJD(1)) * t150; (t14 * t43 + t20 * t55 + t21 * t54) * t265 + ((t148 * t11 + t151 * t12 + t147 * t176) * qJD(4) + t251) * t147 + (-t151 * t1 - t148 * t2 + t147 * (-t10 * t151 - t9 * t148) + (t47 * t147 - t150 * t176) * qJD(4) + (-t151 * t11 + t148 * t12 + t147 * t177) * qJD(1)) * t150;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t17(1), t17(2), t17(4), t17(7), t17(11); t17(2), t17(3), t17(5), t17(8), t17(12); t17(4), t17(5), t17(6), t17(9), t17(13); t17(7), t17(8), t17(9), t17(10), t17(14); t17(11), t17(12), t17(13), t17(14), t17(15);];
Mq = res;
