% Calculate time derivative of joint inertia matrix for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:27:52
% EndTime: 2022-01-23 09:27:58
% DurationCPUTime: 3.56s
% Computational Cost: add. (4968->270), mult. (3832->372), div. (0->0), fcn. (2806->8), ass. (0->156)
t149 = cos(qJ(4));
t147 = sin(qJ(4));
t225 = Icges(6,4) * t147;
t227 = Icges(5,4) * t147;
t260 = t225 + t227 + (Icges(5,2) + Icges(6,2)) * t149;
t224 = Icges(6,4) * t149;
t226 = Icges(5,4) * t149;
t259 = t224 + t226 + (Icges(5,1) + Icges(6,1)) * t147;
t176 = -Icges(6,2) * t147 + t224;
t177 = -Icges(5,2) * t147 + t226;
t178 = Icges(6,1) * t149 - t225;
t179 = Icges(5,1) * t149 - t227;
t258 = (t176 + t177) * t149 + (t178 + t179) * t147;
t255 = t260 * t147 - t259 * t149;
t145 = qJ(1) + pkin(8);
t142 = qJ(3) + t145;
t138 = cos(t142);
t144 = qJD(1) + qJD(3);
t220 = t138 * t144;
t254 = -rSges(6,1) - pkin(4);
t118 = Icges(6,5) * t147 + Icges(6,6) * t149;
t119 = Icges(5,5) * t147 + Icges(5,6) * t149;
t253 = t118 + t119;
t137 = sin(t142);
t209 = qJD(4) * t149;
t215 = t144 * t147;
t250 = -t137 * t209 - t138 * t215;
t174 = Icges(6,5) * t149 - Icges(6,6) * t147;
t175 = Icges(5,5) * t149 - Icges(5,6) * t147;
t249 = t255 * t144 + (t174 + t175) * qJD(4);
t210 = qJD(4) * t147;
t154 = t258 * qJD(4) + t259 * t209 - t260 * t210;
t128 = t137 * rSges(6,3);
t139 = pkin(4) * t149 + pkin(3);
t146 = -qJ(5) - pkin(7);
t217 = t138 * t149;
t218 = t138 * t147;
t54 = rSges(6,1) * t217 - rSges(6,2) * t218 - t137 * t146 + t138 * t139 + t128;
t166 = t177 * t138;
t68 = Icges(5,6) * t137 + t166;
t168 = t179 * t138;
t72 = Icges(5,5) * t137 + t168;
t182 = t147 * t68 - t149 * t72;
t248 = t137 * t182;
t165 = t176 * t138;
t66 = Icges(6,6) * t137 + t165;
t167 = t178 * t138;
t70 = Icges(6,5) * t137 + t167;
t186 = t147 * t66 - t149 * t70;
t247 = t137 * t186;
t67 = -Icges(5,6) * t138 + t137 * t177;
t71 = -Icges(5,5) * t138 + t137 * t179;
t184 = t147 * t67 - t149 * t71;
t246 = t138 * t184;
t65 = -Icges(6,6) * t138 + t137 * t176;
t69 = -Icges(6,5) * t138 + t137 * t178;
t188 = t147 * t65 - t149 * t69;
t245 = t138 * t188;
t211 = pkin(2) * cos(t145) + cos(qJ(1)) * pkin(1);
t219 = t138 * t146;
t222 = t137 * t147;
t244 = rSges(6,2) * t222 + t138 * rSges(6,3) - t219;
t241 = 2 * m(4);
t240 = 2 * m(5);
t239 = 2 * m(6);
t231 = rSges(5,1) * t149;
t104 = (-rSges(5,2) * t147 + t231) * qJD(4);
t236 = m(5) * t104;
t125 = rSges(5,1) * t147 + rSges(5,2) * t149;
t235 = m(5) * t125;
t134 = t138 * pkin(7);
t221 = t137 * t149;
t233 = t134 + t137 * (-pkin(3) + t139) + rSges(6,1) * t221 - t244;
t212 = -t138 * pkin(3) - t137 * pkin(7);
t232 = t212 + t54;
t230 = rSges(6,1) * t149;
t229 = rSges(6,2) * t149;
t129 = t137 * rSges(5,3);
t206 = t137 * t215;
t228 = rSges(5,2) * t206 + rSges(5,3) * t220;
t223 = t137 * t144;
t216 = t144 * t146;
t213 = rSges(5,2) * t222 + t138 * rSges(5,3);
t204 = t137 * t210;
t208 = -rSges(5,1) * t204 + rSges(5,2) * t250;
t207 = rSges(6,2) * t206 + rSges(6,3) * t220 + qJD(5) * t137;
t202 = t138 * t209;
t199 = -pkin(3) - t231;
t124 = rSges(6,1) * t147 + t229;
t198 = -pkin(4) * t147 - t124;
t197 = -t139 - t230;
t84 = t138 * rSges(4,1) - rSges(4,2) * t137;
t82 = -rSges(4,1) * t220 + rSges(4,2) * t223;
t195 = -pkin(2) * sin(t145) - sin(qJ(1)) * pkin(1);
t83 = -rSges(4,1) * t137 - rSges(4,2) * t138;
t189 = t147 * t69 + t149 * t65;
t187 = t147 * t70 + t149 * t66;
t185 = t147 * t71 + t149 * t67;
t183 = t147 * t72 + t149 * t68;
t181 = t197 * t137;
t76 = rSges(5,1) * t217 - rSges(5,2) * t218 + t129;
t171 = t250 * rSges(6,2) - qJD(5) * t138 - t137 * t216 + t254 * t204;
t170 = t195 * qJD(1);
t169 = t211 * qJD(1);
t81 = t83 * t144;
t164 = t175 * t138;
t163 = t174 * t138;
t56 = t76 - t212;
t161 = (t254 * t147 - t229) * qJD(4);
t55 = t137 * t199 + t134 + t213;
t156 = Icges(5,3) * t144 - qJD(4) * t119;
t155 = Icges(6,3) * t144 - qJD(4) * t118;
t53 = t181 + t244;
t153 = (t249 * t137 + (-t182 - t186) * qJD(4) - t258 * t223) * t137 / 0.2e1 - (-t249 * t138 + (-t184 - t188) * qJD(4)) * t138 / 0.2e1 + (-t137 * t255 - t253 * t138 + t185 + t189) * t223 / 0.2e1 + (t137 * t253 - t138 * t255 + t183 + t187) * t220 / 0.2e1 - ((t165 + t166) * t149 + (t167 + t168) * t147) * t220 / 0.2e1;
t26 = (t199 * t138 + (-rSges(5,3) - pkin(7)) * t137) * t144 - t208;
t22 = (t138 * t197 - t128) * t144 - t171;
t113 = pkin(7) * t220;
t25 = -rSges(5,2) * t202 - pkin(3) * t223 + t113 + (-t138 * t210 - t144 * t221) * rSges(5,1) + t228;
t21 = t144 * t181 + (t161 - t216) * t138 + t207;
t103 = (-rSges(6,2) * t147 + t230) * qJD(4);
t80 = t84 + t211;
t79 = t195 + t83;
t78 = t198 * t138;
t77 = t198 * t137;
t74 = rSges(5,1) * t221 - t213;
t64 = Icges(5,3) * t137 + t164;
t63 = -Icges(5,3) * t138 + t137 * t175;
t62 = Icges(6,3) * t137 + t163;
t61 = -Icges(6,3) * t138 + t137 * t174;
t60 = -t169 + t82;
t59 = t81 + t170;
t52 = t56 + t211;
t51 = t195 + t55;
t50 = t54 + t211;
t49 = t195 + t53;
t40 = t137 * t156 + t144 * t164;
t39 = t138 * t156 - t175 * t223;
t38 = t137 * t155 + t144 * t163;
t37 = t138 * t155 - t174 * t223;
t36 = pkin(4) * t250 - t103 * t137 - t124 * t220;
t35 = t124 * t223 - t103 * t138 + (-t202 + t206) * pkin(4);
t24 = -t169 + t26;
t23 = t170 + t25;
t20 = t137 * t64 - t182 * t138;
t19 = t137 * t63 - t246;
t18 = t137 * t62 - t186 * t138;
t17 = t137 * t61 - t245;
t16 = -t138 * t64 - t248;
t15 = -t184 * t137 - t138 * t63;
t14 = -t138 * t62 - t247;
t13 = -t188 * t137 - t138 * t61;
t12 = -t169 + t22;
t11 = t170 + t21;
t2 = ((-t76 + t129) * t144 + t208) * t137 + (-qJD(4) * t125 * t138 + t144 * t74 + t228) * t138;
t1 = (((rSges(6,3) - pkin(7)) * t137 - t232) * t144 + t171) * t137 + (-t113 + (-t219 + t233) * t144 + t138 * t161 + t207) * t138;
t3 = [(t11 * t50 + t12 * t49) * t239 + (t23 * t52 + t24 * t51) * t240 + (t59 * t80 + t60 * t79) * t241 + t154; 0; 0; m(6) * (t11 * t54 + t12 * t53 + t21 * t50 + t22 * t49) + m(5) * (t23 * t56 + t24 * t55 + t25 * t52 + t26 * t51) + m(4) * (t59 * t84 + t60 * t83 + t79 * t82 + t80 * t81) + t154; 0; (t21 * t54 + t22 * t53) * t239 + (t25 * t56 + t26 * t55) * t240 + (t81 * t84 + t82 * t83) * t241 + t154; ((-t144 * t52 - t24) * t138 + (t144 * t51 - t23) * t137) * t235 + (-t137 * t52 - t138 * t51) * t236 + t153 + m(6) * (t11 * t77 + t12 * t78 + t35 * t49 + t36 * t50); m(5) * t2 + m(6) * t1; m(6) * (t21 * t77 + t22 * t78 + t35 * t53 + t36 * t54) + ((-t144 * t56 - t26) * t138 + (t144 * t55 - t25) * t137) * t235 + (-t137 * t56 - t138 * t55) * t236 + t153; ((t137 * t74 + t138 * t76) * t2 + (t137 ^ 2 + t138 ^ 2) * t125 * t104) * t240 + t137 * ((t137 * t39 + (t19 + t248) * t144) * t137 + (t20 * t144 + (t209 * t67 + t210 * t71) * t138 + (-t183 * qJD(4) - t144 * t184 - t40) * t137) * t138) - t138 * ((t138 * t40 + (t16 + t246) * t144) * t138 + (t15 * t144 + (-t209 * t68 - t210 * t72) * t137 + (t185 * qJD(4) - t182 * t144 - t39) * t138) * t137) + (t78 * t35 + t77 * t36 + (t137 * t233 + t138 * t232) * t1) * t239 + t137 * ((t137 * t37 + (t17 + t247) * t144) * t137 + (t18 * t144 + (t209 * t65 + t210 * t69) * t138 + (-t187 * qJD(4) - t144 * t188 - t38) * t137) * t138) - t138 * ((t138 * t38 + (t14 + t245) * t144) * t138 + (t13 * t144 + (-t209 * t66 - t210 * t70) * t137 + (t189 * qJD(4) - t144 * t186 - t37) * t138) * t137) + ((-t13 - t15) * t138 + (t14 + t16) * t137) * t223 + ((-t17 - t19) * t138 + (t18 + t20) * t137) * t220; m(6) * ((t144 * t49 - t11) * t138 + (t144 * t50 + t12) * t137); 0; m(6) * ((t144 * t53 - t21) * t138 + (t144 * t54 + t22) * t137); m(6) * ((t144 * t78 - t36) * t138 + (t144 * t77 + t35) * t137); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
