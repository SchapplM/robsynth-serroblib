% Calculate joint inertia matrix for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:37:03
% EndTime: 2019-03-09 08:37:10
% DurationCPUTime: 3.62s
% Computational Cost: add. (4650->431), mult. (11547->618), div. (0->0), fcn. (13770->8), ass. (0->188)
t263 = Icges(4,1) + Icges(5,1);
t262 = -Icges(4,4) + Icges(5,5);
t261 = Icges(5,4) + Icges(4,5);
t260 = Icges(4,2) + Icges(5,3);
t259 = Icges(4,6) - Icges(5,6);
t256 = rSges(7,1) + pkin(5);
t253 = rSges(7,3) + qJ(6);
t173 = sin(pkin(9));
t174 = cos(pkin(9));
t179 = cos(qJ(1));
t177 = sin(qJ(1));
t178 = cos(qJ(2));
t218 = t177 * t178;
t144 = t173 * t218 + t174 * t179;
t145 = -t173 * t179 + t174 * t218;
t175 = sin(qJ(5));
t234 = cos(qJ(5));
t106 = -t144 * t234 + t145 * t175;
t107 = t144 * t175 + t145 * t234;
t258 = -t106 * t253 - t107 * t256;
t237 = t177 / 0.2e1;
t257 = t179 / 0.2e1;
t176 = sin(qJ(2));
t220 = t176 * t177;
t247 = t144 * t260 + t145 * t262 - t220 * t259;
t87 = Icges(4,5) * t145 - Icges(4,6) * t144 + Icges(4,3) * t220;
t89 = Icges(5,4) * t145 + Icges(5,2) * t220 + Icges(5,6) * t144;
t255 = t87 + t89;
t217 = t178 * t179;
t146 = t173 * t217 - t174 * t177;
t147 = t173 * t177 + t174 * t217;
t219 = t176 * t179;
t88 = Icges(4,5) * t147 - Icges(4,6) * t146 + Icges(4,3) * t219;
t90 = Icges(5,4) * t147 + Icges(5,2) * t219 + Icges(5,6) * t146;
t254 = t88 + t90;
t246 = t146 * t260 + t147 * t262 - t219 * t259;
t245 = t144 * t262 + t145 * t263 + t220 * t261;
t244 = t146 * t262 + t147 * t263 + t219 * t261;
t221 = t174 * t176;
t222 = t173 * t176;
t132 = t175 * t221 - t222 * t234;
t133 = (t173 * t175 + t174 * t234) * t176;
t76 = Icges(7,5) * t133 + Icges(7,6) * t178 + Icges(7,3) * t132;
t77 = Icges(6,5) * t133 - Icges(6,6) * t132 + Icges(6,3) * t178;
t78 = Icges(7,4) * t133 + Icges(7,2) * t178 + Icges(7,6) * t132;
t79 = Icges(6,4) * t133 - Icges(6,2) * t132 + Icges(6,6) * t178;
t80 = Icges(7,1) * t133 + Icges(7,4) * t178 + Icges(7,5) * t132;
t81 = Icges(6,1) * t133 - Icges(6,4) * t132 + Icges(6,5) * t178;
t252 = (t77 + t78) * t178 + (t80 + t81) * t133 + (t76 - t79) * t132;
t47 = Icges(7,5) * t107 - Icges(7,6) * t220 + Icges(7,3) * t106;
t51 = Icges(7,4) * t107 - Icges(7,2) * t220 + Icges(7,6) * t106;
t55 = Icges(7,1) * t107 - Icges(7,4) * t220 + Icges(7,5) * t106;
t12 = t106 * t47 + t107 * t55 - t220 * t51;
t108 = -t146 * t234 + t147 * t175;
t109 = t146 * t175 + t147 * t234;
t48 = Icges(7,5) * t109 - Icges(7,6) * t219 + Icges(7,3) * t108;
t52 = Icges(7,4) * t109 - Icges(7,2) * t219 + Icges(7,6) * t108;
t56 = Icges(7,1) * t109 - Icges(7,4) * t219 + Icges(7,5) * t108;
t13 = t106 * t48 + t107 * t56 - t220 * t52;
t49 = Icges(6,5) * t107 - Icges(6,6) * t106 - Icges(6,3) * t220;
t53 = Icges(6,4) * t107 - Icges(6,2) * t106 - Icges(6,6) * t220;
t57 = Icges(6,1) * t107 - Icges(6,4) * t106 - Icges(6,5) * t220;
t14 = -t106 * t53 + t107 * t57 - t220 * t49;
t50 = Icges(6,5) * t109 - Icges(6,6) * t108 - Icges(6,3) * t219;
t54 = Icges(6,4) * t109 - Icges(6,2) * t108 - Icges(6,6) * t219;
t58 = Icges(6,1) * t109 - Icges(6,4) * t108 - Icges(6,5) * t219;
t15 = -t106 * t54 + t107 * t58 - t220 * t50;
t28 = t106 * t76 + t107 * t80 - t220 * t78;
t29 = -t106 * t79 + t107 * t81 - t220 * t77;
t251 = (t29 + t28) * t178 + ((-t13 - t15) * t179 + (-t12 - t14) * t177) * t176;
t16 = t108 * t47 + t109 * t55 - t219 * t51;
t17 = t108 * t48 + t109 * t56 - t219 * t52;
t18 = -t108 * t53 + t109 * t57 - t219 * t49;
t19 = -t108 * t54 + t109 * t58 - t219 * t50;
t30 = t108 * t76 + t109 * t80 - t219 * t78;
t31 = -t108 * t79 + t109 * t81 - t219 * t77;
t250 = (t30 + t31) * t178 + ((-t17 - t19) * t179 + (-t16 - t18) * t177) * t176;
t20 = t132 * t47 + t133 * t55 + t178 * t51;
t22 = -t132 * t53 + t133 * t57 + t178 * t49;
t249 = -t20 - t22;
t21 = t132 * t48 + t133 * t56 + t178 * t52;
t23 = -t132 * t54 + t133 * t58 + t178 * t50;
t248 = t21 + t23;
t243 = t108 * t253 + t109 * t256;
t171 = t177 ^ 2;
t172 = t179 ^ 2;
t242 = 0.2e1 * t176;
t241 = m(4) / 0.2e1;
t240 = m(5) / 0.2e1;
t239 = m(6) / 0.2e1;
t238 = m(7) / 0.2e1;
t235 = -t179 / 0.2e1;
t233 = pkin(2) * t178;
t232 = t252 * t178;
t231 = -rSges(7,2) * t220 - t258;
t230 = -rSges(7,2) * t219 + t243;
t229 = rSges(7,2) * t178 + t132 * t253 + t133 * t256;
t228 = rSges(5,3) * t144;
t227 = t179 * rSges(3,3);
t225 = rSges(6,1) * t109 - rSges(6,2) * t108;
t224 = Icges(3,4) * t176;
t223 = Icges(3,4) * t178;
t154 = pkin(2) * t176 - qJ(3) * t178;
t216 = rSges(4,3) * t178 - (rSges(4,1) * t174 - rSges(4,2) * t173) * t176 - t154;
t213 = pkin(2) * t217 + qJ(3) * t219;
t215 = t171 * (qJ(3) * t176 + t233) + t179 * t213;
t214 = -(pkin(3) * t174 + qJ(4) * t173) * t176 - t154;
t212 = pkin(1) * t179 + pkin(7) * t177;
t136 = t144 * qJ(4);
t168 = t179 * pkin(7);
t211 = t168 - t136;
t210 = t171 + t172;
t207 = rSges(5,2) * t178 - (rSges(5,1) * t174 + rSges(5,3) * t173) * t176 + t214;
t206 = rSges(5,1) * t147 + rSges(5,2) * t219 + rSges(5,3) * t146;
t205 = rSges(4,1) * t147 - rSges(4,2) * t146 + rSges(4,3) * t219;
t204 = -pkin(4) * t221 - pkin(8) * t178 + t214;
t203 = -pkin(1) - t233;
t116 = -Icges(5,6) * t178 + (Icges(5,5) * t174 + Icges(5,3) * t173) * t176;
t119 = -Icges(4,6) * t178 + (Icges(4,4) * t174 - Icges(4,2) * t173) * t176;
t202 = t119 / 0.2e1 - t116 / 0.2e1;
t120 = -Icges(5,4) * t178 + (Icges(5,1) * t174 + Icges(5,5) * t173) * t176;
t121 = -Icges(4,5) * t178 + (Icges(4,1) * t174 - Icges(4,4) * t173) * t176;
t201 = t121 / 0.2e1 + t120 / 0.2e1;
t200 = pkin(3) * t147 + t146 * qJ(4);
t199 = t240 + t239 + t238;
t83 = rSges(6,1) * t133 - rSges(6,2) * t132 + rSges(6,3) * t178;
t198 = -t83 + t204;
t197 = t177 * (pkin(3) * t145 + t136) + t179 * t200 + t215;
t196 = t212 + t213;
t195 = t204 - t229;
t194 = rSges(3,1) * t178 - rSges(3,2) * t176;
t193 = -rSges(4,1) * t145 + rSges(4,2) * t144;
t192 = -rSges(6,1) * t107 + rSges(6,2) * t106;
t191 = Icges(3,1) * t178 - t224;
t190 = -Icges(3,2) * t176 + t223;
t189 = Icges(3,5) * t178 - Icges(3,6) * t176;
t186 = rSges(3,1) * t217 - rSges(3,2) * t219 + rSges(3,3) * t177;
t185 = -t21 / 0.2e1 - t31 / 0.2e1 - t30 / 0.2e1 - t23 / 0.2e1;
t184 = -t22 / 0.2e1 - t20 / 0.2e1 - t29 / 0.2e1 - t28 / 0.2e1;
t141 = t147 * pkin(4);
t163 = pkin(8) * t220;
t183 = t177 * (pkin(4) * t145 - t163) + t179 * (-pkin(8) * t219 + t141) + t197;
t182 = t163 + (-pkin(3) - pkin(4)) * t145 + t211;
t181 = t196 + t200;
t180 = t141 + t181;
t170 = t176 ^ 2;
t157 = rSges(2,1) * t179 - rSges(2,2) * t177;
t156 = -rSges(2,1) * t177 - rSges(2,2) * t179;
t155 = rSges(3,1) * t176 + rSges(3,2) * t178;
t151 = Icges(3,5) * t176 + Icges(3,6) * t178;
t125 = Icges(3,3) * t177 + t179 * t189;
t124 = -Icges(3,3) * t179 + t177 * t189;
t115 = t186 + t212;
t114 = t227 + t168 + (-pkin(1) - t194) * t177;
t111 = t216 * t179;
t110 = t216 * t177;
t84 = t179 * t186 + (t177 * t194 - t227) * t177;
t72 = t207 * t179;
t71 = t207 * t177;
t66 = t196 + t205;
t65 = t168 + ((-rSges(4,3) - qJ(3)) * t176 + t203) * t177 + t193;
t62 = -rSges(6,3) * t219 + t225;
t60 = -rSges(6,3) * t220 - t192;
t46 = t198 * t179;
t45 = t198 * t177;
t44 = t181 + t206;
t43 = -t228 + (-rSges(5,1) - pkin(3)) * t145 + ((-rSges(5,2) - qJ(3)) * t176 + t203) * t177 + t211;
t42 = t177 * (rSges(4,3) * t220 - t193) + t179 * t205 + t215;
t41 = t195 * t179;
t40 = t195 * t177;
t39 = t178 * t62 + t219 * t83;
t38 = -t178 * t60 - t220 * t83;
t37 = (-rSges(6,3) - pkin(8)) * t219 + t180 + t225;
t36 = ((rSges(6,3) - qJ(3)) * t176 + t203) * t177 + t182 + t192;
t33 = (t177 * t62 - t179 * t60) * t176;
t32 = t177 * (rSges(5,1) * t145 + rSges(5,2) * t220 + t228) + t179 * t206 + t197;
t27 = (-rSges(7,2) - pkin(8)) * t219 + t180 + t243;
t26 = ((rSges(7,2) - qJ(3)) * t176 + t203) * t177 + t182 + t258;
t25 = t178 * t230 + t219 * t229;
t24 = -t178 * t231 - t220 * t229;
t11 = t177 * t60 + t179 * t62 + t183;
t10 = (t177 * t230 - t179 * t231) * t176;
t9 = t177 * t231 + t179 * t230 + t183;
t8 = t177 * t19 - t179 * t18;
t7 = -t16 * t179 + t17 * t177;
t6 = -t14 * t179 + t15 * t177;
t5 = -t12 * t179 + t13 * t177;
t1 = [Icges(2,3) + (t224 + (t173 * t259 - t174 * t261) * t176 + (Icges(4,3) + Icges(5,2) + Icges(3,2)) * t178) * t178 + (Icges(3,1) * t176 + t223 + (t120 + t121) * t174 + (t116 - t119) * t173) * t176 + m(7) * (t26 ^ 2 + t27 ^ 2) + m(6) * (t36 ^ 2 + t37 ^ 2) + m(4) * (t65 ^ 2 + t66 ^ 2) + m(5) * (t43 ^ 2 + t44 ^ 2) + m(3) * (t114 ^ 2 + t115 ^ 2) + m(2) * (t156 ^ 2 + t157 ^ 2) + t252; (t202 * t144 - t201 * t145 + t151 * t257 + t184) * t179 + (-t202 * t146 + t201 * t147 + t151 * t237 - t185) * t177 + m(7) * (t26 * t41 + t27 * t40) + m(6) * (t36 * t46 + t37 * t45) + m(5) * (t43 * t72 + t44 * t71) + m(4) * (t110 * t66 + t111 * t65) + m(3) * (-t114 * t179 - t115 * t177) * t155 + ((Icges(3,6) * t257 - t177 * t190 / 0.2e1 + t89 / 0.2e1 + t87 / 0.2e1) * t179 + (Icges(3,6) * t237 + t190 * t257 - t90 / 0.2e1 - t88 / 0.2e1) * t177) * t178 + ((Icges(3,5) * t177 + t173 * t246 + t174 * t244 + t179 * t191) * t237 + (-Icges(3,5) * t179 + t173 * t247 + t174 * t245 + t177 * t191) * t235) * t176; m(7) * (t40 ^ 2 + t41 ^ 2 + t9 ^ 2) + m(6) * (t11 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(5) * (t32 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(4) * (t110 ^ 2 + t111 ^ 2 + t42 ^ 2) + m(3) * (t155 ^ 2 * t210 + t84 ^ 2) + (-t172 * t124 - t5 - t6 + (t144 * t247 + t145 * t245 + t220 * t255) * t179) * t179 + (t8 + t7 + t171 * t125 + (t146 * t246 + t147 * t244 + t219 * t254) * t177 + (-t177 * t124 + t179 * t125 - t246 * t144 - t244 * t145 - t247 * t146 - t245 * t147 - t219 * t255 - t220 * t254) * t179) * t177; ((t177 * t27 + t179 * t26) * t238 + (t177 * t37 + t179 * t36) * t239 + (t177 * t66 + t179 * t65) * t241 + (t177 * t44 + t179 * t43) * t240) * t242; m(7) * (-t178 * t9 + (t177 * t40 + t179 * t41) * t176) + m(6) * (-t178 * t11 + (t177 * t45 + t179 * t46) * t176) + m(5) * (-t178 * t32 + (t177 * t71 + t179 * t72) * t176) + m(4) * (-t178 * t42 + (t110 * t177 + t111 * t179) * t176); 0.2e1 * (t241 + t199) * (t170 * t210 + t178 ^ 2); m(7) * (t144 * t27 + t146 * t26) + m(6) * (t144 * t37 + t146 * t36) + m(5) * (t144 * t44 + t146 * t43); m(7) * (t144 * t40 + t146 * t41 + t222 * t9) + m(6) * (t11 * t222 + t144 * t45 + t146 * t46) + m(5) * (t144 * t71 + t146 * t72 + t222 * t32); t199 * (t144 * t177 + t146 * t179 - t173 * t178) * t242; 0.2e1 * t199 * (t170 * t173 ^ 2 + t144 ^ 2 + t146 ^ 2); m(7) * (t24 * t26 + t25 * t27) + m(6) * (t36 * t38 + t37 * t39) + (t177 * t184 + t179 * t185) * t176 + t232; m(7) * (t10 * t9 + t24 * t41 + t25 * t40) + m(6) * (t11 * t33 + t38 * t46 + t39 * t45) + ((-t8 / 0.2e1 - t7 / 0.2e1) * t179 + (-t6 / 0.2e1 - t5 / 0.2e1) * t177) * t176 + t250 * t237 + (t177 * t248 + t179 * t249) * t178 / 0.2e1 + t251 * t235; m(6) * (-t33 * t178 + (t177 * t39 + t179 * t38) * t176) + m(7) * (-t10 * t178 + (t177 * t25 + t179 * t24) * t176); m(6) * (t144 * t39 + t146 * t38 + t222 * t33) + m(7) * (t10 * t222 + t144 * t25 + t146 * t24); t232 * t178 + m(7) * (t10 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t33 ^ 2 + t38 ^ 2 + t39 ^ 2) + ((-t248 * t178 - t250) * t179 + (t249 * t178 - t251) * t177) * t176; m(7) * (t106 * t27 + t108 * t26); m(7) * (t106 * t40 + t108 * t41 + t132 * t9); m(7) * (-t132 * t178 + (t106 * t177 + t108 * t179) * t176); m(7) * (t106 * t144 + t108 * t146 + t132 * t222); m(7) * (t10 * t132 + t106 * t25 + t108 * t24); m(7) * (t106 ^ 2 + t108 ^ 2 + t132 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
