% Calculate joint inertia matrix for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:12:14
% EndTime: 2019-03-09 05:12:22
% DurationCPUTime: 4.18s
% Computational Cost: add. (6177->326), mult. (5457->491), div. (0->0), fcn. (5560->10), ass. (0->159)
t257 = Icges(5,4) + Icges(6,6);
t256 = Icges(5,1) + Icges(6,2);
t255 = -Icges(5,2) - Icges(6,3);
t156 = pkin(10) + qJ(3);
t149 = qJ(4) + t156;
t145 = cos(t149);
t254 = t257 * t145;
t144 = sin(t149);
t253 = t257 * t144;
t252 = Icges(6,4) - Icges(5,5);
t251 = Icges(6,5) - Icges(5,6);
t250 = t144 * t255 + t254;
t249 = t145 * t256 - t253;
t163 = sin(qJ(1));
t165 = cos(qJ(1));
t248 = t163 * t250 + t165 * t251;
t247 = -t163 * t251 + t165 * t250;
t246 = t163 * t249 + t165 * t252;
t245 = -t163 * t252 + t165 * t249;
t244 = Icges(6,1) + Icges(5,3);
t243 = t144 * t251 - t145 * t252;
t242 = t145 * t255 - t253;
t241 = t144 * t256 + t254;
t240 = -t163 * t243 + t165 * t244;
t239 = t163 * t244 + t165 * t243;
t238 = t144 * t248 - t145 * t246;
t237 = t144 * t247 - t145 * t245;
t157 = t163 ^ 2;
t212 = t145 * t165;
t164 = cos(qJ(6));
t208 = t164 * t165;
t162 = sin(qJ(6));
t210 = t163 * t162;
t110 = t144 * t208 - t210;
t209 = t163 * t164;
t211 = t162 * t165;
t111 = t144 * t211 + t209;
t62 = rSges(7,1) * t111 + rSges(7,2) * t110 + rSges(7,3) * t212;
t236 = pkin(5) * t163 + pkin(9) * t212 + t62;
t235 = -t144 * t252 - t145 * t251;
t234 = t144 * t242 + t145 * t241;
t158 = t165 ^ 2;
t112 = t144 * t209 + t211;
t113 = t144 * t210 - t208;
t213 = t145 * t163;
t56 = Icges(7,5) * t111 + Icges(7,6) * t110 + Icges(7,3) * t212;
t58 = Icges(7,4) * t111 + Icges(7,2) * t110 + Icges(7,6) * t212;
t60 = Icges(7,1) * t111 + Icges(7,4) * t110 + Icges(7,5) * t212;
t18 = t112 * t58 + t113 * t60 + t213 * t56;
t57 = Icges(7,5) * t113 + Icges(7,6) * t112 + Icges(7,3) * t213;
t59 = Icges(7,4) * t113 + Icges(7,2) * t112 + Icges(7,6) * t213;
t61 = Icges(7,1) * t113 + Icges(7,4) * t112 + Icges(7,5) * t213;
t19 = t112 * t59 + t113 * t61 + t213 * t57;
t9 = t163 * t18 - t165 * t19;
t233 = -t9 + t240 * t158 + (t237 * t163 + (-t238 + t239) * t165) * t163;
t232 = m(6) / 0.2e1;
t231 = m(7) / 0.2e1;
t230 = t163 / 0.2e1;
t229 = -t165 / 0.2e1;
t147 = sin(t156);
t228 = pkin(3) * t147;
t161 = -pkin(7) - qJ(2);
t160 = cos(pkin(10));
t146 = pkin(2) * t160 + pkin(1);
t148 = cos(t156);
t132 = pkin(3) * t148 + t146;
t126 = t165 * t132;
t227 = t165 * (-t146 * t165 + t126) + (t132 - t146) * t157;
t214 = t144 * t165;
t172 = rSges(5,1) * t212 - rSges(5,2) * t214 + rSges(5,3) * t163;
t192 = rSges(5,1) * t145 - rSges(5,2) * t144;
t52 = t163 * (-rSges(5,3) * t165 + t163 * t192) + t165 * t172;
t226 = rSges(4,1) * t148;
t225 = rSges(4,2) * t147;
t24 = t144 * t56 + (-t162 * t60 - t164 * t58) * t145;
t224 = t24 * t163;
t25 = t144 * t57 + (-t162 * t61 - t164 * t59) * t145;
t223 = t25 * t165;
t222 = rSges(3,3) + qJ(2);
t221 = Icges(4,4) * t147;
t220 = Icges(4,4) * t148;
t215 = qJ(5) * t144;
t205 = pkin(4) * t212 + qJ(5) * t214;
t207 = t157 * (pkin(4) * t145 + t215) + t165 * t205;
t123 = pkin(4) * t144 - qJ(5) * t145;
t206 = rSges(6,2) * t144 + rSges(6,3) * t145 - t123;
t204 = rSges(4,3) * t163 + t165 * t226;
t201 = t157 + t158;
t16 = t110 * t58 + t111 * t60 + t212 * t56;
t17 = t110 * t59 + t111 * t61 + t212 * t57;
t8 = t16 * t163 - t165 * t17;
t200 = (t8 + t239 * t157 + ((-t237 + t240) * t163 + t238 * t165) * t165) * t163;
t199 = t232 + t231;
t125 = rSges(5,1) * t144 + rSges(5,2) * t145;
t198 = -t125 - t228;
t155 = -pkin(8) + t161;
t197 = -t155 * t163 + t126;
t171 = rSges(6,1) * t163 - rSges(6,2) * t212 + rSges(6,3) * t214;
t35 = t163 * (-rSges(6,1) * t165 + (-rSges(6,2) * t145 + rSges(6,3) * t144) * t163) + t165 * t171 + t207;
t75 = Icges(7,3) * t144 + (-Icges(7,5) * t162 - Icges(7,6) * t164) * t145;
t76 = Icges(7,6) * t144 + (-Icges(7,4) * t162 - Icges(7,2) * t164) * t145;
t77 = Icges(7,5) * t144 + (-Icges(7,1) * t162 - Icges(7,4) * t164) * t145;
t29 = t110 * t76 + t111 * t77 + t212 * t75;
t3 = t29 * t144 + (t16 * t165 + t163 * t17) * t145;
t30 = t112 * t76 + t113 * t77 + t213 * t75;
t4 = t30 * t144 + (t163 * t19 + t165 * t18) * t145;
t196 = t4 * t229 + t3 * t230 + t144 * (-t223 + t224) / 0.2e1 + t9 * t213 / 0.2e1 + t8 * t212 / 0.2e1;
t78 = rSges(7,3) * t144 + (-rSges(7,1) * t162 - rSges(7,2) * t164) * t145;
t195 = -pkin(9) * t144 - t123 - t78;
t194 = t206 - t228;
t193 = -t225 + t226;
t191 = -t113 * rSges(7,1) - t112 * rSges(7,2);
t186 = -t162 * t77 - t164 * t76;
t185 = Icges(4,1) * t148 - t221;
t183 = -Icges(4,2) * t147 + t220;
t180 = Icges(4,5) * t148 - Icges(4,6) * t147;
t159 = sin(pkin(10));
t170 = rSges(3,1) * t160 - rSges(3,2) * t159 + pkin(1);
t169 = t197 + t205;
t63 = rSges(7,3) * t213 - t191;
t20 = t207 + t236 * t165 + (-pkin(5) * t165 + pkin(9) * t213 + t63) * t163;
t168 = t195 - t228;
t167 = t165 * t233 + t200;
t166 = t224 / 0.2e1 - t223 / 0.2e1 + (t144 * t245 + t145 * t247 + t235 * t163 + t234 * t165 + t29) * t230 + (t144 * t246 + t145 * t248 + t234 * t163 - t235 * t165 + t30) * t229;
t140 = rSges(2,1) * t165 - rSges(2,2) * t163;
t139 = -rSges(2,1) * t163 - rSges(2,2) * t165;
t131 = rSges(4,1) * t147 + rSges(4,2) * t148;
t103 = Icges(4,3) * t163 + t165 * t180;
t102 = -Icges(4,3) * t165 + t163 * t180;
t86 = t163 * t222 + t165 * t170;
t85 = -t163 * t170 + t165 * t222;
t80 = t198 * t165;
t79 = t198 * t163;
t74 = t144 * t75;
t73 = -t163 * t161 + (t146 - t225) * t165 + t204;
t72 = (rSges(4,3) - t161) * t165 + (-t146 - t193) * t163;
t71 = t206 * t165;
t70 = t206 * t163;
t67 = t172 + t197;
t66 = (rSges(5,3) - t155) * t165 + (-t132 - t192) * t163;
t65 = t194 * t165;
t64 = t194 * t163;
t55 = t165 * (-t165 * t225 + t204) + (-t165 * rSges(4,3) + t163 * t193) * t163;
t51 = t195 * t165;
t50 = t195 * t163;
t49 = t169 + t171;
t48 = (rSges(6,1) - t155) * t165 + (-t132 + (rSges(6,2) - pkin(4)) * t145 + (-rSges(6,3) - qJ(5)) * t144) * t163;
t43 = t168 * t165;
t42 = t168 * t163;
t37 = t144 * t62 - t212 * t78;
t36 = -t144 * t63 + t213 * t78;
t34 = t169 + t236;
t33 = (pkin(5) - t155) * t165 + (-t215 - t132 + (-rSges(7,3) - pkin(4) - pkin(9)) * t145) * t163 + t191;
t32 = (t145 * t186 + t74) * t144;
t31 = (-t163 * t62 + t165 * t63) * t145;
t26 = t52 + t227;
t21 = t35 + t227;
t11 = t20 + t227;
t1 = [Icges(3,2) * t160 ^ 2 + t148 * (Icges(4,2) * t148 + t221) + t147 * (Icges(4,1) * t147 + t220) + Icges(2,3) + t74 + (Icges(3,1) * t159 + 0.2e1 * Icges(3,4) * t160) * t159 + t241 * t144 + (t186 - t242) * t145 + m(7) * (t33 ^ 2 + t34 ^ 2) + m(6) * (t48 ^ 2 + t49 ^ 2) + m(5) * (t66 ^ 2 + t67 ^ 2) + m(4) * (t72 ^ 2 + t73 ^ 2) + m(3) * (t85 ^ 2 + t86 ^ 2) + m(2) * (t139 ^ 2 + t140 ^ 2); m(7) * (t163 * t33 - t165 * t34) + m(6) * (t163 * t48 - t165 * t49) + m(5) * (t163 * t66 - t165 * t67) + m(4) * (t163 * t72 - t165 * t73) + m(3) * (t163 * t85 - t165 * t86); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + t199) * t201; ((-Icges(4,6) * t165 + t163 * t183) * t148 + (-Icges(4,5) * t165 + t163 * t185) * t147) * t229 + ((Icges(4,6) * t163 + t165 * t183) * t148 + (Icges(4,5) * t163 + t165 * t185) * t147) * t230 + m(4) * (-t163 * t73 - t165 * t72) * t131 + (t157 / 0.2e1 + t158 / 0.2e1) * (Icges(4,5) * t147 + Icges(4,6) * t148) + t166 + m(7) * (t33 * t43 + t34 * t42) + m(6) * (t48 * t65 + t49 * t64) + m(5) * (t66 * t80 + t67 * t79); m(5) * (t163 * t80 - t165 * t79) + m(6) * (t163 * t65 - t165 * t64) + m(7) * (t163 * t43 - t165 * t42); m(7) * (t11 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(6) * (t21 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(5) * (t26 ^ 2 + t79 ^ 2 + t80 ^ 2) + t163 * t157 * t103 + m(4) * (t131 ^ 2 * t201 + t55 ^ 2) + t200 + (-t158 * t102 + (-t163 * t102 + t165 * t103) * t163 + t233) * t165; m(7) * (t33 * t51 + t34 * t50) + m(6) * (t48 * t71 + t49 * t70) + m(5) * (-t163 * t67 - t165 * t66) * t125 + t166; m(6) * (t163 * t71 - t165 * t70) + m(7) * (t163 * t51 - t165 * t50); m(7) * (t11 * t20 + t42 * t50 + t43 * t51) + m(6) * (t21 * t35 + t64 * t70 + t65 * t71) + m(5) * (t52 * t26 + (-t163 * t79 - t165 * t80) * t125) + t167; m(7) * (t20 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t35 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(5) * (t125 ^ 2 * t201 + t52 ^ 2) + t167; 0.2e1 * ((t163 * t34 + t165 * t33) * t231 + (t163 * t49 + t165 * t48) * t232) * t144; 0; m(7) * (-t145 * t11 + (t163 * t42 + t165 * t43) * t144) + m(6) * (-t145 * t21 + (t163 * t64 + t165 * t65) * t144); m(7) * (-t145 * t20 + (t163 * t50 + t165 * t51) * t144) + m(6) * (-t145 * t35 + (t163 * t70 + t165 * t71) * t144); 0.2e1 * t199 * (t144 ^ 2 * t201 + t145 ^ 2); m(7) * (t33 * t36 + t34 * t37) + t32 + ((t24 / 0.2e1 + t29 / 0.2e1) * t165 + (t25 / 0.2e1 + t30 / 0.2e1) * t163) * t145; m(7) * (t163 * t36 - t165 * t37); m(7) * (t11 * t31 + t36 * t43 + t37 * t42) + t196; m(7) * (t20 * t31 + t36 * t51 + t37 * t50) + t196; m(7) * (-t31 * t145 + (t163 * t37 + t165 * t36) * t144); t144 * t32 + m(7) * (t31 ^ 2 + t36 ^ 2 + t37 ^ 2) + (t165 * t3 + t163 * t4 + t144 * (t163 * t25 + t165 * t24)) * t145;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
