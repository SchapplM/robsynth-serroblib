% Calculate joint inertia matrix for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:00:00
% EndTime: 2019-12-31 21:00:07
% DurationCPUTime: 2.64s
% Computational Cost: add. (4608->370), mult. (7098->539), div. (0->0), fcn. (7663->8), ass. (0->188)
t174 = sin(qJ(2));
t248 = Icges(3,5) * t174;
t247 = t248 / 0.2e1;
t168 = qJ(3) + pkin(8);
t162 = sin(t168);
t163 = cos(t168);
t178 = cos(qJ(1));
t210 = t178 * t163;
t175 = sin(qJ(1));
t177 = cos(qJ(2));
t213 = t175 * t177;
t118 = t162 * t213 + t210;
t211 = t178 * t162;
t119 = t163 * t213 - t211;
t244 = rSges(6,3) + qJ(5);
t245 = rSges(6,1) + pkin(4);
t246 = -t244 * t118 - t245 * t119;
t106 = -Icges(5,3) * t177 + (Icges(5,5) * t163 - Icges(5,6) * t162) * t174;
t107 = -Icges(6,2) * t177 + (Icges(6,4) * t163 + Icges(6,6) * t162) * t174;
t173 = sin(qJ(3));
t176 = cos(qJ(3));
t122 = -Icges(4,3) * t177 + (Icges(4,5) * t176 - Icges(4,6) * t173) * t174;
t243 = t106 + t107 + t122;
t216 = t174 * t175;
t62 = Icges(6,5) * t119 + Icges(6,6) * t216 + Icges(6,3) * t118;
t66 = Icges(6,4) * t119 + Icges(6,2) * t216 + Icges(6,6) * t118;
t70 = Icges(6,1) * t119 + Icges(6,4) * t216 + Icges(6,5) * t118;
t19 = t118 * t62 + t119 * t70 + t216 * t66;
t120 = -t175 * t163 + t177 * t211;
t121 = t175 * t162 + t177 * t210;
t215 = t174 * t178;
t63 = Icges(6,5) * t121 + Icges(6,6) * t215 + Icges(6,3) * t120;
t67 = Icges(6,4) * t121 + Icges(6,2) * t215 + Icges(6,6) * t120;
t71 = Icges(6,1) * t121 + Icges(6,4) * t215 + Icges(6,5) * t120;
t20 = t118 * t63 + t119 * t71 + t216 * t67;
t64 = Icges(5,5) * t119 - Icges(5,6) * t118 + Icges(5,3) * t216;
t68 = Icges(5,4) * t119 - Icges(5,2) * t118 + Icges(5,6) * t216;
t72 = Icges(5,1) * t119 - Icges(5,4) * t118 + Icges(5,5) * t216;
t21 = -t118 * t68 + t119 * t72 + t216 * t64;
t65 = Icges(5,5) * t121 - Icges(5,6) * t120 + Icges(5,3) * t215;
t69 = Icges(5,4) * t121 - Icges(5,2) * t120 + Icges(5,6) * t215;
t73 = Icges(5,1) * t121 - Icges(5,4) * t120 + Icges(5,5) * t215;
t22 = -t118 * t69 + t119 * t73 + t216 * t65;
t208 = t178 * t176;
t138 = -t173 * t213 - t208;
t209 = t178 * t173;
t139 = t176 * t213 - t209;
t85 = Icges(4,5) * t139 + Icges(4,6) * t138 + Icges(4,3) * t216;
t87 = Icges(4,4) * t139 + Icges(4,2) * t138 + Icges(4,6) * t216;
t89 = Icges(4,1) * t139 + Icges(4,4) * t138 + Icges(4,5) * t216;
t31 = t138 * t87 + t139 * t89 + t216 * t85;
t140 = t175 * t176 - t177 * t209;
t214 = t175 * t173;
t141 = t177 * t208 + t214;
t86 = Icges(4,5) * t141 + Icges(4,6) * t140 + Icges(4,3) * t215;
t88 = Icges(4,4) * t141 + Icges(4,2) * t140 + Icges(4,6) * t215;
t90 = Icges(4,1) * t141 + Icges(4,4) * t140 + Icges(4,5) * t215;
t32 = t138 * t88 + t139 * t90 + t216 * t86;
t105 = -Icges(6,6) * t177 + (Icges(6,5) * t163 + Icges(6,3) * t162) * t174;
t109 = -Icges(6,4) * t177 + (Icges(6,1) * t163 + Icges(6,5) * t162) * t174;
t42 = t118 * t105 + t107 * t216 + t119 * t109;
t108 = -Icges(5,6) * t177 + (Icges(5,4) * t163 - Icges(5,2) * t162) * t174;
t110 = -Icges(5,5) * t177 + (Icges(5,1) * t163 - Icges(5,4) * t162) * t174;
t43 = t106 * t216 - t118 * t108 + t119 * t110;
t125 = -Icges(4,6) * t177 + (Icges(4,4) * t176 - Icges(4,2) * t173) * t174;
t128 = -Icges(4,5) * t177 + (Icges(4,1) * t176 - Icges(4,4) * t173) * t174;
t46 = t122 * t216 + t138 * t125 + t139 * t128;
t242 = (-t42 - t43 - t46) * t177 + ((t20 + t22 + t32) * t178 + (t19 + t21 + t31) * t175) * t174;
t23 = t120 * t62 + t121 * t70 + t215 * t66;
t24 = t120 * t63 + t121 * t71 + t215 * t67;
t25 = -t120 * t68 + t121 * t72 + t215 * t64;
t26 = -t120 * t69 + t121 * t73 + t215 * t65;
t33 = t140 * t87 + t141 * t89 + t215 * t85;
t34 = t140 * t88 + t141 * t90 + t215 * t86;
t44 = t120 * t105 + t107 * t215 + t121 * t109;
t45 = t106 * t215 - t120 * t108 + t121 * t110;
t47 = t122 * t215 + t140 * t125 + t141 * t128;
t241 = (-t44 - t45 - t47) * t177 + ((t24 + t26 + t34) * t178 + (t23 + t25 + t33) * t175) * t174;
t27 = -t177 * t66 + (t162 * t62 + t163 * t70) * t174;
t29 = -t177 * t64 + (-t162 * t68 + t163 * t72) * t174;
t37 = -t177 * t85 + (-t173 * t87 + t176 * t89) * t174;
t240 = -t27 - t29 - t37;
t28 = -t177 * t67 + (t162 * t63 + t163 * t71) * t174;
t30 = -t177 * t65 + (-t162 * t69 + t163 * t73) * t174;
t38 = -t177 * t86 + (-t173 * t88 + t176 * t90) * t174;
t239 = t28 + t30 + t38;
t219 = t162 * t174;
t238 = t105 * t219 + (t176 * t128 + (t109 + t110) * t163) * t174;
t170 = t175 ^ 2;
t237 = t177 ^ 2;
t171 = t178 ^ 2;
t236 = m(5) / 0.2e1;
t235 = m(6) / 0.2e1;
t234 = t175 / 0.2e1;
t233 = -t177 / 0.2e1;
t147 = t174 * rSges(3,1) + t177 * rSges(3,2);
t231 = m(3) * t147;
t230 = pkin(2) * t177;
t229 = pkin(7) * t174;
t161 = t176 * pkin(3) + pkin(2);
t228 = -pkin(2) + t161;
t227 = rSges(6,2) * t216 - t246;
t226 = rSges(6,2) * t215 + t244 * t120 + t245 * t121;
t77 = t121 * rSges(5,1) - t120 * rSges(5,2) + rSges(5,3) * t215;
t172 = -qJ(4) - pkin(7);
t212 = t177 * t178;
t182 = pkin(3) * t214 + t161 * t212 - t172 * t215;
t202 = pkin(2) * t212 + pkin(7) * t215;
t92 = t182 - t202;
t225 = -t77 - t92;
t223 = t178 * rSges(3,3);
t104 = (pkin(7) + t172) * t177 + t228 * t174;
t203 = -pkin(3) * t209 - t172 * t216;
t91 = (t177 * t228 - t229) * t175 + t203;
t222 = t104 * t216 + t177 * t91;
t220 = Icges(3,4) * t177;
t217 = t173 * t125;
t112 = -t177 * rSges(5,3) + (rSges(5,1) * t163 - rSges(5,2) * t162) * t174;
t207 = -t104 - t112;
t206 = -t177 * rSges(6,2) + (t244 * t162 + t245 * t163) * t174;
t132 = -t177 * rSges(4,3) + (rSges(4,1) * t176 - rSges(4,2) * t173) * t174;
t150 = t174 * pkin(2) - t177 * pkin(7);
t205 = -t132 - t150;
t204 = t170 * (t229 + t230) + t178 * t202;
t201 = t178 * pkin(1) + t175 * pkin(6);
t200 = t170 + t171;
t199 = t108 * t219 + t174 * t217 + t243 * t177 - t238;
t198 = -t92 - t226;
t197 = -t104 - t206;
t196 = -t150 + t207;
t94 = t141 * rSges(4,1) + t140 * rSges(4,2) + rSges(4,3) * t215;
t166 = t178 * pkin(6);
t195 = t166 - t203;
t194 = -t161 * t177 - pkin(1);
t193 = t175 * t91 + t178 * t92 + t204;
t192 = -t150 + t197;
t191 = rSges(3,1) * t177 - rSges(3,2) * t174;
t190 = -t139 * rSges(4,1) - t138 * rSges(4,2);
t189 = -t119 * rSges(5,1) + t118 * rSges(5,2);
t187 = -Icges(3,2) * t174 + t220;
t186 = Icges(3,5) * t177 - Icges(3,6) * t174;
t183 = rSges(3,1) * t212 - rSges(3,2) * t215 + t175 * rSges(3,3);
t181 = t182 + t201;
t180 = t38 / 0.2e1 + t30 / 0.2e1 + t28 / 0.2e1 + t47 / 0.2e1 + t45 / 0.2e1 + t44 / 0.2e1;
t179 = t43 / 0.2e1 + t42 / 0.2e1 + t37 / 0.2e1 + t29 / 0.2e1 + t27 / 0.2e1 + t46 / 0.2e1;
t169 = t174 ^ 2;
t149 = t178 * rSges(2,1) - t175 * rSges(2,2);
t148 = -t175 * rSges(2,1) - t178 * rSges(2,2);
t144 = Icges(3,6) * t177 + t248;
t124 = Icges(3,3) * t175 + t178 * t186;
t123 = -Icges(3,3) * t178 + t175 * t186;
t102 = t183 + t201;
t101 = t223 + t166 + (-pkin(1) - t191) * t175;
t96 = t205 * t178;
t95 = t205 * t175;
t93 = rSges(4,3) * t216 - t190;
t79 = t178 * t183 + (t175 * t191 - t223) * t175;
t78 = t91 * t215;
t75 = rSges(5,3) * t216 - t189;
t61 = t94 + t201 + t202;
t60 = t166 + (-t230 - pkin(1) + (-rSges(4,3) - pkin(7)) * t174) * t175 + t190;
t59 = t196 * t178;
t58 = t196 * t175;
t57 = -t132 * t215 - t177 * t94;
t56 = t132 * t216 + t177 * t93;
t54 = t181 + t77;
t53 = (-rSges(5,3) * t174 + t194) * t175 + t189 + t195;
t52 = t192 * t178;
t51 = t192 * t175;
t48 = (-t175 * t94 + t178 * t93) * t174;
t41 = t175 * t93 + t178 * t94 + t204;
t40 = t181 + t226;
t39 = (-rSges(6,2) * t174 + t194) * t175 + t195 + t246;
t36 = t177 * t225 + t207 * t215;
t35 = t112 * t216 + t177 * t75 + t222;
t18 = t78 + (t175 * t225 + t178 * t75) * t174;
t17 = t177 * t198 + t197 * t215;
t16 = t177 * t227 + t206 * t216 + t222;
t15 = t175 * t75 + t178 * t77 + t193;
t14 = t78 + (t175 * t198 + t178 * t227) * t174;
t13 = t175 * t227 + t178 * t226 + t193;
t12 = t34 * t175 - t33 * t178;
t11 = t32 * t175 - t31 * t178;
t10 = t26 * t175 - t25 * t178;
t9 = t24 * t175 - t23 * t178;
t8 = t22 * t175 - t21 * t178;
t7 = t20 * t175 - t19 * t178;
t1 = [Icges(2,3) + (Icges(3,1) * t174 - t162 * t108 - t217 + t220) * t174 + (Icges(3,4) * t174 + Icges(3,2) * t177 - t243) * t177 + m(5) * (t53 ^ 2 + t54 ^ 2) + m(6) * (t39 ^ 2 + t40 ^ 2) + m(4) * (t60 ^ 2 + t61 ^ 2) + m(3) * (t101 ^ 2 + t102 ^ 2) + m(2) * (t148 ^ 2 + t149 ^ 2) + t238; m(5) * (t59 * t53 + t58 * t54) + m(6) * (t52 * t39 + t51 * t40) + m(4) * (t96 * t60 + t95 * t61) + (t175 * t187 * t233 - t101 * t231 - t179 + (-Icges(3,6) * t233 + t247 + t144 / 0.2e1) * t178) * t178 + (-t102 * t231 + t177 * (Icges(3,6) * t175 + t178 * t187) / 0.2e1 + t175 * t247 + t144 * t234 + t180) * t175; m(6) * (t13 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(5) * (t15 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(4) * (t41 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(3) * (t147 ^ 2 * t200 + t79 ^ 2) + (-t171 * t123 - t11 - t7 - t8) * t178 + (t170 * t124 + t10 + t12 + t9 + (-t175 * t123 + t178 * t124) * t178) * t175; t199 * t177 + m(5) * (t35 * t53 + t36 * t54) + m(6) * (t16 * t39 + t17 * t40) + m(4) * (t56 * t60 + t57 * t61) + (t175 * t179 + t178 * t180) * t174; m(6) * (t14 * t13 + t16 * t52 + t17 * t51) + m(5) * (t18 * t15 + t35 * t59 + t36 * t58) + m(4) * (t48 * t41 + t56 * t96 + t57 * t95) + ((t9 / 0.2e1 + t12 / 0.2e1 + t10 / 0.2e1) * t178 + (t8 / 0.2e1 + t7 / 0.2e1 + t11 / 0.2e1) * t175) * t174 + t241 * t234 + (t175 * t239 + t178 * t240) * t233 - t242 * t178 / 0.2e1; m(5) * (t18 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(6) * (t14 ^ 2 + t16 ^ 2 + t17 ^ 2) + m(4) * (t48 ^ 2 + t56 ^ 2 + t57 ^ 2) - t199 * t237 + ((-t177 * t239 + t241) * t178 + (t177 * t240 + t242) * t175) * t174; 0.2e1 * ((t175 * t54 + t178 * t53) * t236 + (t175 * t40 + t178 * t39) * t235) * t174; m(6) * (-t177 * t13 + (t175 * t51 + t178 * t52) * t174) + m(5) * (-t177 * t15 + (t175 * t58 + t178 * t59) * t174); m(5) * (-t177 * t18 + (t175 * t36 + t178 * t35) * t174) + m(6) * (-t177 * t14 + (t16 * t178 + t17 * t175) * t174); 0.2e1 * (t236 + t235) * (t169 * t200 + t237); m(6) * (t118 * t40 + t120 * t39); m(6) * (t118 * t51 + t120 * t52 + t13 * t219); m(6) * (t118 * t17 + t120 * t16 + t14 * t219); m(6) * (t118 * t175 + t120 * t178 - t162 * t177) * t174; m(6) * (t169 * t162 ^ 2 + t118 ^ 2 + t120 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
