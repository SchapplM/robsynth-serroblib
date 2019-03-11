% Calculate joint inertia matrix for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:21:22
% EndTime: 2019-03-09 08:21:29
% DurationCPUTime: 2.81s
% Computational Cost: add. (3133->386), mult. (7667->563), div. (0->0), fcn. (8764->8), ass. (0->169)
t236 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t235 = Icges(6,1) + Icges(4,2) + Icges(5,3);
t234 = Icges(4,4) - Icges(6,5) + Icges(5,6);
t233 = Icges(5,4) - Icges(4,5) + Icges(6,6);
t232 = Icges(6,4) - Icges(5,5) + Icges(4,6);
t156 = sin(qJ(1));
t217 = t156 / 0.2e1;
t159 = cos(qJ(1));
t215 = t159 / 0.2e1;
t218 = m(7) / 0.2e1;
t219 = m(6) / 0.2e1;
t190 = t219 + t218;
t231 = 0.2e1 * t190;
t153 = cos(pkin(9));
t158 = cos(qJ(2));
t152 = sin(pkin(9));
t198 = t159 * t152;
t122 = -t156 * t153 + t158 * t198;
t155 = sin(qJ(2));
t201 = t155 * t159;
t199 = t158 * t159;
t123 = t156 * t152 + t153 * t199;
t154 = sin(qJ(6));
t157 = cos(qJ(6));
t79 = -t122 * t154 + t123 * t157;
t80 = t122 * t157 + t123 * t154;
t36 = t80 * rSges(7,1) + t79 * rSges(7,2) + rSges(7,3) * t201;
t230 = t122 * pkin(5) + pkin(8) * t201 + t36;
t200 = t156 * t158;
t120 = t152 * t200 + t153 * t159;
t121 = t153 * t200 - t198;
t202 = t155 * t156;
t229 = t235 * t120 - t234 * t121 - t232 * t202;
t228 = t234 * t120 - t236 * t121 + t233 * t202;
t57 = Icges(5,1) * t202 - Icges(5,4) * t121 + Icges(5,5) * t120;
t61 = Icges(4,5) * t121 - Icges(4,6) * t120 + Icges(4,3) * t202;
t63 = Icges(6,4) * t120 - Icges(6,2) * t202 + Icges(6,6) * t121;
t227 = t57 + t61 - t63;
t226 = -t234 * t122 + t236 * t123 - t233 * t201;
t58 = Icges(5,1) * t201 - Icges(5,4) * t123 + Icges(5,5) * t122;
t62 = Icges(4,5) * t123 - Icges(4,6) * t122 + Icges(4,3) * t201;
t64 = Icges(6,4) * t122 - Icges(6,2) * t201 + Icges(6,6) * t123;
t225 = -t64 + t62 + t58;
t224 = t235 * t122 - t234 * t123 - t232 * t201;
t150 = t156 ^ 2;
t223 = t158 ^ 2;
t151 = t159 ^ 2;
t222 = 0.2e1 * t155;
t221 = m(4) / 0.2e1;
t220 = m(5) / 0.2e1;
t216 = -t159 / 0.2e1;
t214 = pkin(2) * t158;
t213 = pkin(5) * t120;
t212 = -pkin(4) - qJ(3);
t105 = (-t152 * t154 + t153 * t157) * t155;
t106 = (t152 * t157 + t153 * t154) * t155;
t49 = Icges(7,4) * t106 + Icges(7,2) * t105 - Icges(7,6) * t158;
t50 = Icges(7,1) * t106 + Icges(7,4) * t105 - Icges(7,5) * t158;
t211 = t105 * t49 + t106 * t50;
t210 = t120 * rSges(6,1);
t209 = t120 * rSges(5,3);
t208 = t159 * rSges(3,3);
t131 = pkin(2) * t155 - qJ(3) * t158;
t207 = -t131 + t158 * rSges(4,3) - (rSges(4,1) * t153 - rSges(4,2) * t152) * t155;
t206 = Icges(3,4) * t155;
t205 = Icges(3,4) * t158;
t204 = t152 * t155;
t203 = t153 * t155;
t194 = pkin(2) * t199 + qJ(3) * t201;
t197 = t150 * (qJ(3) * t155 + t214) + t159 * t194;
t195 = -(pkin(3) * t153 + qJ(4) * t152) * t155 - t131;
t193 = t159 * pkin(1) + t156 * pkin(7);
t109 = t120 * qJ(4);
t147 = t159 * pkin(7);
t192 = t147 - t109;
t191 = t150 + t151;
t77 = -t120 * t154 + t121 * t157;
t78 = t120 * t157 + t121 * t154;
t29 = Icges(7,5) * t78 + Icges(7,6) * t77 + Icges(7,3) * t202;
t31 = Icges(7,4) * t78 + Icges(7,2) * t77 + Icges(7,6) * t202;
t33 = Icges(7,1) * t78 + Icges(7,4) * t77 + Icges(7,5) * t202;
t10 = t105 * t31 + t106 * t33 - t158 * t29;
t48 = Icges(7,5) * t106 + Icges(7,6) * t105 - Icges(7,3) * t158;
t13 = t48 * t202 + t49 * t77 + t50 * t78;
t189 = t10 / 0.2e1 + t13 / 0.2e1;
t30 = Icges(7,5) * t80 + Icges(7,6) * t79 + Icges(7,3) * t201;
t32 = Icges(7,4) * t80 + Icges(7,2) * t79 + Icges(7,6) * t201;
t34 = Icges(7,1) * t80 + Icges(7,4) * t79 + Icges(7,5) * t201;
t11 = t105 * t32 + t106 * t34 - t158 * t30;
t14 = t48 * t201 + t79 * t49 + t80 * t50;
t188 = t11 / 0.2e1 + t14 / 0.2e1;
t187 = t158 * rSges(5,1) - (-rSges(5,2) * t153 + rSges(5,3) * t152) * t155 + t195;
t186 = t123 * rSges(4,1) - t122 * rSges(4,2) + rSges(4,3) * t201;
t185 = pkin(4) * t158 - qJ(5) * t203 + t195;
t184 = rSges(5,1) * t201 - t123 * rSges(5,2) + t122 * rSges(5,3);
t183 = -pkin(1) - t214;
t182 = t121 * pkin(3) + t109;
t181 = t123 * pkin(3) + t122 * qJ(4);
t180 = pkin(4) * t201 + t123 * qJ(5);
t179 = t156 * t182 + t159 * t181 + t197;
t178 = t220 + t190;
t177 = -t158 * rSges(6,2) - (rSges(6,1) * t152 + rSges(6,3) * t153) * t155 + t185;
t176 = t193 + t194;
t90 = -Icges(4,6) * t158 + (Icges(4,4) * t153 - Icges(4,2) * t152) * t155;
t91 = Icges(6,4) * t158 + (Icges(6,1) * t152 + Icges(6,5) * t153) * t155;
t93 = -Icges(5,5) * t158 + (-Icges(5,6) * t153 + Icges(5,3) * t152) * t155;
t175 = -t93 / 0.2e1 - t91 / 0.2e1 + t90 / 0.2e1;
t87 = Icges(6,6) * t158 + (Icges(6,5) * t152 + Icges(6,3) * t153) * t155;
t92 = -Icges(4,5) * t158 + (Icges(4,1) * t153 - Icges(4,4) * t152) * t155;
t94 = -Icges(5,4) * t158 + (-Icges(5,2) * t153 + Icges(5,6) * t152) * t155;
t174 = t94 / 0.2e1 - t87 / 0.2e1 - t92 / 0.2e1;
t173 = -rSges(7,1) * t78 - rSges(7,2) * t77;
t51 = rSges(7,1) * t106 + rSges(7,2) * t105 - rSges(7,3) * t158;
t172 = -pkin(5) * t204 + pkin(8) * t158 + t185 - t51;
t171 = rSges(3,1) * t158 - rSges(3,2) * t155;
t170 = -t121 * rSges(4,1) + t120 * rSges(4,2);
t169 = Icges(3,1) * t158 - t206;
t168 = -Icges(3,2) * t155 + t205;
t167 = Icges(3,5) * t158 - Icges(3,6) * t155;
t110 = t121 * qJ(5);
t164 = t156 * (pkin(4) * t202 + t110) + t159 * t180 + t179;
t163 = rSges(3,1) * t199 - rSges(3,2) * t201 + t156 * rSges(3,3);
t162 = t122 * rSges(6,1) - rSges(6,2) * t201 + t123 * rSges(6,3);
t161 = t176 + t181;
t160 = t161 + t180;
t149 = t155 ^ 2;
t134 = rSges(2,1) * t159 - t156 * rSges(2,2);
t133 = -t156 * rSges(2,1) - rSges(2,2) * t159;
t132 = rSges(3,1) * t155 + rSges(3,2) * t158;
t128 = Icges(3,5) * t155 + Icges(3,6) * t158;
t100 = Icges(3,3) * t156 + t167 * t159;
t99 = -Icges(3,3) * t159 + t167 * t156;
t86 = t163 + t193;
t85 = t208 + t147 + (-pkin(1) - t171) * t156;
t82 = t207 * t159;
t81 = t207 * t156;
t52 = t159 * t163 + (t171 * t156 - t208) * t156;
t45 = t187 * t159;
t44 = t187 * t156;
t41 = t176 + t186;
t40 = t147 + ((-rSges(4,3) - qJ(3)) * t155 + t183) * t156 + t170;
t39 = t177 * t159;
t38 = t177 * t156;
t35 = rSges(7,3) * t202 - t173;
t28 = t161 + t184;
t27 = -t209 + (rSges(5,2) - pkin(3)) * t121 + ((-rSges(5,1) - qJ(3)) * t155 + t183) * t156 + t192;
t26 = t172 * t159;
t25 = t172 * t156;
t24 = t160 + t162;
t23 = -t210 - t110 + (-rSges(6,3) - pkin(3)) * t121 + ((rSges(6,2) + t212) * t155 + t183) * t156 + t192;
t22 = t156 * (rSges(4,3) * t202 - t170) + t159 * t186 + t197;
t21 = -t158 * t36 - t51 * t201;
t20 = t158 * t35 + t51 * t202;
t19 = -t158 * t48 + t211;
t18 = (-t156 * t36 + t159 * t35) * t155;
t17 = t160 + t230;
t16 = -t213 - t110 + t147 + ((-rSges(7,3) - pkin(8) + t212) * t155 + t183) * t156 + t173 - t182;
t15 = t156 * (rSges(5,1) * t202 - t121 * rSges(5,2) + t209) + t159 * t184 + t179;
t12 = t156 * (-rSges(6,2) * t202 + t121 * rSges(6,3) + t210) + t159 * t162 + t164;
t9 = t30 * t201 + t79 * t32 + t80 * t34;
t8 = t29 * t201 + t79 * t31 + t80 * t33;
t7 = t30 * t202 + t32 * t77 + t34 * t78;
t6 = t29 * t202 + t31 * t77 + t33 * t78;
t5 = t230 * t159 + (pkin(8) * t202 + t213 + t35) * t156 + t164;
t4 = t9 * t156 - t159 * t8;
t3 = t7 * t156 - t159 * t6;
t2 = -t14 * t158 + (t156 * t8 + t159 * t9) * t155;
t1 = -t13 * t158 + (t156 * t6 + t159 * t7) * t155;
t37 = [Icges(2,3) + (-t48 + t206 + (t232 * t152 + t233 * t153) * t155 + (Icges(3,2) + Icges(4,3) + Icges(5,1) + Icges(6,2)) * t158) * t158 + (Icges(3,1) * t155 + t205 + (t87 + t92 - t94) * t153 + (-t90 + t91 + t93) * t152) * t155 + m(7) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2) + m(6) * (t23 ^ 2 + t24 ^ 2) + m(4) * (t40 ^ 2 + t41 ^ 2) + m(3) * (t85 ^ 2 + t86 ^ 2) + m(2) * (t133 ^ 2 + t134 ^ 2) + t211; (t175 * t120 + t174 * t121 + t128 * t215 - t189) * t159 + (-t175 * t122 - t174 * t123 + t128 * t217 + t188) * t156 + m(7) * (t16 * t26 + t17 * t25) + m(6) * (t23 * t39 + t24 * t38) + m(5) * (t27 * t45 + t28 * t44) + m(4) * (t40 * t82 + t41 * t81) + m(3) * (-t156 * t86 - t159 * t85) * t132 + ((Icges(3,6) * t215 - t168 * t156 / 0.2e1 - t63 / 0.2e1 + t61 / 0.2e1 + t57 / 0.2e1) * t159 + (-t62 / 0.2e1 + t64 / 0.2e1 - t58 / 0.2e1 + Icges(3,6) * t217 + t168 * t215) * t156) * t158 + ((Icges(3,5) * t156 + t152 * t224 + t153 * t226 + t169 * t159) * t217 + (-Icges(3,5) * t159 + t152 * t229 - t153 * t228 + t169 * t156) * t216) * t155; m(7) * (t25 ^ 2 + t26 ^ 2 + t5 ^ 2) + m(6) * (t12 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(5) * (t15 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(4) * (t22 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(3) * (t191 * t132 ^ 2 + t52 ^ 2) + (-t151 * t99 - t3 + (t120 * t229 - t121 * t228 + t202 * t227) * t159) * t159 + (t150 * t100 + t4 + (t159 * t100 - t120 * t224 - t121 * t226 - t122 * t229 + t123 * t228 - t201 * t227 - t202 * t225) * t159 + (t122 * t224 + t123 * t226 - t159 * t99 + t201 * t225) * t156) * t156; ((t156 * t17 + t159 * t16) * t218 + (t156 * t28 + t159 * t27) * t220 + (t156 * t24 + t159 * t23) * t219 + (t156 * t41 + t159 * t40) * t221) * t222; m(7) * (-t158 * t5 + (t156 * t25 + t159 * t26) * t155) + m(6) * (-t158 * t12 + (t156 * t38 + t159 * t39) * t155) + m(5) * (-t158 * t15 + (t156 * t44 + t159 * t45) * t155) + m(4) * (-t158 * t22 + (t156 * t81 + t159 * t82) * t155); 0.2e1 * (t221 + t178) * (t191 * t149 + t223); m(7) * (t120 * t17 + t122 * t16) + m(5) * (t120 * t28 + t122 * t27) + m(6) * (t120 * t24 + t122 * t23); m(7) * (t120 * t25 + t122 * t26 + t5 * t204) + m(6) * (t12 * t204 + t120 * t38 + t122 * t39) + m(5) * (t120 * t44 + t122 * t45 + t15 * t204); t178 * (t120 * t156 + t122 * t159 - t152 * t158) * t222; 0.2e1 * t178 * (t149 * t152 ^ 2 + t120 ^ 2 + t122 ^ 2); m(7) * (t121 * t17 + t123 * t16) + m(6) * (t121 * t24 + t123 * t23); m(7) * (t121 * t25 + t123 * t26 + t5 * t203) + m(6) * (t12 * t203 + t121 * t38 + t123 * t39); t190 * (t121 * t156 + t123 * t159 - t153 * t158) * t222; (t149 * t152 * t153 + t120 * t121 + t122 * t123) * t231; (t149 * t153 ^ 2 + t121 ^ 2 + t123 ^ 2) * t231; m(7) * (t16 * t20 + t17 * t21) - t19 * t158 + (t189 * t156 + t188 * t159) * t155; m(7) * (t18 * t5 + t20 * t26 + t21 * t25) + t1 * t216 + t2 * t217 - t158 * (-t10 * t159 + t11 * t156) / 0.2e1 + (t4 * t215 + t3 * t217) * t155; m(7) * (-t18 * t158 + (t156 * t21 + t159 * t20) * t155); m(7) * (t120 * t21 + t122 * t20 + t18 * t204); m(7) * (t121 * t21 + t123 * t20 + t18 * t203); t223 * t19 + m(7) * (t18 ^ 2 + t20 ^ 2 + t21 ^ 2) + (t159 * t2 + t156 * t1 - t158 * (t10 * t156 + t11 * t159)) * t155;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t37(1) t37(2) t37(4) t37(7) t37(11) t37(16); t37(2) t37(3) t37(5) t37(8) t37(12) t37(17); t37(4) t37(5) t37(6) t37(9) t37(13) t37(18); t37(7) t37(8) t37(9) t37(10) t37(14) t37(19); t37(11) t37(12) t37(13) t37(14) t37(15) t37(20); t37(16) t37(17) t37(18) t37(19) t37(20) t37(21);];
Mq  = res;
