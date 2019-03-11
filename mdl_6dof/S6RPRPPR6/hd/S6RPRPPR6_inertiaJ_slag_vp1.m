% Calculate joint inertia matrix for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:53:17
% EndTime: 2019-03-09 02:53:22
% DurationCPUTime: 2.01s
% Computational Cost: add. (3523->324), mult. (3791->476), div. (0->0), fcn. (3843->10), ass. (0->159)
t224 = Icges(4,3) + Icges(5,3);
t134 = qJ(3) + pkin(9);
t125 = sin(t134);
t127 = cos(t134);
t141 = sin(qJ(3));
t143 = cos(qJ(3));
t223 = Icges(4,5) * t141 + Icges(5,5) * t125 + Icges(4,6) * t143 + Icges(5,6) * t127;
t222 = Icges(4,5) * t143 / 0.2e1 - Icges(4,6) * t141 / 0.2e1;
t221 = -Icges(5,6) * t125 / 0.2e1 + t222;
t142 = sin(qJ(1));
t219 = -t142 / 0.2e1;
t144 = cos(qJ(1));
t204 = t144 / 0.2e1;
t218 = rSges(6,3) + qJ(5);
t189 = t125 * t144;
t113 = pkin(4) * t189;
t138 = cos(pkin(10));
t179 = t142 * t138;
t137 = sin(pkin(10));
t185 = t137 * t144;
t88 = t125 * t185 + t179;
t180 = t142 * t137;
t184 = t138 * t144;
t89 = -t125 * t184 + t180;
t164 = -t89 * rSges(6,1) - t88 * rSges(6,2);
t165 = t127 * t218;
t129 = t144 * qJ(2);
t139 = -qJ(4) - pkin(7);
t177 = t144 * t141 * pkin(3) + t142 * t139;
t167 = t129 + t177;
t24 = -t142 * pkin(1) - t144 * t165 + t113 + t164 + t167;
t183 = t141 * t142;
t118 = pkin(3) * t183;
t176 = t144 * pkin(1) + t142 * qJ(2);
t145 = -t139 * t144 + t118 + t176;
t190 = t125 * t142;
t86 = -t125 * t180 + t184;
t87 = t125 * t179 + t185;
t210 = t87 * rSges(6,1) + t86 * rSges(6,2) + pkin(4) * t190;
t25 = -t142 * t165 + t145 + t210;
t217 = m(6) * (t142 * t24 - t144 * t25);
t140 = -pkin(8) - qJ(5);
t133 = pkin(10) + qJ(6);
t124 = sin(t133);
t126 = cos(t133);
t181 = t142 * t126;
t72 = t124 * t189 + t181;
t182 = t142 * t124;
t188 = t126 * t144;
t73 = -t125 * t188 + t182;
t163 = -t73 * rSges(7,1) - t72 * rSges(7,2);
t121 = pkin(5) * t138 + pkin(4);
t191 = t121 * t125;
t17 = (-pkin(5) * t137 - pkin(1)) * t142 + (t191 + (-rSges(7,3) + t140) * t127) * t144 + t163 + t167;
t115 = pkin(5) * t185;
t187 = t127 * t142;
t70 = -t125 * t182 + t188;
t71 = t124 * t144 + t125 * t181;
t32 = t71 * rSges(7,1) + t70 * rSges(7,2) - rSges(7,3) * t187;
t213 = t121 * t190 + t140 * t187 + t115 + t32;
t18 = t145 + t213;
t216 = m(7) * (t142 * t17 - t144 * t18);
t215 = t223 * t142 + t224 * t144;
t214 = t224 * t142 - t223 * t144;
t212 = (rSges(4,1) * t141 + rSges(4,2) * t143) * t144;
t211 = -qJ(5) - t140;
t135 = t142 ^ 2;
t136 = t144 ^ 2;
t55 = Icges(6,6) * t125 + (Icges(6,4) * t138 - Icges(6,2) * t137) * t127;
t209 = t55 / 0.2e1;
t56 = Icges(6,5) * t125 + (Icges(6,1) * t138 - Icges(6,4) * t137) * t127;
t208 = t56 / 0.2e1;
t206 = t142 / 0.2e1;
t108 = rSges(4,1) * t143 - rSges(4,2) * t141;
t203 = m(4) * t108;
t202 = pkin(3) * t143;
t50 = Icges(7,3) * t125 + (Icges(7,5) * t126 - Icges(7,6) * t124) * t127;
t52 = Icges(7,5) * t125 + (Icges(7,1) * t126 - Icges(7,4) * t124) * t127;
t201 = t127 * t126 * t52 + t125 * t50;
t53 = rSges(7,3) * t125 + (rSges(7,1) * t126 - rSges(7,2) * t124) * t127;
t200 = (-pkin(4) + t121) * t127 + t211 * t125 + t53;
t186 = t127 * t144;
t74 = t144 * (-t142 * pkin(7) - t177);
t199 = t144 * (qJ(5) * t186 - t113) + t74;
t51 = Icges(7,6) * t125 + (Icges(7,4) * t126 - Icges(7,2) * t124) * t127;
t197 = t124 * t51;
t130 = t144 * rSges(5,3);
t178 = t142 * t143;
t119 = pkin(3) * t178;
t97 = pkin(4) * t127 + qJ(5) * t125;
t196 = t142 * t97 + t119;
t193 = Icges(5,4) * t125;
t192 = Icges(5,4) * t127;
t114 = t135 + t136;
t175 = m(6) / 0.2e1 + m(7) / 0.2e1;
t26 = Icges(7,5) * t71 + Icges(7,6) * t70 - Icges(7,3) * t187;
t28 = Icges(7,4) * t71 + Icges(7,2) * t70 - Icges(7,6) * t187;
t30 = Icges(7,1) * t71 + Icges(7,4) * t70 - Icges(7,5) * t187;
t10 = t125 * t26 + (-t124 * t28 + t126 * t30) * t127;
t13 = -t187 * t50 + t51 * t70 + t52 * t71;
t174 = -t10 / 0.2e1 - t13 / 0.2e1;
t27 = Icges(7,5) * t73 + Icges(7,6) * t72 + Icges(7,3) * t186;
t29 = Icges(7,4) * t73 + Icges(7,2) * t72 + Icges(7,6) * t186;
t31 = Icges(7,1) * t73 + Icges(7,4) * t72 + Icges(7,5) * t186;
t11 = t125 * t27 + (-t124 * t29 + t126 * t31) * t127;
t14 = t186 * t50 + t72 * t51 + t73 * t52;
t173 = t11 / 0.2e1 + t14 / 0.2e1;
t172 = (m(5) + m(6) + m(7)) * t114;
t170 = Icges(5,5) * t127 / 0.2e1 + t221;
t169 = -rSges(5,1) * t190 - rSges(5,2) * t187 - t130;
t168 = rSges(4,1) * t183 + rSges(4,2) * t178 + t144 * rSges(4,3);
t166 = -t97 - t202;
t161 = rSges(5,1) * t125 + rSges(5,2) * t127;
t19 = t125 * t32 + t187 * t53;
t33 = rSges(7,3) * t186 - t163;
t20 = -t125 * t33 + t186 * t53;
t155 = t20 * t142 - t144 * t19;
t22 = t142 * t200 + t196;
t23 = (t166 - t200) * t144;
t154 = t22 * t142 - t144 * t23;
t57 = rSges(6,3) * t125 + (rSges(6,1) * t138 - rSges(6,2) * t137) * t127;
t34 = t142 * t57 + t196;
t35 = (t166 - t57) * t144;
t152 = t34 * t142 - t144 * t35;
t150 = Icges(5,1) * t125 + t192;
t148 = Icges(5,2) * t127 + t193;
t109 = rSges(2,1) * t144 - t142 * rSges(2,2);
t107 = -t142 * rSges(2,1) - rSges(2,2) * t144;
t98 = rSges(5,1) * t127 - rSges(5,2) * t125;
t90 = t118 + (-pkin(7) - t139) * t144;
t85 = -rSges(3,2) * t144 + t142 * rSges(3,3) + t176;
t84 = rSges(3,3) * t144 + t129 + (rSges(3,2) - pkin(1)) * t142;
t59 = (-t98 - t202) * t144;
t58 = t142 * t98 + t119;
t49 = pkin(7) * t144 + t168 + t176;
t48 = t129 + t212 + (-rSges(4,3) - pkin(1) - pkin(7)) * t142;
t44 = t145 - t169;
t43 = t161 * t144 + (-rSges(5,3) - pkin(1)) * t142 + t167;
t42 = -t142 * t168 + (t142 * rSges(4,3) - t212) * t144;
t41 = Icges(6,1) * t89 + Icges(6,4) * t88 + Icges(6,5) * t186;
t40 = Icges(6,1) * t87 + Icges(6,4) * t86 - Icges(6,5) * t187;
t39 = Icges(6,4) * t89 + Icges(6,2) * t88 + Icges(6,6) * t186;
t38 = Icges(6,4) * t87 + Icges(6,2) * t86 - Icges(6,6) * t187;
t37 = Icges(6,5) * t89 + Icges(6,6) * t88 + Icges(6,3) * t186;
t36 = Icges(6,5) * t87 + Icges(6,6) * t86 - Icges(6,3) * t187;
t21 = t74 - t161 * t136 + (t169 - t90 + t130) * t142;
t16 = (-t127 * t197 + t201) * t125;
t15 = (-t142 * t33 - t144 * t32) * t127;
t12 = t144 * (rSges(6,3) * t186 - t164) + (t218 * t187 - t210 - t90) * t142 + t199;
t9 = t186 * t27 + t72 * t29 + t73 * t31;
t8 = t186 * t26 + t72 * t28 + t73 * t30;
t7 = -t187 * t27 + t29 * t70 + t31 * t71;
t6 = -t187 * t26 + t28 * t70 + t30 * t71;
t5 = (t113 + t33 + (t211 * t127 - t191) * t144) * t144 + (-t90 + t115 - t213) * t142 + t199;
t4 = t9 * t142 + t144 * t8;
t3 = t7 * t142 + t144 * t6;
t2 = t14 * t125 + (-t142 * t8 + t144 * t9) * t127;
t1 = t13 * t125 + (-t142 * t6 + t144 * t7) * t127;
t45 = [Icges(4,1) * t143 ^ 2 + Icges(3,1) + Icges(2,3) + (-t192 + (Icges(6,5) * t138 - Icges(6,6) * t137) * t127 + (Icges(5,2) + Icges(6,3)) * t125) * t125 + (Icges(5,1) * t127 - t137 * t55 + t138 * t56 - t193 - t197) * t127 + m(7) * (t17 ^ 2 + t18 ^ 2) + m(6) * (t24 ^ 2 + t25 ^ 2) + m(5) * (t43 ^ 2 + t44 ^ 2) + m(4) * (t48 ^ 2 + t49 ^ 2) + m(3) * (t84 ^ 2 + t85 ^ 2) + m(2) * (t107 ^ 2 + t109 ^ 2) + t201 + (-0.2e1 * Icges(4,4) * t143 + Icges(4,2) * t141) * t141; t216 + t217 + m(5) * (t142 * t43 - t144 * t44) + m(4) * (t142 * t48 - t144 * t49) + m(3) * (t142 * t84 - t144 * t85); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * t114 + t172; m(7) * (t17 * t22 + t18 * t23) + m(6) * (t24 * t34 + t25 * t35) + m(5) * (t43 * t58 + t44 * t59) + ((-t137 * t39 + t138 * t41) * t206 + (-t137 * t38 + t138 * t40) * t204) * t127 + (-t49 * t203 + t87 * t208 + t86 * t209 - t174 + (t36 / 0.2e1 + t148 * t219) * t125 + (Icges(5,5) * t204 - t150 * t206) * t127 + (t170 + t221) * t144) * t144 + (t48 * t203 + t89 * t208 + t88 * t209 + t173 + (t37 / 0.2e1 + Icges(5,6) * t219 + t148 * t204) * t125 + (Icges(5,5) * t206 + t150 * t204) * t127 + (t170 + t222) * t142) * t142; m(5) * (t58 * t142 - t144 * t59) + m(6) * t152 + m(7) * t154 + t114 * t203; m(7) * (t22 ^ 2 + t23 ^ 2 + t5 ^ 2) + m(6) * (t12 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(5) * (t21 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(4) * (t108 ^ 2 * t114 + t42 ^ 2) + (t3 + (-t36 * t187 + t86 * t38 + t87 * t40) * t144 + t215 * t136) * t144 + (t4 + (t37 * t186 + t88 * t39 + t89 * t41) * t142 + t214 * t135 + (t215 * t142 + t214 * t144 + t36 * t186 - t37 * t187 + t88 * t38 + t39 * t86 + t89 * t40 + t41 * t87) * t144) * t142; m(7) * (t142 * t18 + t144 * t17) + m(6) * (t142 * t25 + t144 * t24) + m(5) * (t142 * t44 + t144 * t43); 0; m(7) * (t142 * t23 + t144 * t22) + m(6) * (t142 * t35 + t144 * t34) + m(5) * (t142 * t59 + t144 * t58); t172; 0.2e1 * (-t216 / 0.2e1 - t217 / 0.2e1) * t127; -0.2e1 * t175 * t114 * t127; m(7) * (t125 * t5 - t127 * t154) + m(6) * (t125 * t12 - t127 * t152); 0; 0.2e1 * t175 * (t114 * t127 ^ 2 + t125 ^ 2); m(7) * (t17 * t20 + t18 * t19) + t16 + (t142 * t174 + t144 * t173) * t127; m(7) * t155; t1 * t204 + t2 * t206 + t125 * (t10 * t144 + t11 * t142) / 0.2e1 + m(7) * (t15 * t5 + t19 * t23 + t20 * t22) + (t4 * t204 + t3 * t219) * t127; m(7) * (t19 * t142 + t144 * t20); m(7) * (t15 * t125 - t127 * t155); t125 * t16 + m(7) * (t15 ^ 2 + t19 ^ 2 + t20 ^ 2) + (-t142 * t1 + t144 * t2 + t125 * (-t10 * t142 + t11 * t144)) * t127;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t45(1) t45(2) t45(4) t45(7) t45(11) t45(16); t45(2) t45(3) t45(5) t45(8) t45(12) t45(17); t45(4) t45(5) t45(6) t45(9) t45(13) t45(18); t45(7) t45(8) t45(9) t45(10) t45(14) t45(19); t45(11) t45(12) t45(13) t45(14) t45(15) t45(20); t45(16) t45(17) t45(18) t45(19) t45(20) t45(21);];
Mq  = res;
