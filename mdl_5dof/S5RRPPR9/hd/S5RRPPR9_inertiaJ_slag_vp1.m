% Calculate joint inertia matrix for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR9_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR9_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:42
% EndTime: 2019-12-31 19:40:46
% DurationCPUTime: 1.76s
% Computational Cost: add. (1336->253), mult. (3155->374), div. (0->0), fcn. (3195->6), ass. (0->125)
t198 = Icges(3,1) + Icges(4,1);
t197 = Icges(5,1) + Icges(4,3);
t194 = Icges(3,6) + Icges(5,5);
t191 = Icges(3,5) + Icges(4,4) + Icges(5,6);
t118 = cos(qJ(2));
t196 = (Icges(5,4) - Icges(4,5)) * t118;
t115 = sin(qJ(2));
t195 = (Icges(3,4) - Icges(4,5)) * t115;
t193 = t115 * t197 - t196;
t192 = t118 * t198 - t195;
t190 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t189 = -t191 * t118 + (-Icges(4,6) + t194) * t115;
t116 = sin(qJ(1));
t188 = -t116 / 0.2e1;
t175 = t116 / 0.2e1;
t119 = cos(qJ(1));
t187 = -t119 / 0.2e1;
t186 = t119 / 0.2e1;
t185 = t115 / 0.2e1;
t154 = t118 * t119;
t160 = t115 * t119;
t114 = sin(qJ(5));
t117 = cos(qJ(5));
t157 = t116 * t117;
t76 = -t114 * t160 - t157;
t155 = t117 * t119;
t158 = t116 * t114;
t77 = t115 * t155 - t158;
t34 = rSges(6,1) * t77 + rSges(6,2) * t76 + rSges(6,3) * t154;
t183 = pkin(4) * t160 + pkin(7) * t154 + t34;
t182 = t116 * t189 + t119 * t190;
t181 = t116 * t190 - t119 * t189;
t112 = t116 ^ 2;
t113 = t119 ^ 2;
t180 = m(4) / 0.2e1;
t179 = m(5) / 0.2e1;
t178 = m(6) / 0.2e1;
t177 = -pkin(2) - pkin(3);
t90 = rSges(3,1) * t115 + rSges(3,2) * t118;
t176 = m(3) * t90;
t45 = Icges(6,3) * t115 + (-Icges(6,5) * t117 + Icges(6,6) * t114) * t118;
t52 = Icges(6,6) * t115 + (-Icges(6,4) * t117 + Icges(6,2) * t114) * t118;
t174 = t114 * t118 * t52 + t115 * t45;
t167 = pkin(2) * t154 + qJ(3) * t160;
t173 = t112 * (pkin(2) * t118 + qJ(3) * t115) + t119 * t167;
t88 = pkin(2) * t115 - qJ(3) * t118;
t172 = -rSges(4,1) * t115 + rSges(4,3) * t118 - t88;
t59 = Icges(6,5) * t115 + (-Icges(6,1) * t117 + Icges(6,4) * t114) * t118;
t171 = t117 * t59;
t170 = t119 * rSges(4,2);
t169 = t119 * rSges(3,3);
t168 = -rSges(5,3) - qJ(4);
t165 = Icges(3,4) * t118;
t164 = Icges(5,4) * t115;
t159 = t116 * qJ(4);
t156 = t116 * t118;
t153 = t119 * qJ(4);
t151 = pkin(1) * t119 + pkin(6) * t116;
t150 = t112 + t113;
t149 = t179 + t178;
t28 = Icges(6,5) * t77 + Icges(6,6) * t76 + Icges(6,3) * t154;
t30 = Icges(6,4) * t77 + Icges(6,2) * t76 + Icges(6,6) * t154;
t32 = Icges(6,1) * t77 + Icges(6,4) * t76 + Icges(6,5) * t154;
t11 = t115 * t28 + (t114 * t30 - t117 * t32) * t118;
t14 = t154 * t45 + t52 * t76 + t59 * t77;
t148 = t11 / 0.2e1 + t14 / 0.2e1;
t74 = -t115 * t158 + t155;
t75 = t114 * t119 + t115 * t157;
t27 = Icges(6,5) * t75 + Icges(6,6) * t74 + Icges(6,3) * t156;
t29 = Icges(6,4) * t75 + Icges(6,2) * t74 + Icges(6,6) * t156;
t31 = Icges(6,1) * t75 + Icges(6,4) * t74 + Icges(6,5) * t156;
t10 = t115 * t27 + (t114 * t29 - t117 * t31) * t118;
t13 = t156 * t45 + t52 * t74 + t59 * t75;
t147 = t13 / 0.2e1 + t10 / 0.2e1;
t146 = rSges(4,1) * t154 + rSges(4,2) * t116 + rSges(4,3) * t160;
t145 = -pkin(3) * t115 - t88;
t104 = pkin(3) * t154;
t144 = t116 * (pkin(3) * t156 + t153) + t119 * (t104 - t159) + t173;
t143 = t151 + t167;
t142 = rSges(5,1) * t118 + rSges(5,2) * t115 + t145;
t141 = t191 * t185 + (-Icges(4,6) / 0.2e1 + t194 / 0.2e1) * t118;
t140 = rSges(5,1) * t160 - rSges(5,2) * t154;
t139 = -t75 * rSges(6,1) - t74 * rSges(6,2);
t138 = t104 + t143;
t137 = rSges(3,1) * t118 - rSges(3,2) * t115;
t66 = t115 * rSges(6,3) + (-rSges(6,1) * t117 + rSges(6,2) * t114) * t118;
t130 = pkin(4) * t118 - pkin(7) * t115 + t145 - t66;
t126 = -Icges(3,2) * t115 + t165;
t124 = -Icges(5,2) * t118 + t164;
t120 = rSges(3,1) * t154 - rSges(3,2) * t160 + rSges(3,3) * t116;
t109 = t119 * pkin(6);
t93 = rSges(2,1) * t119 - rSges(2,2) * t116;
t91 = -rSges(2,1) * t116 - rSges(2,2) * t119;
t42 = t172 * t119;
t41 = t172 * t116;
t40 = t120 + t151;
t39 = t169 + t109 + (-pkin(1) - t137) * t116;
t38 = t142 * t119;
t37 = t142 * t116;
t36 = t143 + t146;
t35 = t170 + t109 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t118 + (-rSges(4,3) - qJ(3)) * t115) * t116;
t33 = rSges(6,3) * t156 - t139;
t26 = t119 * t120 + (t116 * t137 - t169) * t116;
t25 = t116 * t168 + t138 + t140;
t24 = t109 + t168 * t119 + (-pkin(1) + (-rSges(5,1) - qJ(3)) * t115 + (rSges(5,2) + t177) * t118) * t116;
t23 = t130 * t119;
t22 = t130 * t116;
t21 = t115 * t34 - t154 * t66;
t20 = -t115 * t33 + t156 * t66;
t19 = (-t118 * t171 + t174) * t115;
t18 = t119 * t146 + (-t170 + (rSges(4,1) * t118 + rSges(4,3) * t115) * t116) * t116 + t173;
t17 = t138 - t159 + t183;
t16 = -t153 + t109 + (-pkin(1) + (-pkin(4) - qJ(3)) * t115 + (-rSges(6,3) - pkin(7) + t177) * t118) * t116 + t139;
t15 = (-t116 * t34 + t119 * t33) * t118;
t12 = t119 * t140 + (rSges(5,1) * t115 - rSges(5,2) * t118) * t112 + t144;
t9 = t154 * t28 + t30 * t76 + t32 * t77;
t8 = t154 * t27 + t29 * t76 + t31 * t77;
t7 = t156 * t28 + t30 * t74 + t32 * t75;
t6 = t156 * t27 + t29 * t74 + t31 * t75;
t5 = t183 * t119 + (t33 + (pkin(4) * t115 + pkin(7) * t118) * t116) * t116 + t144;
t4 = t116 * t9 - t119 * t8;
t3 = t116 * t7 - t119 * t6;
t2 = t14 * t115 + (t116 * t8 + t119 * t9) * t118;
t1 = t13 * t115 + (t116 * t6 + t119 * t7) * t118;
t43 = [Icges(2,3) + m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t24 ^ 2 + t25 ^ 2) + m(4) * (t35 ^ 2 + t36 ^ 2) + m(3) * (t39 ^ 2 + t40 ^ 2) + m(2) * (t91 ^ 2 + t93 ^ 2) + t174 + (t164 - t171 + (Icges(3,2) + t197) * t118 + t195) * t118 + (t165 + (Icges(5,2) + t198) * t115 + t196) * t115; m(6) * (t16 * t23 + t17 * t22) + m(5) * (t24 * t38 + t25 * t37) + m(4) * (t35 * t42 + t36 * t41) + (-t39 * t176 + t141 * t119 + (Icges(4,6) * t187 + t126 * t188 + t175 * t193 + t186 * t194) * t118 - t147) * t119 + (-t40 * t176 + t141 * t116 + (Icges(4,6) * t188 + t126 * t186 + t175 * t194 + t187 * t193) * t118 + t148) * t116 + (t124 * t187 * t116 + t192 * t119 * t188 + (t116 * t191 + t119 * t124) * t175 + (t116 * t192 + t119 * t191) * t186) * t115; m(6) * (t22 ^ 2 + t23 ^ 2 + t5 ^ 2) + m(5) * (t12 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(4) * (t18 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(3) * (t150 * t90 ^ 2 + t26 ^ 2) + (t113 * t182 - t3) * t119 + (t4 + t181 * t112 + (t116 * t182 + t119 * t181) * t119) * t116; 0.2e1 * ((t116 * t17 + t119 * t16) * t178 + (t116 * t25 + t119 * t24) * t179 + (t116 * t36 + t119 * t35) * t180) * t115; m(6) * (-t118 * t5 + (t116 * t22 + t119 * t23) * t115) + m(5) * (-t118 * t12 + (t116 * t37 + t119 * t38) * t115) + m(4) * (-t118 * t18 + (t116 * t41 + t119 * t42) * t115); 0.2e1 * (t180 + t149) * (t115 ^ 2 * t150 + t118 ^ 2); m(6) * (-t116 * t16 + t119 * t17) + m(5) * (-t116 * t24 + t119 * t25); m(6) * (-t116 * t23 + t119 * t22) + m(5) * (-t116 * t38 + t119 * t37); 0; 0.2e1 * t149 * t150; m(6) * (t16 * t20 + t17 * t21) + t19 + (t116 * t147 + t119 * t148) * t118; (-t10 * t119 + t11 * t116) * t185 + m(6) * (t15 * t5 + t20 * t23 + t21 * t22) + t1 * t187 + t2 * t175 + (t3 * t175 + t186 * t4) * t118; m(6) * (-t15 * t118 + (t116 * t21 + t119 * t20) * t115); m(6) * (-t116 * t20 + t119 * t21); m(6) * (t15 ^ 2 + t20 ^ 2 + t21 ^ 2) + t115 * t19 + (t119 * t2 + t116 * t1 + t115 * (t10 * t116 + t11 * t119)) * t118;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t43(1), t43(2), t43(4), t43(7), t43(11); t43(2), t43(3), t43(5), t43(8), t43(12); t43(4), t43(5), t43(6), t43(9), t43(13); t43(7), t43(8), t43(9), t43(10), t43(14); t43(11), t43(12), t43(13), t43(14), t43(15);];
Mq = res;
