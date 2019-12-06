% Calculate joint inertia matrix for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:18
% EndTime: 2019-12-05 18:42:22
% DurationCPUTime: 1.10s
% Computational Cost: add. (2683->215), mult. (1866->300), div. (0->0), fcn. (1632->10), ass. (0->118)
t176 = Icges(4,3) + Icges(5,3);
t111 = qJ(3) + pkin(9);
t104 = sin(t111);
t105 = cos(t111);
t114 = sin(qJ(3));
t116 = cos(qJ(3));
t175 = Icges(4,5) * t116 + Icges(5,5) * t105 - Icges(4,6) * t114 - Icges(5,6) * t104;
t174 = Icges(4,5) * t114 + Icges(5,5) * t104 + Icges(4,6) * t116 + Icges(5,6) * t105;
t112 = qJ(1) + qJ(2);
t107 = sin(t112);
t108 = cos(t112);
t173 = t107 * t108;
t148 = Icges(5,4) * t104;
t69 = Icges(5,2) * t105 + t148;
t147 = Icges(5,4) * t105;
t70 = Icges(5,1) * t104 + t147;
t150 = Icges(4,4) * t114;
t84 = Icges(4,2) * t116 + t150;
t149 = Icges(4,4) * t116;
t85 = Icges(4,1) * t114 + t149;
t172 = t104 * t69 - t105 * t70 + t114 * t84 - t116 * t85;
t171 = -t175 * t107 + t176 * t108;
t170 = t176 * t107 + t175 * t108;
t102 = t107 ^ 2;
t103 = t108 ^ 2;
t106 = qJ(5) + t111;
t98 = sin(t106);
t99 = cos(t106);
t132 = Icges(6,5) * t99 - Icges(6,6) * t98;
t35 = Icges(6,3) * t108 - t132 * t107;
t36 = Icges(6,3) * t107 + t132 * t108;
t169 = t108 * (t103 * t35 + t36 * t173) + t107 * (t102 * t36 + t35 * t173);
t86 = t114 * rSges(4,1) + t116 * rSges(4,2);
t168 = m(4) * t86;
t67 = t98 * rSges(6,1) + t99 * rSges(6,2);
t167 = m(6) * t67;
t166 = t107 / 0.2e1;
t165 = t108 / 0.2e1;
t164 = rSges(6,1) * t99;
t163 = rSges(6,2) * t98;
t162 = pkin(3) * t114;
t115 = sin(qJ(1));
t161 = t115 * pkin(1);
t117 = cos(qJ(1));
t160 = t117 * pkin(1);
t101 = t116 * pkin(3) + pkin(2);
t159 = pkin(2) - t101;
t113 = -qJ(4) - pkin(7);
t158 = t108 * rSges(6,3) + t107 * t163;
t146 = t107 * t114;
t157 = rSges(4,2) * t146 + t108 * rSges(4,3);
t156 = rSges(4,1) * t116;
t155 = rSges(5,1) * t105;
t154 = rSges(5,2) * t104;
t153 = Icges(6,4) * t98;
t152 = Icges(6,4) * t99;
t78 = pkin(4) * t105 + t101;
t151 = -t101 + t78;
t145 = t108 * t113;
t144 = t103 + t102;
t133 = -Icges(6,2) * t98 + t152;
t134 = Icges(6,1) * t99 - t153;
t65 = Icges(6,2) * t99 + t153;
t66 = Icges(6,1) * t98 + t152;
t135 = t65 * t98 - t66 * t99;
t64 = Icges(6,5) * t98 + Icges(6,6) * t99;
t143 = (t107 * t64 - t135 * t108 + t99 * (Icges(6,6) * t107 + t133 * t108) + t98 * (Icges(6,5) * t107 + t134 * t108)) * t166 + (t135 * t107 + t108 * t64 + t99 * (Icges(6,6) * t108 - t133 * t107) + t98 * (Icges(6,5) * t108 - t134 * t107)) * t165;
t142 = -pkin(2) - t156;
t141 = -t78 - t164;
t140 = pkin(4) * t104 + t67;
t73 = -t108 * rSges(3,1) + t107 * rSges(3,2);
t139 = -t107 * rSges(6,3) + t108 * t163;
t138 = -t101 - t155;
t72 = -t107 * rSges(3,1) - t108 * rSges(3,2);
t125 = Icges(4,1) * t116 - t150;
t124 = Icges(5,1) * t105 - t148;
t123 = -Icges(4,2) * t114 + t149;
t122 = -Icges(5,2) * t104 + t147;
t119 = t104 * t70 + t105 * t69 + t114 * t85 + t116 * t84 + t99 * t65 + t98 * t66 + Icges(3,3);
t100 = t108 * pkin(7);
t32 = t142 * t107 + t100 + t157;
t110 = -pkin(8) + t113;
t91 = t107 * t110;
t20 = t141 * t108 + t139 + t91;
t19 = t141 * t107 - t108 * t110 + t158;
t118 = t143 + (t104 * (Icges(5,5) * t107 + t108 * t124) + t105 * (Icges(5,6) * t107 + t108 * t122) + t114 * (Icges(4,5) * t107 + t108 * t125) + t116 * (Icges(4,6) * t107 + t108 * t123) - t172 * t108 + t174 * t107) * t166 + (t104 * (Icges(5,5) * t108 - t107 * t124) + t105 * (Icges(5,6) * t108 - t107 * t122) + t114 * (Icges(4,5) * t108 - t107 * t125) + t116 * (Icges(4,6) * t108 - t107 * t123) + t174 * t108 + t172 * t107) * t165;
t81 = t108 * t154;
t93 = t107 * t113;
t26 = -t107 * rSges(5,3) + t138 * t108 + t81 + t93;
t90 = t108 * t114 * rSges(4,2);
t33 = t90 + t142 * t108 + (-rSges(4,3) - pkin(7)) * t107;
t80 = t107 * t154;
t25 = t108 * rSges(5,3) + t138 * t107 - t145 + t80;
t92 = pkin(3) * t146;
t88 = -t117 * rSges(2,1) + t115 * rSges(2,2);
t87 = -t115 * rSges(2,1) - t117 * rSges(2,2);
t71 = t104 * rSges(5,1) + t105 * rSges(5,2);
t62 = t73 - t160;
t61 = t72 - t161;
t50 = (-t71 - t162) * t108;
t49 = t107 * t71 + t92;
t42 = -t107 * t164 + t158;
t41 = t159 * t107 - t100 - t145;
t34 = t108 * (t108 * t164 - t139);
t31 = t108 * (-t107 * pkin(7) - t159 * t108 - t93);
t30 = t33 - t160;
t29 = t32 - t161;
t28 = (-t140 - t162) * t108;
t27 = t140 * t107 + t92;
t24 = t26 - t160;
t23 = t25 - t161;
t18 = t108 * (t107 * rSges(4,3) + t108 * t156 - t90) - t107 * (-t107 * t156 + t157);
t17 = t20 - t160;
t16 = t19 - t161;
t13 = -t107 * t42 + t34;
t4 = t31 + t108 * (t108 * t155 - t81) + (t107 * t155 - t41 - t80) * t107;
t3 = t31 + t34 + (t151 * t108 - t91 + t93) * t108 + (-t41 - t42 + t151 * t107 + (t110 - t113) * t108) * t107;
t1 = [Icges(2,3) + m(2) * (t87 ^ 2 + t88 ^ 2) + m(3) * (t61 ^ 2 + t62 ^ 2) + m(4) * (t29 ^ 2 + t30 ^ 2) + m(5) * (t23 ^ 2 + t24 ^ 2) + m(6) * (t16 ^ 2 + t17 ^ 2) + t119; m(3) * (t72 * t61 + t73 * t62) + m(4) * (t32 * t29 + t33 * t30) + m(5) * (t25 * t23 + t26 * t24) + m(6) * (t19 * t16 + t20 * t17) + t119; m(6) * (t19 ^ 2 + t20 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2) + m(4) * (t32 ^ 2 + t33 ^ 2) + m(3) * (t72 ^ 2 + t73 ^ 2) + t119; t118 + m(5) * (t50 * t23 + t49 * t24) + m(6) * (t28 * t16 + t27 * t17) + (t107 * t30 - t108 * t29) * t168; t118 + m(6) * (t28 * t19 + t27 * t20) + m(5) * (t50 * t25 + t49 * t26) + (t107 * t33 - t108 * t32) * t168; m(6) * (t27 ^ 2 + t28 ^ 2 + t3 ^ 2) + m(5) * (t4 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(4) * (t144 * t86 ^ 2 + t18 ^ 2) + t169 + t170 * t107 * t102 + (t171 * t103 + (t171 * t107 + t170 * t108) * t107) * t108; m(5) * (t107 * t23 + t108 * t24) + m(6) * (t107 * t16 + t108 * t17); m(6) * (t107 * t19 + t108 * t20) + m(5) * (t107 * t25 + t108 * t26); m(6) * (t107 * t28 + t108 * t27) + m(5) * (t107 * t50 + t108 * t49); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t144; (t107 * t17 - t108 * t16) * t167 + t143; (t107 * t20 - t108 * t19) * t167 + t143; m(6) * (t13 * t3 + (t107 * t27 - t108 * t28) * t67) + t169; 0; m(6) * (t144 * t67 ^ 2 + t13 ^ 2) + t169;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
