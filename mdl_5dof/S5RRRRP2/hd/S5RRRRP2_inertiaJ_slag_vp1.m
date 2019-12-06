% Calculate joint inertia matrix for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:47:33
% EndTime: 2019-12-05 18:47:38
% DurationCPUTime: 1.50s
% Computational Cost: add. (2961->198), mult. (2286->271), div. (0->0), fcn. (2028->8), ass. (0->116)
t201 = Icges(5,4) + Icges(6,4);
t200 = Icges(5,1) + Icges(6,1);
t199 = Icges(5,5) + Icges(6,5);
t198 = Icges(5,2) + Icges(6,2);
t197 = Icges(5,6) + Icges(6,6);
t118 = qJ(3) + qJ(4);
t114 = cos(t118);
t196 = t201 * t114;
t112 = sin(t118);
t195 = t201 * t112;
t192 = t112 * t200 + t196;
t193 = t114 * t198 + t195;
t194 = (t192 + t196) * t114 + (-t193 - t195 + (-t198 + t200) * t114) * t112;
t191 = Icges(5,3) + Icges(6,3);
t190 = -t112 * t197 + t114 * t199;
t187 = 0.2e1 * t112 * t199 + 0.2e1 * t114 * t197;
t119 = qJ(1) + qJ(2);
t113 = sin(t119);
t115 = cos(t119);
t186 = -t113 * t190 + t115 * t191;
t185 = t113 * t191 + t115 * t190;
t124 = -pkin(8) - pkin(7);
t117 = -qJ(5) + t124;
t158 = t112 * t115;
t183 = rSges(6,2) * t158 + (-rSges(6,3) + t117) * t113;
t181 = t113 * t115;
t110 = t113 ^ 2;
t111 = t115 ^ 2;
t120 = sin(qJ(3));
t122 = cos(qJ(3));
t96 = rSges(4,1) * t120 + rSges(4,2) * t122;
t180 = m(4) * t96;
t81 = rSges(5,1) * t112 + rSges(5,2) * t114;
t179 = m(5) * t81;
t178 = t113 / 0.2e1;
t177 = t115 / 0.2e1;
t176 = pkin(3) * t120;
t121 = sin(qJ(1));
t175 = t121 * pkin(1);
t123 = cos(qJ(1));
t174 = t123 * pkin(1);
t109 = pkin(3) * t122 + pkin(2);
t173 = pkin(2) - t109;
t103 = t113 * t124;
t155 = t114 * t115;
t86 = pkin(4) * t114 + t109;
t166 = t109 - t86;
t172 = (rSges(6,1) * t155 - t115 * t166 + t103 - t183) * t115;
t157 = t113 * t114;
t159 = t112 * t113;
t169 = rSges(6,2) * t159 + rSges(6,3) * t115;
t171 = -(-t117 + t124) * t115 - t166 * t113 + rSges(6,1) * t157 - t169;
t80 = rSges(6,1) * t112 + rSges(6,2) * t114;
t42 = pkin(4) * t159 + t113 * t80;
t170 = rSges(4,1) * t122;
t168 = rSges(5,2) * t159 + rSges(5,3) * t115;
t156 = t113 * t120;
t167 = rSges(4,2) * t156 + rSges(4,3) * t115;
t165 = Icges(4,4) * t120;
t164 = Icges(4,4) * t122;
t154 = t115 * t124;
t153 = t111 + t110;
t152 = t185 * t113 * t110 + (t186 * t111 + (t113 * t186 + t115 * t185) * t113) * t115;
t151 = -pkin(2) - t170;
t150 = -pkin(4) * t112 - t80;
t149 = -rSges(6,1) * t114 - t86;
t148 = rSges(5,2) * t158 - t113 * rSges(5,3);
t83 = -rSges(3,1) * t115 + rSges(3,2) * t113;
t146 = -rSges(5,1) * t114 - t109;
t82 = -rSges(3,1) * t113 - rSges(3,2) * t115;
t94 = Icges(4,2) * t122 + t165;
t95 = Icges(4,1) * t120 + t164;
t137 = t120 * t94 - t122 * t95;
t136 = Icges(4,1) * t122 - t165;
t133 = -Icges(4,2) * t120 + t164;
t130 = Icges(4,5) * t122 - Icges(4,6) * t120;
t127 = t112 * t192 + t114 * t193 + t120 * t95 + t122 * t94 + Icges(3,3);
t126 = (t187 * t177 + t178 * t194) * t115 + (-t177 * t194 + t187 * t178) * t113;
t108 = t115 * pkin(7);
t36 = t113 * t151 + t108 + t167;
t25 = t115 * t149 + t183;
t30 = t115 * t146 + t103 + t148;
t24 = t113 * t149 - t115 * t117 + t169;
t100 = t115 * t120 * rSges(4,2);
t37 = t100 + t151 * t115 + (-rSges(4,3) - pkin(7)) * t113;
t29 = t113 * t146 - t154 + t168;
t93 = Icges(4,5) * t120 + Icges(4,6) * t122;
t125 = t126 + (t113 * t93 - t137 * t115 + t120 * (Icges(4,5) * t113 + t115 * t136) + t122 * (Icges(4,6) * t113 + t115 * t133)) * t178 + (t137 * t113 + t115 * t93 + t120 * (Icges(4,5) * t115 - t113 * t136) + t122 * (Icges(4,6) * t115 - t113 * t133)) * t177;
t102 = pkin(3) * t156;
t98 = -rSges(2,1) * t123 + rSges(2,2) * t121;
t97 = -rSges(2,1) * t121 - rSges(2,2) * t123;
t72 = t83 - t174;
t71 = t82 - t175;
t61 = Icges(4,3) * t113 + t115 * t130;
t60 = Icges(4,3) * t115 - t113 * t130;
t59 = (-t81 - t176) * t115;
t58 = t113 * t81 + t102;
t57 = -rSges(5,1) * t157 + t168;
t43 = t150 * t115;
t41 = t113 * t173 - t108 - t154;
t40 = t115 * (rSges(5,1) * t155 - t148);
t38 = t115 * (-t113 * pkin(7) - t115 * t173 - t103);
t35 = (t150 - t176) * t115;
t34 = t102 + t42;
t33 = t37 - t174;
t32 = t36 - t175;
t27 = t30 - t174;
t26 = t29 - t175;
t21 = t25 - t174;
t20 = t24 - t175;
t19 = t115 * (t113 * rSges(4,3) + t115 * t170 - t100) - t113 * (-t113 * t170 + t167);
t16 = -t113 * t57 + t40;
t7 = t38 + t40 + (-t41 - t57) * t113;
t6 = t113 * t171 + t172;
t1 = t38 + (-t41 + t171) * t113 + t172;
t2 = [Icges(2,3) + m(2) * (t97 ^ 2 + t98 ^ 2) + m(3) * (t71 ^ 2 + t72 ^ 2) + m(4) * (t32 ^ 2 + t33 ^ 2) + m(5) * (t26 ^ 2 + t27 ^ 2) + m(6) * (t20 ^ 2 + t21 ^ 2) + t127; m(3) * (t71 * t82 + t72 * t83) + m(4) * (t32 * t36 + t33 * t37) + m(5) * (t26 * t29 + t27 * t30) + m(6) * (t20 * t24 + t21 * t25) + t127; m(6) * (t24 ^ 2 + t25 ^ 2) + m(5) * (t29 ^ 2 + t30 ^ 2) + m(4) * (t36 ^ 2 + t37 ^ 2) + m(3) * (t82 ^ 2 + t83 ^ 2) + t127; t125 + m(6) * (t20 * t35 + t21 * t34) + m(5) * (t26 * t59 + t27 * t58) + (t113 * t33 - t115 * t32) * t180; t125 + m(6) * (t24 * t35 + t25 * t34) + m(5) * (t29 * t59 + t30 * t58) + (t113 * t37 - t115 * t36) * t180; m(6) * (t1 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(5) * (t58 ^ 2 + t59 ^ 2 + t7 ^ 2) + m(4) * (t153 * t96 ^ 2 + t19 ^ 2) + t113 * (t110 * t61 + t181 * t60) + t115 * (t111 * t60 + t181 * t61) + t152; m(6) * (t20 * t43 + t21 * t42) + (t113 * t27 - t115 * t26) * t179 + t126; m(6) * (t24 * t43 + t25 * t42) + (t113 * t30 - t115 * t29) * t179 + t126; m(6) * (t1 * t6 + t34 * t42 + t35 * t43) + m(5) * (t16 * t7 + (t113 * t58 - t115 * t59) * t81) + t152; m(5) * (t153 * t81 ^ 2 + t16 ^ 2) + m(6) * (t42 ^ 2 + t43 ^ 2 + t6 ^ 2) + t152; m(6) * (t113 * t20 + t115 * t21); m(6) * (t113 * t24 + t115 * t25); m(6) * (t113 * t35 + t115 * t34); m(6) * (t113 * t43 + t115 * t42); m(6) * t153;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
