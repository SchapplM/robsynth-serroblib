% Calculate joint inertia matrix for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:32
% EndTime: 2019-12-05 16:43:36
% DurationCPUTime: 1.30s
% Computational Cost: add. (2239->163), mult. (1700->229), div. (0->0), fcn. (1522->6), ass. (0->91)
t169 = Icges(5,4) + Icges(6,4);
t168 = Icges(5,1) + Icges(6,1);
t167 = Icges(5,2) + Icges(6,2);
t92 = qJ(3) + qJ(4);
t88 = cos(t92);
t166 = t169 * t88;
t87 = sin(t92);
t165 = t169 * t87;
t164 = Icges(5,5) + Icges(6,5);
t163 = Icges(5,6) + Icges(6,6);
t162 = -t167 * t87 + t166;
t161 = t168 * t88 - t165;
t91 = pkin(8) + qJ(2);
t85 = sin(t91);
t86 = cos(t91);
t160 = t162 * t85 - t163 * t86;
t159 = t162 * t86 + t163 * t85;
t158 = t161 * t85 - t164 * t86;
t157 = t161 * t86 + t164 * t85;
t156 = Icges(5,3) + Icges(6,3);
t155 = -t163 * t87 + t164 * t88;
t154 = t167 * t88 + t165;
t153 = t168 * t87 + t166;
t152 = t155 * t86 + t156 * t85;
t151 = -t155 * t85 + t156 * t86;
t150 = -t157 * t88 + t159 * t87;
t149 = -t158 * t88 + t160 * t87;
t148 = t85 * pkin(6);
t133 = t86 * t88;
t134 = t86 * t87;
t94 = cos(qJ(3));
t82 = t94 * pkin(3) + pkin(2);
t68 = pkin(4) * t88 + t82;
t147 = rSges(6,1) * t133 - rSges(6,2) * t134 + t85 * rSges(6,3) + t86 * t68;
t146 = t163 * t88 + t164 * t87;
t83 = t85 ^ 2;
t84 = t86 ^ 2;
t130 = t83 + t84;
t145 = -rSges(6,1) * t88 + rSges(6,2) * t87 - t68;
t144 = t153 * t88 - t154 * t87;
t143 = t151 * t84 + (t150 * t85 + (-t149 + t152) * t86) * t85;
t93 = sin(qJ(3));
t102 = Icges(4,5) * t94 - Icges(4,6) * t93;
t142 = t85 * (Icges(4,3) * t85 + t102 * t86);
t141 = t85 / 0.2e1;
t140 = -t86 / 0.2e1;
t95 = -pkin(7) - pkin(6);
t139 = (t152 * t83 + ((-t150 + t151) * t85 + t149 * t86) * t86) * t85;
t138 = pkin(3) * t93;
t137 = rSges(4,1) * t94;
t136 = rSges(4,2) * t93;
t135 = t86 * rSges(4,3);
t69 = t86 * t82;
t81 = t86 * pkin(6);
t132 = t85 * (t81 + (-pkin(2) + t82) * t85) + t86 * (-t86 * pkin(2) - t148 + t69);
t119 = rSges(5,1) * t88 - rSges(5,2) * t87;
t99 = rSges(5,1) * t133 - rSges(5,2) * t134 + t85 * rSges(5,3);
t16 = t85 * (-t86 * rSges(5,3) + t119 * t85) + t86 * t99;
t131 = t85 * rSges(4,3) + t86 * t137;
t128 = Icges(4,4) * t93;
t127 = Icges(4,4) * t94;
t67 = t87 * rSges(5,1) + t88 * rSges(5,2);
t122 = -t67 - t138;
t121 = -t88 * rSges(6,2) + (-rSges(6,1) - pkin(4)) * t87;
t6 = (-t69 + t147) * t86 + (-t86 * rSges(6,3) + (-t82 - t145) * t85) * t85;
t120 = -t136 + t137;
t109 = t143 * t86 + t139;
t108 = Icges(4,1) * t94 - t128;
t105 = -Icges(4,2) * t93 + t127;
t97 = t121 - t138;
t96 = (t144 * t86 + t146 * t85 + t157 * t87 + t159 * t88) * t141 + (t144 * t85 - t146 * t86 + t158 * t87 + t160 * t88) * t140;
t90 = -qJ(5) + t95;
t76 = t93 * rSges(4,1) + t94 * rSges(4,2);
t58 = t86 * rSges(3,1) - t85 * rSges(3,2);
t57 = -t85 * rSges(3,1) - t86 * rSges(3,2);
t49 = t122 * t86;
t48 = t122 * t85;
t35 = t121 * t86;
t34 = t121 * t85;
t27 = t148 + (pkin(2) - t136) * t86 + t131;
t26 = t135 + t81 + (-pkin(2) - t120) * t85;
t25 = t97 * t86;
t24 = t97 * t85;
t23 = -t85 * t95 + t69 + t99;
t22 = (rSges(5,3) - t95) * t86 + (-t119 - t82) * t85;
t19 = -t85 * t90 + t147;
t18 = (rSges(6,3) - t90) * t86 + t145 * t85;
t17 = t86 * (-t86 * t136 + t131) + (t120 * t85 - t135) * t85;
t7 = t16 + t132;
t1 = t6 + t132;
t2 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; t94 * (Icges(4,2) * t94 + t128) + t93 * (Icges(4,1) * t93 + t127) + Icges(3,3) + t154 * t88 + t153 * t87 + m(6) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t26 ^ 2 + t27 ^ 2) + m(3) * (t57 ^ 2 + t58 ^ 2); m(4) * t17 + m(5) * t7 + m(6) * t1; (t94 * (Icges(4,6) * t85 + t105 * t86) + t93 * (Icges(4,5) * t85 + t108 * t86)) * t141 + (t94 * (-Icges(4,6) * t86 + t105 * t85) + t93 * (-Icges(4,5) * t86 + t108 * t85)) * t140 + m(6) * (t25 * t18 + t24 * t19) + m(5) * (t49 * t22 + t48 * t23) + m(4) * (-t26 * t86 - t27 * t85) * t76 + (t83 / 0.2e1 + t84 / 0.2e1) * (Icges(4,5) * t93 + Icges(4,6) * t94) + t96; m(6) * (t1 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t48 ^ 2 + t49 ^ 2 + t7 ^ 2) + m(4) * (t130 * t76 ^ 2 + t17 ^ 2) + t83 * t142 + t139 + (t86 * t142 - t130 * (-Icges(4,3) * t86 + t102 * t85) + t143) * t86; m(5) * t16 + m(6) * t6; m(6) * (t35 * t18 + t34 * t19) + m(5) * (-t22 * t86 - t23 * t85) * t67 + t96; m(6) * (t6 * t1 + t34 * t24 + t35 * t25) + m(5) * (t16 * t7 + (-t48 * t85 - t49 * t86) * t67) + t109; m(5) * (t130 * t67 ^ 2 + t16 ^ 2) + m(6) * (t34 ^ 2 + t35 ^ 2 + t6 ^ 2) + t109; 0; m(6) * (t85 * t18 - t86 * t19); m(6) * (-t86 * t24 + t85 * t25); m(6) * (-t86 * t34 + t85 * t35); m(6) * t130;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
