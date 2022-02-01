% Calculate joint inertia matrix for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:33
% EndTime: 2022-01-23 09:29:34
% DurationCPUTime: 1.30s
% Computational Cost: add. (2263->175), mult. (1730->238), div. (0->0), fcn. (1546->8), ass. (0->97)
t175 = Icges(5,4) + Icges(6,4);
t174 = Icges(5,1) + Icges(6,1);
t173 = Icges(5,2) + Icges(6,2);
t95 = qJ(3) + qJ(4);
t90 = cos(t95);
t172 = t175 * t90;
t89 = sin(t95);
t171 = t175 * t89;
t170 = Icges(5,5) + Icges(6,5);
t169 = Icges(5,6) + Icges(6,6);
t168 = -t173 * t89 + t172;
t167 = t174 * t90 - t171;
t94 = qJ(1) + pkin(8);
t87 = sin(t94);
t88 = cos(t94);
t166 = t168 * t87 - t169 * t88;
t165 = t168 * t88 + t169 * t87;
t164 = t167 * t87 - t170 * t88;
t163 = t167 * t88 + t170 * t87;
t162 = Icges(5,3) + Icges(6,3);
t161 = -t169 * t89 + t170 * t90;
t160 = t173 * t90 + t171;
t159 = t174 * t89 + t172;
t158 = t161 * t88 + t162 * t87;
t157 = -t161 * t87 + t162 * t88;
t156 = -t163 * t90 + t165 * t89;
t155 = -t164 * t90 + t166 * t89;
t154 = t87 * pkin(6);
t153 = t169 * t90 + t170 * t89;
t85 = t87 ^ 2;
t86 = t88 ^ 2;
t135 = t85 + t86;
t138 = t88 * t90;
t139 = t88 * t89;
t98 = cos(qJ(3));
t84 = pkin(3) * t98 + pkin(2);
t68 = pkin(4) * t90 + t84;
t152 = rSges(6,1) * t138 - rSges(6,2) * t139 + rSges(6,3) * t87 + t68 * t88;
t151 = -rSges(6,1) * t90 + rSges(6,2) * t89 - t68;
t150 = t159 * t90 - t160 * t89;
t149 = t157 * t86 + (t156 * t87 + (-t155 + t158) * t88) * t87;
t96 = sin(qJ(3));
t107 = Icges(4,5) * t98 - Icges(4,6) * t96;
t148 = t87 * (Icges(4,3) * t87 + t107 * t88);
t147 = t87 / 0.2e1;
t146 = -t88 / 0.2e1;
t100 = -pkin(7) - pkin(6);
t145 = (t158 * t85 + ((-t156 + t157) * t87 + t155 * t88) * t88) * t87;
t144 = pkin(3) * t96;
t97 = sin(qJ(1));
t143 = t97 * pkin(1);
t142 = rSges(4,1) * t98;
t141 = rSges(4,2) * t96;
t140 = t88 * rSges(4,3);
t69 = t88 * t84;
t83 = t88 * pkin(6);
t137 = t87 * (t83 + (-pkin(2) + t84) * t87) + t88 * (-t88 * pkin(2) - t154 + t69);
t104 = rSges(5,1) * t138 - rSges(5,2) * t139 + rSges(5,3) * t87;
t124 = rSges(5,1) * t90 - rSges(5,2) * t89;
t16 = t87 * (-t88 * rSges(5,3) + t124 * t87) + t88 * t104;
t136 = rSges(4,3) * t87 + t142 * t88;
t134 = Icges(4,4) * t96;
t133 = Icges(4,4) * t98;
t67 = rSges(5,1) * t89 + rSges(5,2) * t90;
t127 = -t67 - t144;
t126 = -t90 * rSges(6,2) + (-rSges(6,1) - pkin(4)) * t89;
t6 = (-t69 + t152) * t88 + (-t88 * rSges(6,3) + (-t84 - t151) * t87) * t87;
t125 = -t141 + t142;
t114 = t149 * t88 + t145;
t113 = Icges(4,1) * t98 - t134;
t110 = -Icges(4,2) * t96 + t133;
t102 = t126 - t144;
t101 = (t150 * t88 + t153 * t87 + t163 * t89 + t165 * t90) * t147 + (t150 * t87 - t153 * t88 + t164 * t89 + t166 * t90) * t146;
t99 = cos(qJ(1));
t93 = -qJ(5) + t100;
t92 = t99 * pkin(1);
t78 = rSges(2,1) * t99 - rSges(2,2) * t97;
t77 = -rSges(2,1) * t97 - rSges(2,2) * t99;
t76 = rSges(4,1) * t96 + rSges(4,2) * t98;
t57 = rSges(3,1) * t88 - rSges(3,2) * t87 + t92;
t56 = -rSges(3,1) * t87 - rSges(3,2) * t88 - t143;
t49 = t127 * t88;
t48 = t127 * t87;
t35 = t126 * t88;
t34 = t126 * t87;
t27 = t102 * t88;
t26 = t102 * t87;
t25 = t154 + t92 + (pkin(2) - t141) * t88 + t136;
t24 = t140 - t143 + t83 + (-pkin(2) - t125) * t87;
t21 = -t87 * t100 + t104 + t69 + t92;
t20 = -t143 + (rSges(5,3) - t100) * t88 + (-t124 - t84) * t87;
t19 = -t87 * t93 + t152 + t92;
t18 = -t143 + (rSges(6,3) - t93) * t88 + t151 * t87;
t17 = t88 * (-t141 * t88 + t136) + (t125 * t87 - t140) * t87;
t7 = t16 + t137;
t1 = t6 + t137;
t2 = [t98 * (Icges(4,2) * t98 + t134) + t96 * (Icges(4,1) * t96 + t133) + Icges(2,3) + Icges(3,3) + t160 * t90 + t159 * t89 + m(6) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + m(3) * (t56 ^ 2 + t57 ^ 2) + m(2) * (t77 ^ 2 + t78 ^ 2); 0; m(3) + m(4) + m(5) + m(6); (t98 * (-Icges(4,6) * t88 + t110 * t87) + t96 * (-Icges(4,5) * t88 + t113 * t87)) * t146 + (t98 * (Icges(4,6) * t87 + t110 * t88) + t96 * (Icges(4,5) * t87 + t113 * t88)) * t147 + m(6) * (t18 * t27 + t19 * t26) + m(5) * (t20 * t49 + t21 * t48) + m(4) * (-t24 * t88 - t25 * t87) * t76 + (t86 / 0.2e1 + t85 / 0.2e1) * (Icges(4,5) * t96 + Icges(4,6) * t98) + t101; m(4) * t17 + m(5) * t7 + m(6) * t1; m(6) * (t1 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(5) * (t48 ^ 2 + t49 ^ 2 + t7 ^ 2) + t85 * t148 + m(4) * (t135 * t76 ^ 2 + t17 ^ 2) + t145 + (t88 * t148 - t135 * (-Icges(4,3) * t88 + t107 * t87) + t149) * t88; m(6) * (t18 * t35 + t19 * t34) + m(5) * (-t20 * t88 - t21 * t87) * t67 + t101; m(5) * t16 + m(6) * t6; m(6) * (t1 * t6 + t26 * t34 + t27 * t35) + m(5) * (t16 * t7 + (-t48 * t87 - t49 * t88) * t67) + t114; m(5) * (t135 * t67 ^ 2 + t16 ^ 2) + m(6) * (t34 ^ 2 + t35 ^ 2 + t6 ^ 2) + t114; m(6) * (t18 * t87 - t19 * t88); 0; m(6) * (-t26 * t88 + t27 * t87); m(6) * (-t34 * t88 + t35 * t87); m(6) * t135;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
