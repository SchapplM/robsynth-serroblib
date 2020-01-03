% Calculate joint inertia matrix for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP5_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP5_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:38
% EndTime: 2019-12-31 17:16:42
% DurationCPUTime: 1.30s
% Computational Cost: add. (1283->157), mult. (1644->230), div. (0->0), fcn. (1489->6), ass. (0->93)
t174 = Icges(4,4) - Icges(5,5);
t173 = Icges(4,1) + Icges(5,1);
t172 = Icges(4,2) + Icges(5,3);
t90 = qJ(2) + qJ(3);
t83 = cos(t90);
t171 = t174 * t83;
t82 = sin(t90);
t170 = t174 * t82;
t169 = -t172 * t82 + t171;
t168 = Icges(5,4) + Icges(4,5);
t167 = Icges(4,6) - Icges(5,6);
t166 = t173 * t83 - t170;
t92 = sin(qJ(1));
t94 = cos(qJ(1));
t165 = -t167 * t94 + t169 * t92;
t164 = t167 * t92 + t169 * t94;
t163 = t166 * t92 - t168 * t94;
t162 = t166 * t94 + t168 * t92;
t156 = rSges(5,1) + pkin(3);
t161 = t156 * t83;
t160 = -t167 * t82 + t168 * t83;
t159 = Icges(5,2) + Icges(4,3);
t158 = t172 * t83 + t170;
t157 = t173 * t82 + t171;
t155 = t159 * t94 - t160 * t92;
t154 = t159 * t92 + t160 * t94;
t153 = -t162 * t83 + t164 * t82;
t152 = -t163 * t83 + t165 * t82;
t88 = t92 ^ 2;
t151 = t92 * pkin(5);
t149 = t92 * rSges(5,2);
t148 = t167 * t83 + t168 * t82;
t89 = t94 ^ 2;
t129 = t88 + t89;
t147 = rSges(5,3) + qJ(4);
t146 = t157 * t83 - t158 * t82;
t122 = qJ(4) * t82;
t135 = t83 * t94;
t136 = t82 * t94;
t145 = rSges(5,3) * t136 + t94 * t122 + t156 * t135 + t149;
t144 = t155 * t89 + (t153 * t92 + (-t152 + t154) * t94) * t92;
t91 = sin(qJ(2));
t93 = cos(qJ(2));
t100 = Icges(3,5) * t93 - Icges(3,6) * t91;
t143 = t92 * (Icges(3,3) * t92 + t100 * t94);
t142 = t92 / 0.2e1;
t141 = -t94 / 0.2e1;
t140 = (t154 * t88 + ((-t153 + t155) * t92 + t152 * t94) * t94) * t92;
t139 = pkin(2) * t91;
t138 = rSges(3,1) * t93;
t137 = rSges(3,2) * t91;
t134 = t94 * rSges(3,3);
t80 = t93 * pkin(2) + pkin(1);
t74 = t94 * t80;
t87 = t94 * pkin(5);
t133 = t92 * (t87 + (-pkin(1) + t80) * t92) + t94 * (-t94 * pkin(1) - t151 + t74);
t116 = rSges(4,1) * t83 - rSges(4,2) * t82;
t97 = rSges(4,1) * t135 - rSges(4,2) * t136 + t92 * rSges(4,3);
t18 = t92 * (-t94 * rSges(4,3) + t116 * t92) + t94 * t97;
t132 = t147 * t83 - t156 * t82;
t130 = t92 * rSges(3,3) + t94 * t138;
t128 = Icges(3,4) * t91;
t127 = Icges(3,4) * t93;
t65 = t82 * rSges(4,1) + t83 * rSges(4,2);
t120 = -t65 - t139;
t95 = -pkin(6) - pkin(5);
t119 = -t92 * t95 + t74;
t7 = (t145 - t149) * t94 + (rSges(5,3) * t82 + t122 + t161) * t88;
t118 = t132 - t139;
t117 = -t137 + t138;
t107 = t144 * t94 + t140;
t106 = Icges(3,1) * t93 - t128;
t103 = -Icges(3,2) * t91 + t127;
t96 = (t146 * t94 + t148 * t92 + t162 * t82 + t164 * t83) * t142 + (t146 * t92 - t148 * t94 + t163 * t82 + t165 * t83) * t141;
t72 = t94 * rSges(2,1) - t92 * rSges(2,2);
t71 = -t92 * rSges(2,1) - t94 * rSges(2,2);
t70 = t91 * rSges(3,1) + t93 * rSges(3,2);
t35 = t120 * t94;
t34 = t120 * t92;
t27 = t151 + (pkin(1) - t137) * t94 + t130;
t26 = t134 + t87 + (-pkin(1) - t117) * t92;
t25 = t132 * t94;
t24 = t132 * t92;
t23 = t119 + t97;
t22 = (rSges(4,3) - t95) * t94 + (-t116 - t80) * t92;
t21 = t118 * t94;
t20 = t118 * t92;
t19 = t94 * (-t94 * t137 + t130) + (t117 * t92 - t134) * t92;
t17 = t119 + t145;
t16 = (rSges(5,2) - t95) * t94 + (-t147 * t82 - t161 - t80) * t92;
t6 = t18 + t133;
t5 = t7 + t133;
t1 = [t93 * (Icges(3,2) * t93 + t128) + t91 * (Icges(3,1) * t91 + t127) + Icges(2,3) + t158 * t83 + t157 * t82 + m(2) * (t71 ^ 2 + t72 ^ 2) + m(3) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2); (t93 * (Icges(3,6) * t92 + t103 * t94) + t91 * (Icges(3,5) * t92 + t106 * t94)) * t142 + (t93 * (-Icges(3,6) * t94 + t103 * t92) + t91 * (-Icges(3,5) * t94 + t106 * t92)) * t141 + m(4) * (t35 * t22 + t34 * t23) + m(5) * (t21 * t16 + t20 * t17) + m(3) * (-t26 * t94 - t27 * t92) * t70 + (t88 / 0.2e1 + t89 / 0.2e1) * (Icges(3,5) * t91 + Icges(3,6) * t93) + t96; m(5) * (t20 ^ 2 + t21 ^ 2 + t5 ^ 2) + m(4) * (t34 ^ 2 + t35 ^ 2 + t6 ^ 2) + t88 * t143 + m(3) * (t129 * t70 ^ 2 + t19 ^ 2) + t140 + (t94 * t143 - t129 * (-Icges(3,3) * t94 + t100 * t92) + t144) * t94; m(5) * (t25 * t16 + t24 * t17) + m(4) * (-t22 * t94 - t23 * t92) * t65 + t96; m(5) * (t24 * t20 + t25 * t21 + t7 * t5) + m(4) * (t18 * t6 + (-t34 * t92 - t35 * t94) * t65) + t107; m(4) * (t129 * t65 ^ 2 + t18 ^ 2) + m(5) * (t24 ^ 2 + t25 ^ 2 + t7 ^ 2) + t107; m(5) * (t16 * t94 + t17 * t92) * t82; m(5) * (-t83 * t5 + (t20 * t92 + t21 * t94) * t82); m(5) * (-t83 * t7 + (t24 * t92 + t25 * t94) * t82); m(5) * (t129 * t82 ^ 2 + t83 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
