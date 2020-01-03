% Calculate joint inertia matrix for
% S4RRRP4
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP4_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP4_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:10
% EndTime: 2019-12-31 17:15:13
% DurationCPUTime: 1.23s
% Computational Cost: add. (1239->155), mult. (1556->224), div. (0->0), fcn. (1406->6), ass. (0->90)
t168 = Icges(4,4) + Icges(5,4);
t167 = Icges(4,1) + Icges(5,1);
t166 = Icges(4,2) + Icges(5,2);
t89 = qJ(2) + qJ(3);
t80 = cos(t89);
t165 = t168 * t80;
t79 = sin(t89);
t164 = t168 * t79;
t163 = Icges(4,5) + Icges(5,5);
t162 = Icges(4,6) + Icges(5,6);
t161 = -t166 * t79 + t165;
t160 = t167 * t80 - t164;
t91 = sin(qJ(1));
t93 = cos(qJ(1));
t159 = t161 * t91 - t162 * t93;
t158 = t161 * t93 + t162 * t91;
t157 = t160 * t91 - t163 * t93;
t156 = t160 * t93 + t163 * t91;
t155 = -t162 * t79 + t163 * t80;
t154 = Icges(4,3) + Icges(5,3);
t153 = t166 * t80 + t164;
t152 = t167 * t79 + t165;
t151 = t154 * t91 + t155 * t93;
t150 = t154 * t93 - t155 * t91;
t149 = -t156 * t80 + t158 * t79;
t148 = -t157 * t80 + t159 * t79;
t147 = t91 * pkin(5);
t133 = t80 * t93;
t134 = t79 * t93;
t92 = cos(qJ(2));
t78 = t92 * pkin(2) + pkin(1);
t66 = pkin(3) * t80 + t78;
t146 = rSges(5,1) * t133 - rSges(5,2) * t134 + t91 * rSges(5,3) + t93 * t66;
t145 = t162 * t80 + t163 * t79;
t87 = t91 ^ 2;
t88 = t93 ^ 2;
t128 = t87 + t88;
t144 = -rSges(5,1) * t80 + rSges(5,2) * t79 - t66;
t143 = t152 * t80 - t153 * t79;
t142 = t150 * t88 + (t149 * t91 + (-t148 + t151) * t93) * t91;
t90 = sin(qJ(2));
t101 = Icges(3,5) * t92 - Icges(3,6) * t90;
t141 = t91 * (Icges(3,3) * t91 + t101 * t93);
t140 = t91 / 0.2e1;
t139 = -t93 / 0.2e1;
t94 = -pkin(6) - pkin(5);
t138 = (t151 * t87 + ((-t149 + t150) * t91 + t148 * t93) * t93) * t91;
t137 = pkin(2) * t90;
t136 = rSges(3,1) * t92;
t135 = rSges(3,2) * t90;
t132 = t93 * rSges(3,3);
t74 = t93 * t78;
t85 = t93 * pkin(5);
t131 = t91 * (t85 + (-pkin(1) + t78) * t91) + t93 * (-t93 * pkin(1) - t147 + t74);
t118 = rSges(4,1) * t80 - rSges(4,2) * t79;
t98 = rSges(4,1) * t133 - rSges(4,2) * t134 + t91 * rSges(4,3);
t16 = t91 * (-t93 * rSges(4,3) + t118 * t91) + t93 * t98;
t130 = t91 * rSges(3,3) + t93 * t136;
t127 = Icges(3,4) * t90;
t126 = Icges(3,4) * t92;
t65 = t79 * rSges(4,1) + t80 * rSges(4,2);
t121 = -t65 - t137;
t120 = -t80 * rSges(5,2) + (-rSges(5,1) - pkin(3)) * t79;
t6 = (-t74 + t146) * t93 + (-t93 * rSges(5,3) + (-t78 - t144) * t91) * t91;
t119 = -t135 + t136;
t108 = t142 * t93 + t138;
t107 = Icges(3,1) * t92 - t127;
t104 = -Icges(3,2) * t90 + t126;
t96 = t120 - t137;
t95 = (t143 * t93 + t145 * t91 + t156 * t79 + t158 * t80) * t140 + (t143 * t91 - t145 * t93 + t157 * t79 + t159 * t80) * t139;
t86 = -qJ(4) + t94;
t73 = t93 * rSges(2,1) - t91 * rSges(2,2);
t72 = -t91 * rSges(2,1) - t93 * rSges(2,2);
t71 = t90 * rSges(3,1) + t92 * rSges(3,2);
t37 = t121 * t93;
t36 = t121 * t91;
t35 = t120 * t93;
t34 = t120 * t91;
t27 = t147 + (pkin(1) - t135) * t93 + t130;
t26 = t132 + t85 + (-pkin(1) - t119) * t91;
t25 = t96 * t93;
t24 = t96 * t91;
t23 = -t91 * t94 + t74 + t98;
t22 = (rSges(4,3) - t94) * t93 + (-t118 - t78) * t91;
t19 = -t91 * t86 + t146;
t18 = (rSges(5,3) - t86) * t93 + t144 * t91;
t17 = t93 * (-t93 * t135 + t130) + (t119 * t91 - t132) * t91;
t7 = t16 + t131;
t3 = t6 + t131;
t1 = [t92 * (Icges(3,2) * t92 + t127) + t90 * (Icges(3,1) * t90 + t126) + Icges(2,3) + t153 * t80 + t152 * t79 + m(2) * (t72 ^ 2 + t73 ^ 2) + m(3) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2); (t92 * (Icges(3,6) * t91 + t104 * t93) + t90 * (Icges(3,5) * t91 + t107 * t93)) * t140 + (t92 * (-Icges(3,6) * t93 + t104 * t91) + t90 * (-Icges(3,5) * t93 + t107 * t91)) * t139 + m(4) * (t37 * t22 + t36 * t23) + m(5) * (t25 * t18 + t24 * t19) + m(3) * (-t26 * t93 - t27 * t91) * t71 + (t87 / 0.2e1 + t88 / 0.2e1) * (Icges(3,5) * t90 + Icges(3,6) * t92) + t95; m(5) * (t24 ^ 2 + t25 ^ 2 + t3 ^ 2) + m(4) * (t36 ^ 2 + t37 ^ 2 + t7 ^ 2) + t87 * t141 + m(3) * (t128 * t71 ^ 2 + t17 ^ 2) + t138 + (t93 * t141 - t128 * (-Icges(3,3) * t93 + t101 * t91) + t142) * t93; m(5) * (t35 * t18 + t34 * t19) + m(4) * (-t22 * t93 - t23 * t91) * t65 + t95; m(5) * (t34 * t24 + t35 * t25 + t6 * t3) + m(4) * (t16 * t7 + (-t36 * t91 - t37 * t93) * t65) + t108; m(4) * (t128 * t65 ^ 2 + t16 ^ 2) + m(5) * (t34 ^ 2 + t35 ^ 2 + t6 ^ 2) + t108; m(5) * (t91 * t18 - t93 * t19); m(5) * (-t93 * t24 + t91 * t25); m(5) * (-t93 * t34 + t91 * t35); m(5) * t128;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
