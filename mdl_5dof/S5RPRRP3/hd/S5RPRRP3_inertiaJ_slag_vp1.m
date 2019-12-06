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
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:03:10
% EndTime: 2019-12-05 18:03:14
% DurationCPUTime: 1.13s
% Computational Cost: add. (2263->176), mult. (1730->229), div. (0->0), fcn. (1546->8), ass. (0->97)
t171 = Icges(5,4) + Icges(6,4);
t170 = Icges(5,1) + Icges(6,1);
t169 = Icges(5,5) + Icges(6,5);
t168 = Icges(5,2) + Icges(6,2);
t167 = Icges(5,6) + Icges(6,6);
t95 = qJ(3) + qJ(4);
t91 = cos(t95);
t166 = t171 * t91;
t90 = sin(t95);
t165 = t171 * t90;
t162 = t170 * t90 + t166;
t163 = t168 * t91 + t165;
t164 = (t162 + t166) * t91 + (-t163 - t165 + (-t168 + t170) * t91) * t90;
t161 = Icges(5,3) + Icges(6,3);
t160 = -t167 * t90 + t169 * t91;
t157 = 0.2e1 * t167 * t91 + 0.2e1 * t169 * t90;
t94 = qJ(1) + pkin(8);
t88 = sin(t94);
t100 = -pkin(7) - pkin(6);
t93 = -qJ(5) + t100;
t156 = (rSges(6,3) - t93) * t88;
t89 = cos(t94);
t154 = -t160 * t88 + t161 * t89;
t153 = t160 * t89 + t161 * t88;
t150 = t88 * t89;
t86 = t88 ^ 2;
t87 = t89 ^ 2;
t149 = t88 / 0.2e1;
t148 = t89 / 0.2e1;
t96 = sin(qJ(3));
t147 = pkin(3) * t96;
t97 = sin(qJ(1));
t146 = t97 * pkin(1);
t99 = cos(qJ(1));
t145 = t99 * pkin(1);
t98 = cos(qJ(3));
t85 = t98 * pkin(3) + pkin(2);
t144 = pkin(2) - t85;
t143 = rSges(4,1) * t98;
t142 = rSges(5,1) * t91;
t141 = rSges(6,1) * t91;
t140 = t88 * rSges(5,3);
t139 = t88 * t90;
t138 = t88 * t91;
t137 = t88 * t96;
t119 = -rSges(6,2) * t90 + t141;
t68 = pkin(4) * t91 + t85;
t134 = t68 - t85;
t81 = t88 * t100;
t136 = (t81 + (t134 + t119) * t89 + t156) * t89;
t133 = rSges(6,2) * t139 + t89 * rSges(6,3);
t135 = -(t100 - t93) * t89 + t134 * t88 + rSges(6,1) * t138 - t133;
t132 = rSges(5,2) * t139 + t89 * rSges(5,3);
t66 = t90 * rSges(6,1) + t91 * rSges(6,2);
t32 = pkin(4) * t139 + t88 * t66;
t131 = rSges(4,2) * t137 + t89 * rSges(4,3);
t130 = t87 + t86;
t129 = Icges(4,4) * t96;
t128 = Icges(4,4) * t98;
t123 = t154 * t89 * t87 + (t153 * t86 + (t153 * t89 + t154 * t88) * t89) * t88;
t122 = -pkin(4) * t90 - t66;
t121 = -rSges(4,2) * t96 + t143;
t120 = -rSges(5,2) * t90 + t142;
t110 = Icges(4,1) * t98 - t129;
t107 = -Icges(4,2) * t96 + t128;
t104 = Icges(4,5) * t98 - Icges(4,6) * t96;
t101 = (t157 * t148 + t149 * t164) * t89 + (-t148 * t164 + t157 * t149) * t88;
t80 = pkin(3) * t137;
t78 = -t99 * rSges(2,1) + t97 * rSges(2,2);
t77 = -t97 * rSges(2,1) - t99 * rSges(2,2);
t76 = t96 * rSges(4,1) + t98 * rSges(4,2);
t67 = t90 * rSges(5,1) + t91 * rSges(5,2);
t58 = -t89 * rSges(3,1) + t88 * rSges(3,2) - t145;
t57 = -t88 * rSges(3,1) - t89 * rSges(3,2) - t146;
t51 = Icges(4,3) * t88 + t104 * t89;
t50 = Icges(4,3) * t89 - t104 * t88;
t49 = (-t67 - t147) * t89;
t48 = t88 * t67 + t80;
t47 = -rSges(5,1) * t138 + t132;
t33 = t122 * t89;
t31 = (-pkin(6) - t100) * t89 + t144 * t88;
t30 = t89 * (t120 * t89 + t140);
t28 = t89 * (-t88 * pkin(6) - t144 * t89 - t81);
t27 = (t122 - t147) * t89;
t26 = t80 + t32;
t25 = -t145 + (-rSges(4,3) - pkin(6)) * t88 + (-pkin(2) - t121) * t89;
t24 = -t146 + t89 * pkin(6) + (-pkin(2) - t143) * t88 + t131;
t21 = -t140 - t145 + t81 + (-t120 - t85) * t89;
t20 = -t146 - t89 * t100 + (-t85 - t142) * t88 + t132;
t19 = -t145 - t156 + (-t119 - t68) * t89;
t18 = -t146 - t89 * t93 + (-t68 - t141) * t88 + t133;
t17 = -t88 * (-t88 * t143 + t131) + (t88 * rSges(4,3) + t121 * t89) * t89;
t16 = -t88 * t47 + t30;
t7 = t28 + t30 + (-t31 - t47) * t88;
t6 = t135 * t88 + t136;
t1 = t28 + (-t31 + t135) * t88 + t136;
t2 = [t98 * (Icges(4,2) * t98 + t129) + t96 * (Icges(4,1) * t96 + t128) + Icges(2,3) + Icges(3,3) + t163 * t91 + t162 * t90 + m(2) * (t77 ^ 2 + t78 ^ 2) + m(3) * (t57 ^ 2 + t58 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(6) * (t18 ^ 2 + t19 ^ 2); 0; m(3) + m(4) + m(5) + m(6); (t98 * (Icges(4,6) * t89 - t107 * t88) + t96 * (Icges(4,5) * t89 - t110 * t88)) * t148 + (t98 * (Icges(4,6) * t88 + t107 * t89) + t96 * (Icges(4,5) * t88 + t110 * t89)) * t149 + m(5) * (t49 * t20 + t48 * t21) + m(6) * (t27 * t18 + t26 * t19) + m(4) * (-t24 * t89 + t25 * t88) * t76 + (t87 / 0.2e1 + t86 / 0.2e1) * (Icges(4,5) * t96 + Icges(4,6) * t98) + t101; m(4) * t17 + m(5) * t7 + m(6) * t1; m(6) * (t1 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(5) * (t48 ^ 2 + t49 ^ 2 + t7 ^ 2) + m(4) * (t130 * t76 ^ 2 + t17 ^ 2) + t88 * (t50 * t150 + t86 * t51) + t89 * (t51 * t150 + t87 * t50) + t123; m(6) * (t33 * t18 + t32 * t19) + m(5) * (-t20 * t89 + t21 * t88) * t67 + t101; m(5) * t16 + m(6) * t6; m(6) * (t6 * t1 + t32 * t26 + t33 * t27) + m(5) * (t16 * t7 + (t48 * t88 - t49 * t89) * t67) + t123; m(5) * (t130 * t67 ^ 2 + t16 ^ 2) + m(6) * (t32 ^ 2 + t33 ^ 2 + t6 ^ 2) + t123; m(6) * (t88 * t18 + t89 * t19); 0; m(6) * (t89 * t26 + t88 * t27); m(6) * (t89 * t32 + t88 * t33); m(6) * t130;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
