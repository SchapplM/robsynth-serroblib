% Calculate joint inertia matrix for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:11
% EndTime: 2022-01-20 11:02:11
% DurationCPUTime: 0.70s
% Computational Cost: add. (2367->170), mult. (1556->251), div. (0->0), fcn. (1342->10), ass. (0->93)
t99 = qJ(1) + qJ(2);
t93 = sin(t99);
t94 = cos(t99);
t144 = t93 * t94;
t146 = qJ(3) + rSges(4,3);
t100 = sin(pkin(9));
t101 = cos(pkin(9));
t145 = -rSges(4,1) * t101 + rSges(4,2) * t100 - pkin(2);
t88 = t93 ^ 2;
t89 = t94 ^ 2;
t143 = t93 / 0.2e1;
t142 = -t94 / 0.2e1;
t98 = pkin(9) + qJ(4);
t90 = sin(t98);
t91 = cos(t98);
t60 = t90 * rSges(5,1) + t91 * rSges(5,2);
t141 = m(5) * t60;
t92 = qJ(5) + t98;
t84 = sin(t92);
t85 = cos(t92);
t53 = t84 * rSges(6,1) + t85 * rSges(6,2);
t140 = m(6) * t53;
t86 = t101 * pkin(3) + pkin(2);
t139 = rSges(5,1) * t91;
t138 = rSges(6,1) * t85;
t137 = rSges(5,2) * t90;
t136 = rSges(6,2) * t84;
t103 = sin(qJ(1));
t135 = t103 * pkin(1);
t102 = -pkin(7) - qJ(3);
t97 = pkin(8) - t102;
t134 = t94 * t97;
t107 = t93 * rSges(6,3) + (-t136 + t138) * t94;
t132 = t94 * rSges(6,3) + t93 * t136;
t12 = t93 * (t93 * t138 - t132) + t94 * t107;
t63 = pkin(4) * t91 + t86;
t133 = t94 * t63 + t97 * t93;
t131 = t94 * rSges(5,3) + t93 * t137;
t130 = t88 + t89;
t127 = Icges(5,4) * t90;
t126 = Icges(5,4) * t91;
t125 = Icges(6,4) * t84;
t124 = Icges(6,4) * t85;
t111 = -Icges(6,2) * t84 + t124;
t113 = Icges(6,1) * t85 - t125;
t51 = Icges(6,2) * t85 + t125;
t52 = Icges(6,1) * t84 + t124;
t116 = -t51 * t84 + t52 * t85;
t50 = Icges(6,5) * t84 + Icges(6,6) * t85;
t123 = (t116 * t94 + t85 * (Icges(6,6) * t93 + t111 * t94) + t84 * (Icges(6,5) * t93 + t113 * t94) + t93 * t50) * t143 + (t116 * t93 + t85 * (-Icges(6,6) * t94 + t111 * t93) + t84 * (-Icges(6,5) * t94 + t113 * t93) - t94 * t50) * t142;
t109 = Icges(6,5) * t85 - Icges(6,6) * t84;
t30 = -Icges(6,3) * t94 + t109 * t93;
t31 = Icges(6,3) * t93 + t109 * t94;
t122 = -t94 * (-t31 * t144 + t89 * t30) + t93 * (-t30 * t144 + t88 * t31);
t121 = -pkin(4) * t90 - t53;
t62 = t94 * rSges(3,1) - t93 * rSges(3,2);
t61 = -t93 * rSges(3,1) - t94 * rSges(3,2);
t56 = Icges(5,2) * t91 + t127;
t57 = Icges(5,1) * t90 + t126;
t115 = -t56 * t90 + t57 * t91;
t114 = Icges(5,1) * t91 - t127;
t112 = -Icges(5,2) * t90 + t126;
t110 = Icges(5,5) * t91 - Icges(5,6) * t90;
t108 = t93 * rSges(5,3) + (-t137 + t139) * t94;
t55 = Icges(5,5) * t90 + Icges(5,6) * t91;
t106 = t123 + (t115 * t94 + t91 * (Icges(5,6) * t93 + t112 * t94) + t90 * (Icges(5,5) * t93 + t114 * t94) + t93 * t55) * t143 + (t115 * t93 + t91 * (-Icges(5,6) * t94 + t112 * t93) + t90 * (-Icges(5,5) * t94 + t114 * t93) - t94 * t55) * t142;
t105 = Icges(4,2) * t101 ^ 2 + t85 * t51 + t84 * t52 + t91 * t56 + t90 * t57 + Icges(3,3) + (Icges(4,1) * t100 + 0.2e1 * Icges(4,4) * t101) * t100;
t17 = t107 + t133;
t25 = -t145 * t94 + t146 * t93;
t24 = t145 * t93 + t146 * t94;
t66 = t94 * t86;
t21 = -t93 * t102 + t108 + t66;
t16 = t134 + (-t63 - t138) * t93 + t132;
t20 = -t94 * t102 + (-t86 - t139) * t93 + t131;
t104 = cos(qJ(1));
t96 = t104 * pkin(1);
t71 = t104 * rSges(2,1) - t103 * rSges(2,2);
t70 = -t103 * rSges(2,1) - t104 * rSges(2,2);
t48 = t62 + t96;
t47 = t61 - t135;
t37 = Icges(5,3) * t93 + t110 * t94;
t36 = -Icges(5,3) * t94 + t110 * t93;
t29 = t121 * t94;
t28 = t121 * t93;
t23 = t25 + t96;
t22 = t24 - t135;
t19 = t21 + t96;
t18 = t20 - t135;
t15 = t17 + t96;
t14 = t16 - t135;
t13 = t93 * (t93 * t139 - t131) + t94 * t108;
t3 = t94 * (-t66 + t133) + (-t134 + (t63 - t86) * t93) * t93 + t12;
t1 = [Icges(2,3) + m(6) * (t14 ^ 2 + t15 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(3) * (t47 ^ 2 + t48 ^ 2) + m(2) * (t70 ^ 2 + t71 ^ 2) + t105; m(6) * (t16 * t14 + t17 * t15) + m(5) * (t20 * t18 + t21 * t19) + m(4) * (t24 * t22 + t25 * t23) + m(3) * (t61 * t47 + t62 * t48) + t105; m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + m(3) * (t61 ^ 2 + t62 ^ 2) + t105; m(6) * (t93 * t14 - t94 * t15) + m(5) * (t93 * t18 - t94 * t19) + m(4) * (t93 * t22 - t94 * t23); m(6) * (t93 * t16 - t94 * t17) + m(5) * (t93 * t20 - t94 * t21) + m(4) * (t93 * t24 - t94 * t25); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t130; m(6) * (t29 * t14 + t28 * t15) + (-t18 * t94 - t19 * t93) * t141 + t106; m(6) * (t29 * t16 + t28 * t17) + (-t20 * t94 - t21 * t93) * t141 + t106; m(6) * (-t28 * t94 + t29 * t93); m(5) * (t130 * t60 ^ 2 + t13 ^ 2) + t93 * (-t36 * t144 + t88 * t37) - t94 * (-t37 * t144 + t89 * t36) + m(6) * (t28 ^ 2 + t29 ^ 2 + t3 ^ 2) + t122; (-t14 * t94 - t15 * t93) * t140 + t123; (-t16 * t94 - t17 * t93) * t140 + t123; 0; m(6) * (t12 * t3 + (-t28 * t93 - t29 * t94) * t53) + t122; m(6) * (t130 * t53 ^ 2 + t12 ^ 2) + t122;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
