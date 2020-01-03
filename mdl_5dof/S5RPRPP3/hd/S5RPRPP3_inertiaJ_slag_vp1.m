% Calculate joint inertia matrix for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:20
% EndTime: 2019-12-31 18:12:23
% DurationCPUTime: 0.97s
% Computational Cost: add. (1066->157), mult. (1292->230), div. (0->0), fcn. (1114->6), ass. (0->79)
t149 = Icges(4,1) + Icges(6,3);
t148 = Icges(4,5) + Icges(6,5);
t147 = Icges(5,5) + Icges(6,4);
t146 = Icges(6,2) + Icges(5,3);
t71 = pkin(7) + qJ(3);
t66 = cos(t71);
t145 = (Icges(5,6) - Icges(6,6)) * t66;
t65 = sin(t71);
t144 = (Icges(4,4) - Icges(6,6)) * t65;
t143 = t149 * t66 - t144;
t142 = t146 * t65 - t145;
t141 = (-Icges(5,4) + t148) * t66 + (-Icges(4,6) + t147) * t65;
t140 = Icges(5,1) + Icges(6,1) + Icges(4,3);
t128 = Icges(4,6) / 0.2e1;
t139 = t128 - Icges(5,5) / 0.2e1 - Icges(6,4) / 0.2e1;
t133 = -Icges(5,4) / 0.2e1;
t138 = t133 + Icges(4,5) / 0.2e1 + Icges(6,5) / 0.2e1;
t77 = sin(qJ(1));
t137 = -t77 / 0.2e1;
t136 = t77 / 0.2e1;
t78 = cos(qJ(1));
t135 = -t78 / 0.2e1;
t134 = t78 / 0.2e1;
t127 = rSges(6,1) + pkin(4);
t106 = rSges(6,3) + qJ(5);
t72 = t77 ^ 2;
t73 = t78 ^ 2;
t115 = t72 + t73;
t124 = t140 * t78 - t141 * t77;
t123 = t140 * t77 + t141 * t78;
t122 = m(5) / 0.2e1;
t120 = t65 * t78;
t119 = t66 * t78;
t108 = qJ(4) * t65;
t116 = pkin(3) * t119 + t108 * t78;
t118 = t72 * (pkin(3) * t66 + t108) + t78 * t116;
t47 = pkin(3) * t65 - qJ(4) * t66;
t117 = rSges(5,2) * t65 + rSges(5,3) * t66 - t47;
t113 = Icges(4,4) * t66;
t112 = Icges(5,6) * t65;
t107 = rSges(3,3) + qJ(2);
t105 = t122 + m(6) / 0.2e1;
t75 = cos(pkin(7));
t62 = pkin(2) * t75 + pkin(1);
t76 = -pkin(6) - qJ(2);
t104 = t62 * t78 - t76 * t77;
t103 = t66 * rSges(6,2) - t106 * t65 - t47;
t101 = rSges(6,2) * t120 + t106 * t119 + t127 * t77;
t8 = t103 * t77;
t9 = t103 * t78;
t100 = t77 * t8 + t78 * t9;
t99 = rSges(4,1) * t66 - rSges(4,2) * t65;
t91 = -Icges(4,2) * t65 + t113;
t87 = -Icges(5,2) * t66 + t112;
t83 = rSges(4,1) * t119 - rSges(4,2) * t120 + rSges(4,3) * t77;
t82 = rSges(5,1) * t77 - rSges(5,2) * t119 + rSges(5,3) * t120;
t81 = t104 + t116;
t74 = sin(pkin(7));
t80 = rSges(3,1) * t75 - rSges(3,2) * t74 + pkin(1);
t3 = (-t76 + t127) * t78 + (-t62 + (-rSges(6,2) - qJ(4)) * t65 + (-pkin(3) - t106) * t66) * t77;
t4 = t81 + t101;
t79 = m(6) * (t3 * t78 + t4 * t77);
t64 = t66 ^ 2;
t63 = t65 ^ 2;
t52 = rSges(2,1) * t78 - rSges(2,2) * t77;
t51 = -rSges(2,1) * t77 - rSges(2,2) * t78;
t49 = rSges(4,1) * t65 + rSges(4,2) * t66;
t15 = t107 * t77 + t78 * t80;
t14 = t107 * t78 - t77 * t80;
t13 = t117 * t78;
t12 = t117 * t77;
t11 = t104 + t83;
t10 = (rSges(4,3) - t76) * t78 + (-t62 - t99) * t77;
t7 = t78 * t83 + (-t78 * rSges(4,3) + t99 * t77) * t77;
t6 = t81 + t82;
t5 = (rSges(5,1) - t76) * t78 + (-t62 + (rSges(5,2) - pkin(3)) * t66 + (-rSges(5,3) - qJ(4)) * t65) * t77;
t2 = t78 * t82 + (-t78 * rSges(5,1) + (-rSges(5,2) * t66 + rSges(5,3) * t65) * t77) * t77 + t118;
t1 = t101 * t78 + ((rSges(6,2) * t65 + t106 * t66) * t77 - t127 * t78) * t77 + t118;
t16 = [Icges(3,2) * t75 ^ 2 + Icges(2,3) + (Icges(3,1) * t74 + 0.2e1 * Icges(3,4) * t75) * t74 + m(2) * (t51 ^ 2 + t52 ^ 2) + m(3) * (t14 ^ 2 + t15 ^ 2) + m(4) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t5 ^ 2 + t6 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2) + (t112 + (Icges(4,2) + t146) * t66 + t144) * t66 + (t113 + (Icges(5,2) + t149) * t65 + t145) * t65; m(3) * (t14 * t77 - t15 * t78) + m(4) * (t10 * t77 - t11 * t78) + m(5) * (t5 * t77 - t6 * t78) + m(6) * (t3 * t77 - t4 * t78); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t105) * t115; m(5) * (t12 * t6 + t13 * t5) + m(6) * (t3 * t9 + t4 * t8) + m(4) * (-t10 * t78 - t11 * t77) * t49 + ((t136 * t142 + t91 * t137 + t139 * t78) * t78 + (t91 * t134 + t135 * t142 + t139 * t77) * t77) * t66 + ((t87 * t136 + t137 * t143 + t138 * t78) * t78 + (t134 * t143 + t87 * t135 + t138 * t77) * t77) * t65 + t115 * ((t128 - t147 / 0.2e1) * t66 + (t133 + t148 / 0.2e1) * t65); m(5) * (-t12 * t78 + t13 * t77) + m(6) * (t77 * t9 - t78 * t8); m(5) * (t12 ^ 2 + t13 ^ 2 + t2 ^ 2) + m(6) * (t1 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(4) * (t115 * t49 ^ 2 + t7 ^ 2) + t124 * t78 * t73 + (t123 * t72 + (t123 * t78 + t124 * t77) * t78) * t77; 0.2e1 * ((t5 * t78 + t6 * t77) * t122 + t79 / 0.2e1) * t65; 0; m(5) * (-t66 * t2 + (t12 * t77 + t13 * t78) * t65) + m(6) * (-t66 * t1 + t100 * t65); 0.2e1 * t105 * (t115 * t63 + t64); t66 * t79; 0; m(6) * (t65 * t1 + t100 * t66); m(6) * (-0.1e1 + t115) * t66 * t65; m(6) * (t115 * t64 + t63);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t16(1), t16(2), t16(4), t16(7), t16(11); t16(2), t16(3), t16(5), t16(8), t16(12); t16(4), t16(5), t16(6), t16(9), t16(13); t16(7), t16(8), t16(9), t16(10), t16(14); t16(11), t16(12), t16(13), t16(14), t16(15);];
Mq = res;
