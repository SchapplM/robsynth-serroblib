% Calculate joint inertia matrix for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:08
% EndTime: 2019-12-31 20:53:11
% DurationCPUTime: 1.04s
% Computational Cost: add. (1667->146), mult. (1782->210), div. (0->0), fcn. (1554->6), ass. (0->79)
t107 = cos(qJ(3));
t103 = t107 ^ 2;
t105 = sin(qJ(3));
t175 = Icges(5,4) - Icges(4,5) - Icges(6,5);
t172 = t175 * t105;
t174 = Icges(6,4) + Icges(5,5) - Icges(4,6);
t173 = t174 * t107;
t166 = -t172 - t173;
t104 = qJ(1) + qJ(2);
t99 = sin(t104);
t97 = t99 ^ 2;
t100 = cos(t104);
t98 = t100 ^ 2;
t171 = Icges(5,1) + Icges(6,1) + Icges(4,3);
t168 = t174 * t105 - t175 * t107;
t167 = rSges(6,1) + pkin(4);
t164 = rSges(6,3) + qJ(5);
t163 = t168 * t100 + t171 * t99;
t162 = t171 * t100 - t168 * t99;
t161 = 0.2e1 * t105;
t160 = m(5) / 0.2e1;
t72 = t105 * rSges(4,1) + t107 * rSges(4,2);
t158 = m(4) * t72;
t156 = pkin(3) * t107;
t106 = sin(qJ(1));
t155 = t106 * pkin(1);
t137 = t100 * t107;
t138 = t100 * t105;
t150 = pkin(3) * t137 + qJ(4) * t138;
t154 = t97 * (qJ(4) * t105 + t156) + t100 * t150;
t70 = t105 * pkin(3) - t107 * qJ(4);
t153 = t105 * rSges(5,2) + t107 * rSges(5,3) - t70;
t146 = t105 * t99;
t152 = rSges(4,2) * t146 + t100 * rSges(4,3);
t145 = t107 * t99;
t151 = t100 * rSges(5,1) + rSges(5,2) * t145;
t149 = t100 * pkin(2) + t99 * pkin(7);
t148 = t167 * t100;
t147 = t97 + t98;
t56 = t100 * rSges(3,1) - t99 * rSges(3,2);
t136 = t149 + t150;
t135 = rSges(6,2) * t138 + t164 * t137 + t167 * t99;
t134 = t107 * rSges(6,2) - t164 * t105 - t70;
t55 = -t99 * rSges(3,1) - t100 * rSges(3,2);
t26 = t134 * t99;
t27 = t134 * t100;
t133 = t100 * t27 + t26 * t99;
t114 = rSges(4,1) * t137 - rSges(4,2) * t138 + t99 * rSges(4,3);
t113 = t99 * rSges(5,1) - rSges(5,2) * t137 + rSges(5,3) * t138;
t94 = t100 * pkin(7);
t5 = t94 + (-pkin(2) + (-rSges(6,2) - qJ(4)) * t105 + (-pkin(3) - t164) * t107) * t99 + t148;
t3 = t5 - t155;
t108 = cos(qJ(1));
t101 = t108 * pkin(1);
t6 = t135 + t136;
t4 = t101 + t6;
t112 = m(6) * (t100 * t3 + t4 * t99);
t111 = m(6) * (t100 * t5 + t6 * t99);
t29 = t114 + t149;
t28 = t94 + (-rSges(4,1) * t107 - pkin(2)) * t99 + t152;
t23 = t113 + t136;
t110 = Icges(3,3) + (Icges(4,2) + Icges(6,2) + Icges(5,3)) * t103 + ((Icges(4,1) + Icges(5,2) + Icges(6,3)) * t105 + (2 * Icges(4,4) + 2 * Icges(5,6) - 2 * Icges(6,6)) * t107) * t105;
t109 = t166 * t97 + (-t173 / 0.2e1 - t172 / 0.2e1 + t166 / 0.2e1) * t98;
t22 = t94 + (-t156 - pkin(2) + (-rSges(5,3) - qJ(4)) * t105) * t99 + t151;
t102 = t105 ^ 2;
t75 = t108 * rSges(2,1) - t106 * rSges(2,2);
t73 = -t106 * rSges(2,1) - t108 * rSges(2,2);
t54 = t101 + t56;
t53 = t55 - t155;
t31 = t153 * t100;
t30 = t153 * t99;
t25 = t101 + t29;
t24 = t28 - t155;
t15 = t101 + t23;
t14 = t22 - t155;
t13 = t99 * (rSges(4,1) * t145 - t152) + t100 * t114;
t2 = t99 * (rSges(5,3) * t146 - t151) + t100 * t113 + t154;
t1 = t135 * t100 + ((rSges(6,2) * t105 + t164 * t107) * t99 - t148) * t99 + t154;
t7 = [Icges(2,3) + m(6) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t14 ^ 2 + t15 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + m(3) * (t53 ^ 2 + t54 ^ 2) + m(2) * (t73 ^ 2 + t75 ^ 2) + t110; m(6) * (t5 * t3 + t6 * t4) + m(5) * (t22 * t14 + t23 * t15) + m(4) * (t28 * t24 + t29 * t25) + m(3) * (t55 * t53 + t56 * t54) + t110; m(5) * (t22 ^ 2 + t23 ^ 2) + m(6) * (t5 ^ 2 + t6 ^ 2) + m(4) * (t28 ^ 2 + t29 ^ 2) + m(3) * (t55 ^ 2 + t56 ^ 2) + t110; t109 + (-t100 * t24 - t25 * t99) * t158 + m(6) * (t26 * t4 + t27 * t3) + m(5) * (t31 * t14 + t30 * t15); t109 + m(5) * (t31 * t22 + t30 * t23) + m(6) * (t26 * t6 + t27 * t5) + (-t100 * t28 - t29 * t99) * t158; m(6) * (t1 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(5) * (t2 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(4) * (t147 * t72 ^ 2 + t13 ^ 2) + t163 * t99 * t97 + (t162 * t98 + (t163 * t100 + t162 * t99) * t99) * t100; (t112 / 0.2e1 + (t100 * t14 + t15 * t99) * t160) * t161; ((t100 * t22 + t23 * t99) * t160 + t111 / 0.2e1) * t161; m(6) * (-t107 * t1 + t133 * t105) + m(5) * (-t107 * t2 + (t100 * t31 + t30 * t99) * t105); 0.2e1 * (t160 + m(6) / 0.2e1) * (t147 * t102 + t103); t107 * t112; t107 * t111; m(6) * (t105 * t1 + t133 * t107); m(6) * (-0.1e1 + t147) * t107 * t105; m(6) * (t147 * t103 + t102);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
