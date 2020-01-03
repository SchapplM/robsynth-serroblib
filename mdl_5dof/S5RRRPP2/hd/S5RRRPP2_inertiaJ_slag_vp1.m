% Calculate joint inertia matrix for
% S5RRRPP2
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:29
% EndTime: 2019-12-31 20:51:32
% DurationCPUTime: 1.03s
% Computational Cost: add. (1628->146), mult. (1726->208), div. (0->0), fcn. (1506->6), ass. (0->74)
t103 = cos(qJ(3));
t167 = t103 ^ 2;
t101 = sin(qJ(3));
t166 = Icges(5,4) + Icges(4,5) - Icges(6,5);
t163 = t166 * t101;
t165 = Icges(4,6) - Icges(5,6) + Icges(6,6);
t164 = t165 * t103;
t158 = -t163 - t164;
t100 = qJ(1) + qJ(2);
t96 = sin(t100);
t94 = t96 ^ 2;
t97 = cos(t100);
t95 = t97 ^ 2;
t161 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t159 = t165 * t101 - t166 * t103;
t148 = rSges(6,1) + pkin(4);
t156 = t159 * t96 + t161 * t97;
t155 = -t159 * t97 + t161 * t96;
t154 = 0.2e1 * t101;
t153 = m(5) / 0.2e1;
t152 = m(6) / 0.2e1;
t74 = t101 * rSges(4,1) + t103 * rSges(4,2);
t149 = m(4) * t74;
t102 = sin(qJ(1));
t147 = t102 * pkin(1);
t130 = qJ(4) * t101;
t138 = t103 * t97;
t143 = pkin(3) * t138 + t97 * t130;
t146 = t94 * (pkin(3) * t103 + t130) + t97 * t143;
t71 = t101 * pkin(3) - t103 * qJ(4);
t145 = -t101 * rSges(5,1) + t103 * rSges(5,3) - t71;
t144 = t96 * t101 * rSges(4,2) + t97 * rSges(4,3);
t142 = t97 * pkin(2) + t96 * pkin(7);
t141 = t94 + t95;
t140 = rSges(4,1) * t103;
t139 = t101 * t97;
t137 = -rSges(6,3) - qJ(5);
t129 = rSges(6,2) * t139 + t148 * t138;
t128 = rSges(5,1) * t138 + t96 * rSges(5,2) + rSges(5,3) * t139;
t56 = t97 * rSges(3,1) - t96 * rSges(3,2);
t127 = t142 + t143;
t126 = t103 * rSges(6,2) - t148 * t101 - t71;
t55 = -t96 * rSges(3,1) - t97 * rSges(3,2);
t116 = rSges(4,1) * t138 - rSges(4,2) * t139 + t96 * rSges(4,3);
t23 = t127 + t128;
t29 = t116 + t142;
t92 = t97 * pkin(7);
t28 = t92 + (-pkin(2) - t140) * t96 + t144;
t106 = Icges(3,3) + (Icges(4,2) + Icges(6,2) + Icges(5,3)) * t167 + ((Icges(4,1) + Icges(5,1) + Icges(6,1)) * t101 + (2 * Icges(4,4) - 2 * Icges(6,4) - 2 * Icges(5,5)) * t103) * t101;
t6 = t137 * t96 + t127 + t129;
t105 = -t158 * t95 + (-t158 / 0.2e1 + t164 / 0.2e1 + t163 / 0.2e1) * t94;
t89 = t97 * rSges(5,2);
t22 = t89 + t92 + (-pkin(2) + (-rSges(5,1) - pkin(3)) * t103 + (-rSges(5,3) - qJ(4)) * t101) * t96;
t5 = t92 + t137 * t97 + (-pkin(2) + (-rSges(6,2) - qJ(4)) * t101 + (-pkin(3) - t148) * t103) * t96;
t104 = cos(qJ(1));
t98 = t104 * pkin(1);
t76 = t104 * rSges(2,1) - t102 * rSges(2,2);
t75 = -t102 * rSges(2,1) - t104 * rSges(2,2);
t54 = t56 + t98;
t53 = t55 - t147;
t31 = t145 * t97;
t30 = t145 * t96;
t27 = t126 * t97;
t26 = t126 * t96;
t25 = t29 + t98;
t24 = t28 - t147;
t15 = t98 + t23;
t14 = t22 - t147;
t13 = t96 * (t96 * t140 - t144) + t97 * t116;
t4 = t98 + t6;
t3 = t5 - t147;
t2 = t97 * t128 + (-t89 + (rSges(5,1) * t103 + rSges(5,3) * t101) * t96) * t96 + t146;
t1 = t129 * t97 + (rSges(6,2) * t101 + t148 * t103) * t94 + t146;
t7 = [Icges(2,3) + m(6) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t14 ^ 2 + t15 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + m(3) * (t53 ^ 2 + t54 ^ 2) + m(2) * (t75 ^ 2 + t76 ^ 2) + t106; m(6) * (t5 * t3 + t6 * t4) + m(5) * (t22 * t14 + t23 * t15) + m(4) * (t28 * t24 + t29 * t25) + m(3) * (t55 * t53 + t56 * t54) + t106; m(5) * (t22 ^ 2 + t23 ^ 2) + m(6) * (t5 ^ 2 + t6 ^ 2) + m(4) * (t28 ^ 2 + t29 ^ 2) + m(3) * (t55 ^ 2 + t56 ^ 2) + t106; t105 + m(6) * (t26 * t4 + t27 * t3) + m(5) * (t31 * t14 + t30 * t15) + (-t24 * t97 - t25 * t96) * t149; t105 + m(5) * (t31 * t22 + t30 * t23) + m(6) * (t26 * t6 + t27 * t5) + (-t28 * t97 - t29 * t96) * t149; m(5) * (t2 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t1 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(4) * (t141 * t74 ^ 2 + t13 ^ 2) + t155 * t96 * t94 + (t156 * t95 + (t155 * t97 + t156 * t96) * t96) * t97; ((t3 * t97 + t4 * t96) * t152 + (t14 * t97 + t15 * t96) * t153) * t154; ((t22 * t97 + t23 * t96) * t153 + (t5 * t97 + t6 * t96) * t152) * t154; m(5) * (-t103 * t2 + (t30 * t96 + t31 * t97) * t101) + m(6) * (-t103 * t1 + (t26 * t96 + t27 * t97) * t101); 0.2e1 * (t153 + t152) * (t141 * t101 ^ 2 + t167); m(6) * (-t96 * t3 + t97 * t4); m(6) * (-t96 * t5 + t97 * t6); m(6) * (t97 * t26 - t96 * t27); 0; m(6) * t141;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
