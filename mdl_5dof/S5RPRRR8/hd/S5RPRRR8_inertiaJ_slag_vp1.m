% Calculate joint inertia matrix for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR8_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR8_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:43
% EndTime: 2019-12-31 19:05:45
% DurationCPUTime: 0.96s
% Computational Cost: add. (1717->152), mult. (2802->235), div. (0->0), fcn. (3234->8), ass. (0->78)
t143 = sin(qJ(3));
t144 = cos(qJ(3));
t96 = sin(qJ(1));
t98 = cos(qJ(1));
t73 = -t96 * t143 - t98 * t144;
t72 = t73 ^ 2;
t74 = t98 * t143 - t96 * t144;
t71 = t74 ^ 2;
t97 = cos(qJ(4));
t127 = Icges(5,4) * t97;
t95 = sin(qJ(4));
t105 = Icges(5,2) * t95 - t127;
t128 = Icges(5,4) * t95;
t107 = -Icges(5,1) * t97 + t128;
t94 = qJ(4) + qJ(5);
t88 = sin(t94);
t89 = cos(t94);
t117 = (-t72 / 0.2e1 - t71 / 0.2e1) * (0.2e1 * Icges(6,5) * t88 + 0.2e1 * Icges(6,6) * t89);
t152 = -t97 / 0.2e1;
t153 = -t95 / 0.2e1;
t150 = Icges(5,5) * t153 + Icges(5,6) * t152;
t156 = t117 + (t73 * t150 + t97 * (-Icges(5,6) * t73 + t105 * t74) / 0.2e1 + t95 * (-Icges(5,5) * t73 + t107 * t74) / 0.2e1) * t73 + ((Icges(5,6) * t74 + t105 * t73) * t152 + (Icges(5,5) * t74 + t107 * t73) * t153 + t74 * t150) * t74;
t149 = t74 * t73;
t151 = -rSges(6,1) * t89 + rSges(6,2) * t88;
t37 = t74 * rSges(6,3) + t151 * t73;
t87 = t97 * pkin(4) + pkin(3);
t99 = -pkin(8) - pkin(7);
t25 = t73 * t87 + t74 * t99 - t37;
t140 = rSges(5,2) * t95;
t119 = -pkin(3) + t140;
t142 = rSges(5,1) * t97;
t131 = t74 * rSges(5,3) - t73 * t142;
t65 = t74 * pkin(7);
t27 = -t119 * t73 - t131 - t65;
t133 = t73 * rSges(5,3) + t74 * t142;
t64 = t73 * pkin(7);
t26 = t119 * t74 - t133 - t64;
t102 = -Icges(6,5) * t89 + Icges(6,6) * t88;
t31 = -Icges(6,3) * t73 + t102 * t74;
t32 = Icges(6,3) * t74 + t102 * t73;
t147 = t74 * (-t31 * t149 + t71 * t32) - t73 * (-t32 * t149 + t72 * t31);
t80 = -t95 * rSges(5,1) - t97 * rSges(5,2);
t146 = m(5) * t80;
t70 = -t88 * rSges(6,1) - t89 * rSges(6,2);
t145 = m(6) * t70;
t132 = t73 * t99 - t74 * t87;
t130 = t71 + t72;
t129 = t98 * pkin(1) + t96 * qJ(2);
t122 = t98 * pkin(2) + t129;
t118 = pkin(4) * t95 - t70;
t91 = t98 * qJ(2);
t115 = t91 + (-pkin(1) - pkin(2)) * t96;
t48 = -t74 * rSges(4,1) + t73 * rSges(4,2);
t49 = t73 * rSges(4,1) + t74 * rSges(4,2);
t109 = -t73 * t96 + t74 * t98;
t103 = -Icges(5,5) * t97 + Icges(5,6) * t95;
t101 = -t73 * rSges(6,3) + t151 * t74;
t24 = t101 + t132;
t100 = -t89 ^ 2 * Icges(6,2) + t97 * (-Icges(5,2) * t97 - t128) + t95 * (-Icges(5,1) * t95 - t127) - Icges(4,3) + (-Icges(6,1) * t88 - 0.2e1 * Icges(6,4) * t89) * t88;
t82 = t98 * rSges(2,1) - t96 * rSges(2,2);
t81 = -t96 * rSges(2,1) - t98 * rSges(2,2);
t55 = t98 * rSges(3,1) + t96 * rSges(3,3) + t129;
t54 = t98 * rSges(3,3) + t91 + (-rSges(3,1) - pkin(1)) * t96;
t47 = -t49 + t122;
t46 = -t48 + t115;
t45 = t118 * t73;
t44 = t118 * t74;
t39 = Icges(5,3) * t74 + t103 * t73;
t38 = -Icges(5,3) * t73 + t103 * t74;
t30 = t74 * t101;
t23 = t122 - t27;
t22 = t115 - t26;
t19 = -t25 + t122;
t18 = -t24 + t115;
t9 = t74 * (t74 * t140 - t133) + t73 * (t73 * t140 + t131);
t4 = t73 * t37 + t30;
t3 = t74 * (t74 * pkin(3) + t132 + t64) + t30 + (t73 * pkin(3) - t25 - t65) * t73;
t1 = [Icges(3,2) + Icges(2,3) + m(6) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t46 ^ 2 + t47 ^ 2) + m(3) * (t54 ^ 2 + t55 ^ 2) + m(2) * (t81 ^ 2 + t82 ^ 2) - t100; m(6) * (t96 * t18 - t98 * t19) + m(5) * (t96 * t22 - t98 * t23) + m(4) * (t96 * t46 - t98 * t47) + m(3) * (t96 * t54 - t98 * t55); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t96 ^ 2 + t98 ^ 2); m(6) * (t24 * t18 + t25 * t19) + m(5) * (t26 * t22 + t27 * t23) + m(4) * (t48 * t46 + t49 * t47) + t100; m(4) * (t48 * t96 - t49 * t98) + m(5) * (t26 * t96 - t27 * t98) + m(6) * (t24 * t96 - t25 * t98); m(6) * (t24 ^ 2 + t25 ^ 2) + m(5) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t48 ^ 2 + t49 ^ 2) - t100; m(6) * (t45 * t18 + t44 * t19) + (-t22 * t73 - t23 * t74) * t146 + t156; m(6) * (-t44 * t98 + t45 * t96) + t109 * t146; m(6) * (t45 * t24 + t44 * t25) + (-t26 * t73 - t27 * t74) * t146 - t156; m(5) * (t130 * t80 ^ 2 + t9 ^ 2) + t74 * (-t38 * t149 + t71 * t39) - t73 * (-t39 * t149 + t72 * t38) + m(6) * (t3 ^ 2 + t44 ^ 2 + t45 ^ 2) + t147; (-t18 * t73 - t19 * t74) * t145 + t117; t109 * t145; (-t24 * t73 - t25 * t74) * t145 - t117; m(6) * (t4 * t3 + (-t44 * t74 - t45 * t73) * t70) + t147; m(6) * (t130 * t70 ^ 2 + t4 ^ 2) + t147;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
