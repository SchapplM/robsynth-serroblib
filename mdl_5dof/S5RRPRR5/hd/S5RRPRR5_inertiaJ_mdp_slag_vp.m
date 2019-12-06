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
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRPRR5_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:34:05
% EndTime: 2019-12-05 18:34:06
% DurationCPUTime: 0.30s
% Computational Cost: add. (382->80), mult. (680->99), div. (0->0), fcn. (715->8), ass. (0->57)
t120 = sin(qJ(2));
t107 = pkin(1) * t120 + qJ(3);
t116 = sin(pkin(9));
t117 = cos(pkin(9));
t139 = t116 ^ 2 + t117 ^ 2;
t141 = t139 * t107;
t118 = sin(qJ(5));
t121 = cos(qJ(5));
t119 = sin(qJ(4));
t122 = cos(qJ(4));
t94 = (-pkin(7) - t107) * t116;
t113 = t117 * pkin(7);
t95 = t107 * t117 + t113;
t137 = -t119 * t95 + t122 * t94;
t99 = t116 * t122 + t117 * t119;
t149 = pkin(8) * t99;
t73 = t137 - t149;
t132 = -t119 * t94 - t122 * t95;
t98 = t116 * t119 - t122 * t117;
t96 = t98 * pkin(8);
t74 = -t132 - t96;
t156 = (-t118 * t74 + t121 * t73) * MDP(23) + (-t118 * t73 - t121 * t74) * MDP(24);
t101 = (-pkin(7) - qJ(3)) * t116;
t102 = qJ(3) * t117 + t113;
t135 = t122 * t101 - t102 * t119;
t75 = t135 - t149;
t130 = -t101 * t119 - t122 * t102;
t76 = -t130 - t96;
t155 = (-t118 * t76 + t121 * t75) * MDP(23) + (-t118 * t75 - t121 * t76) * MDP(24);
t82 = t118 * t99 + t121 * t98;
t83 = -t118 * t98 + t121 * t99;
t154 = t82 * MDP(23) + t83 * MDP(24);
t153 = t98 * MDP(16) + t99 * MDP(17);
t152 = t117 * MDP(7) - t116 * MDP(8);
t151 = 2 * MDP(9);
t150 = pkin(4) * t98;
t123 = cos(qJ(2));
t148 = pkin(1) * t123;
t146 = t83 * MDP(20) - t82 * MDP(21);
t144 = pkin(2) * MDP(10);
t140 = t139 * qJ(3);
t109 = -pkin(2) - t148;
t138 = MDP(10) * t109;
t108 = -pkin(3) * t117 - pkin(2);
t136 = t99 * MDP(13) - t98 * MDP(14) + t146;
t134 = t139 * MDP(10);
t133 = MDP(4) + (MDP(11) * t99 - 0.2e1 * MDP(12) * t98) * t99 + (MDP(18) * t83 - 0.2e1 * MDP(19) * t82) * t83;
t131 = 0.2e1 * t152;
t100 = t108 - t148;
t129 = (MDP(5) * t123 - MDP(6) * t120) * pkin(1);
t128 = (MDP(23) * t121 - MDP(24) * t118) * pkin(4);
t127 = -t152 + t153 + t154;
t126 = 0.2e1 * t153;
t125 = 0.2e1 * t154;
t86 = t108 + t150;
t85 = t100 + t150;
t1 = [MDP(1) + t141 * t151 + t85 * t125 + t107 ^ 2 * t134 + (-t131 + t138) * t109 + t100 * t126 + 0.2e1 * t129 + t133; (t140 + t141) * MDP(9) + (-t109 * pkin(2) + qJ(3) * t141) * MDP(10) + t129 + t133 + t152 * (pkin(2) - t109) + t154 * (t85 + t86) + t153 * (t100 + t108); t140 * t151 + t86 * t125 + qJ(3) ^ 2 * t134 + t108 * t126 + (t131 + t144) * pkin(2) + t133; t127 + t138; t127 - t144; MDP(10); t137 * MDP(16) + t132 * MDP(17) + t136 + t156; t135 * MDP(16) + t130 * MDP(17) + t136 + t155; 0; MDP(15) + MDP(22) + 0.2e1 * t128; t146 + t156; t146 + t155; 0; MDP(22) + t128; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
