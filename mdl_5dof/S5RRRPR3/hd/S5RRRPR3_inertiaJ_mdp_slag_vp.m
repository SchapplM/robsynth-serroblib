% Calculate joint inertia matrix for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRRPR3_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:59
% EndTime: 2019-12-05 18:43:00
% DurationCPUTime: 0.23s
% Computational Cost: add. (423->89), mult. (753->120), div. (0->0), fcn. (784->8), ass. (0->53)
t125 = sin(qJ(5));
t128 = cos(qJ(5));
t123 = sin(pkin(9));
t124 = cos(pkin(9));
t126 = sin(qJ(3));
t129 = cos(qJ(3));
t107 = t123 * t129 + t124 * t126;
t150 = t107 * pkin(8);
t127 = sin(qJ(2));
t116 = t127 * pkin(1) + pkin(7);
t104 = (-qJ(4) - t116) * t126;
t122 = t129 * qJ(4);
t105 = t129 * t116 + t122;
t84 = t124 * t104 - t123 * t105;
t75 = t84 - t150;
t106 = -t123 * t126 + t124 * t129;
t103 = t106 * pkin(8);
t85 = t123 * t104 + t124 * t105;
t76 = t103 + t85;
t154 = (-t125 * t76 + t128 * t75) * MDP(21) + (-t125 * t75 - t128 * t76) * MDP(22);
t112 = (-qJ(4) - pkin(7)) * t126;
t113 = t129 * pkin(7) + t122;
t93 = t124 * t112 - t123 * t113;
t79 = t93 - t150;
t94 = t123 * t112 + t124 * t113;
t80 = t103 + t94;
t153 = (-t125 * t80 + t128 * t79) * MDP(21) + (-t125 * t79 - t128 * t80) * MDP(22);
t91 = -t128 * t106 + t125 * t107;
t92 = t125 * t106 + t128 * t107;
t145 = t91 * MDP(21) + t92 * MDP(22);
t135 = t129 * MDP(12) - t126 * MDP(13);
t152 = 2 * MDP(14);
t151 = pkin(3) * t123;
t130 = cos(qJ(2));
t149 = t130 * pkin(1);
t147 = t85 * t106 - t84 * t107;
t146 = t94 * t106 - t93 * t107;
t144 = t92 * MDP(18) - t91 * MDP(19);
t142 = MDP(15) * pkin(3);
t115 = t124 * pkin(3) + pkin(4);
t141 = (t128 * t115 - t125 * t151) * MDP(21);
t140 = (-t125 * t115 - t128 * t151) * MDP(22);
t118 = -t129 * pkin(3) - pkin(2);
t137 = t129 * MDP(10) + t126 * MDP(9) + (t106 * t123 - t107 * t124) * pkin(3) * MDP(14) + t144;
t136 = MDP(4) + (MDP(16) * t92 - 0.2e1 * MDP(17) * t91) * t92 + (MDP(7) * t126 + 0.2e1 * MDP(8) * t129) * t126;
t96 = -t106 * pkin(4) + t118;
t134 = -t126 * MDP(12) - t129 * MDP(13);
t133 = (t130 * MDP(5) - t127 * MDP(6)) * pkin(1);
t132 = 0.2e1 * t145;
t117 = -pkin(2) - t149;
t111 = t118 - t149;
t95 = t96 - t149;
t1 = [MDP(1) + t147 * t152 + (t111 ^ 2 + t84 ^ 2 + t85 ^ 2) * MDP(15) + t95 * t132 + t136 - 0.2e1 * t135 * t117 + 0.2e1 * t133; (t146 + t147) * MDP(14) + (t111 * t118 + t84 * t93 + t85 * t94) * MDP(15) + t133 + t136 + t135 * (pkin(2) - t117) + t145 * (t95 + t96); t146 * t152 + (t118 ^ 2 + t93 ^ 2 + t94 ^ 2) * MDP(15) + t96 * t132 + 0.2e1 * t135 * pkin(2) + t136; t134 * t116 + (t123 * t85 + t124 * t84) * t142 + t137 + t154; t134 * pkin(7) + (t123 * t94 + t124 * t93) * t142 + t137 + t153; MDP(11) + MDP(20) + (t123 ^ 2 + t124 ^ 2) * MDP(15) * pkin(3) ^ 2 + 0.2e1 * t141 + 0.2e1 * t140; t111 * MDP(15) + t145; t118 * MDP(15) + t145; 0; MDP(15); t144 + t154; t144 + t153; MDP(20) + t140 + t141; 0; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
