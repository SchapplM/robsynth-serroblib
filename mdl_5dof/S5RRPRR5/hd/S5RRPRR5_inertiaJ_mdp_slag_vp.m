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
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRR5_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:57
% EndTime: 2022-01-20 11:02:58
% DurationCPUTime: 0.27s
% Computational Cost: add. (374->79), mult. (667->102), div. (0->0), fcn. (707->8), ass. (0->53)
t118 = sin(qJ(2));
t106 = t118 * pkin(1) + qJ(3);
t114 = sin(pkin(9));
t115 = cos(pkin(9));
t134 = t114 ^ 2 + t115 ^ 2;
t136 = t134 * t106;
t116 = sin(qJ(5));
t119 = cos(qJ(5));
t117 = sin(qJ(4));
t120 = cos(qJ(4));
t93 = (-pkin(7) - t106) * t114;
t111 = t115 * pkin(7);
t94 = t115 * t106 + t111;
t133 = -t117 * t94 + t120 * t93;
t98 = t120 * t114 + t117 * t115;
t142 = t98 * pkin(8);
t72 = t133 - t142;
t129 = -t117 * t93 - t120 * t94;
t97 = t117 * t114 - t120 * t115;
t95 = t97 * pkin(8);
t73 = -t129 - t95;
t148 = (-t116 * t73 + t119 * t72) * MDP(22) + (-t116 * t72 - t119 * t73) * MDP(23);
t100 = (-pkin(7) - qJ(3)) * t114;
t101 = t115 * qJ(3) + t111;
t131 = t120 * t100 - t117 * t101;
t74 = t131 - t142;
t128 = -t117 * t100 - t120 * t101;
t75 = -t128 - t95;
t147 = (-t116 * t75 + t119 * t74) * MDP(22) + (-t116 * t74 - t119 * t75) * MDP(23);
t81 = t116 * t98 + t119 * t97;
t82 = -t116 * t97 + t119 * t98;
t146 = t81 * MDP(22) + t82 * MDP(23);
t145 = t97 * MDP(15) + t98 * MDP(16);
t144 = 2 * MDP(8);
t143 = t97 * pkin(4);
t121 = cos(qJ(2));
t141 = t121 * pkin(1);
t140 = t82 * MDP(19) - t81 * MDP(20);
t137 = t115 * MDP(7);
t135 = t134 * qJ(3);
t107 = -t115 * pkin(3) - pkin(2);
t132 = t98 * MDP(12) - t97 * MDP(13) + t140;
t130 = MDP(4) + (MDP(10) * t98 - 0.2e1 * MDP(11) * t97) * t98 + (MDP(17) * t82 - 0.2e1 * MDP(18) * t81) * t82;
t99 = t107 - t141;
t127 = -t137 + t145 + t146;
t126 = (t121 * MDP(5) - t118 * MDP(6)) * pkin(1);
t125 = (MDP(22) * t119 - MDP(23) * t116) * pkin(4);
t124 = 0.2e1 * t145;
t123 = 0.2e1 * t146;
t108 = -pkin(2) - t141;
t85 = t107 + t143;
t84 = t99 + t143;
t1 = [MDP(1) - 0.2e1 * t108 * t137 + t136 * t144 + (t134 * t106 ^ 2 + t108 ^ 2) * MDP(9) + t99 * t124 + t84 * t123 + 0.2e1 * t126 + t130; (t135 + t136) * MDP(8) + (-t108 * pkin(2) + qJ(3) * t136) * MDP(9) + (pkin(2) - t108) * t137 + t126 + t130 + t146 * (t84 + t85) + t145 * (t107 + t99); 0.2e1 * pkin(2) * t137 + t135 * t144 + (t134 * qJ(3) ^ 2 + pkin(2) ^ 2) * MDP(9) + t85 * t123 + t107 * t124 + t130; t108 * MDP(9) + t127; -pkin(2) * MDP(9) + t127; MDP(9); t133 * MDP(15) + t129 * MDP(16) + t132 + t148; t131 * MDP(15) + t128 * MDP(16) + t132 + t147; 0; MDP(14) + MDP(21) + 0.2e1 * t125; t140 + t148; t140 + t147; 0; MDP(21) + t125; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
