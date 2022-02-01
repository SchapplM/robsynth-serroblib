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
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRPR3_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:43:16
% EndTime: 2022-01-20 11:43:17
% DurationCPUTime: 0.28s
% Computational Cost: add. (481->101), mult. (849->131), div. (0->0), fcn. (874->8), ass. (0->58)
t125 = sin(qJ(5));
t128 = cos(qJ(5));
t123 = sin(pkin(9));
t124 = cos(pkin(9));
t126 = sin(qJ(3));
t129 = cos(qJ(3));
t108 = t123 * t129 + t124 * t126;
t154 = t108 * pkin(8);
t127 = sin(qJ(2));
t116 = t127 * pkin(1) + pkin(7);
t105 = (-qJ(4) - t116) * t126;
t122 = t129 * qJ(4);
t106 = t129 * t116 + t122;
t85 = t124 * t105 - t123 * t106;
t76 = t85 - t154;
t107 = t123 * t126 - t124 * t129;
t104 = t107 * pkin(8);
t86 = t123 * t105 + t124 * t106;
t77 = -t104 + t86;
t160 = (-t125 * t77 + t128 * t76) * MDP(23) + (-t125 * t76 - t128 * t77) * MDP(24);
t111 = (-qJ(4) - pkin(7)) * t126;
t112 = t129 * pkin(7) + t122;
t94 = t124 * t111 - t123 * t112;
t80 = t94 - t154;
t95 = t123 * t111 + t124 * t112;
t81 = -t104 + t95;
t159 = (-t125 * t81 + t128 * t80) * MDP(23) + (-t125 * t80 - t128 * t81) * MDP(24);
t92 = t128 * t107 + t125 * t108;
t93 = -t125 * t107 + t128 * t108;
t158 = t92 * MDP(23) + t93 * MDP(24);
t147 = MDP(17) * pkin(3);
t157 = t107 * MDP(14) + t108 * MDP(15);
t136 = t129 * MDP(12) - t126 * MDP(13);
t156 = 2 * MDP(16);
t155 = pkin(3) * t123;
t130 = cos(qJ(2));
t153 = t130 * pkin(1);
t151 = -t86 * t107 - t85 * t108;
t150 = -t95 * t107 - t94 * t108;
t149 = t93 * MDP(20) - t92 * MDP(21);
t115 = t124 * pkin(3) + pkin(4);
t146 = (t128 * t115 - t125 * t155) * MDP(23);
t144 = (-t125 * t115 - t128 * t155) * MDP(24);
t118 = -t129 * pkin(3) - pkin(2);
t110 = t118 - t153;
t143 = t110 * MDP(17);
t142 = t118 * MDP(17);
t139 = t157 + t158;
t138 = t129 * MDP(10) + t126 * MDP(9) + (-t107 * t123 - t108 * t124) * pkin(3) * MDP(16) + t149;
t137 = MDP(4) + (MDP(18) * t93 - 0.2e1 * MDP(19) * t92) * t93 + (MDP(7) * t126 + 0.2e1 * MDP(8) * t129) * t126;
t97 = t107 * pkin(4) + t118;
t135 = -t126 * MDP(12) - t129 * MDP(13);
t134 = 0.2e1 * t157;
t133 = (t130 * MDP(5) - t127 * MDP(6)) * pkin(1);
t132 = 0.2e1 * t158;
t117 = -pkin(2) - t153;
t96 = t97 - t153;
t1 = [MDP(1) + t151 * t156 + (t85 ^ 2 + t86 ^ 2) * MDP(17) + t96 * t132 + (t134 + t143) * t110 + t137 - 0.2e1 * t136 * t117 + 0.2e1 * t133; (t150 + t151) * MDP(16) + (t110 * t118 + t85 * t94 + t86 * t95) * MDP(17) + t133 + t137 + t136 * (pkin(2) - t117) + t158 * (t96 + t97) + t157 * (t110 + t118); t150 * t156 + (t94 ^ 2 + t95 ^ 2) * MDP(17) + t97 * t132 + (t134 + t142) * t118 + 0.2e1 * t136 * pkin(2) + t137; t85 * MDP(14) - t86 * MDP(15) + t135 * t116 + (t123 * t86 + t124 * t85) * t147 + t138 + t160; t94 * MDP(14) - t95 * MDP(15) + t135 * pkin(7) + (t123 * t95 + t124 * t94) * t147 + t138 + t159; MDP(11) + MDP(22) + 0.2e1 * t146 + 0.2e1 * t144 + (0.2e1 * t124 * MDP(14) - 0.2e1 * t123 * MDP(15) + (t123 ^ 2 + t124 ^ 2) * t147) * pkin(3); t139 + t143; t139 + t142; 0; MDP(17); t149 + t160; t149 + t159; MDP(22) + t144 + t146; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
