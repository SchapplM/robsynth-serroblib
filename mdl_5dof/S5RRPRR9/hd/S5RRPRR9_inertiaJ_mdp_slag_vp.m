% Calculate joint inertia matrix for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR9_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:47:36
% EndTime: 2021-01-15 21:47:38
% DurationCPUTime: 0.40s
% Computational Cost: add. (522->116), mult. (1002->179), div. (0->0), fcn. (1081->8), ass. (0->63)
t147 = -qJ(3) - pkin(6);
t112 = sin(pkin(9));
t106 = pkin(2) * t112 + pkin(7);
t115 = sin(qJ(4));
t118 = cos(qJ(4));
t124 = MDP(20) * t115 + MDP(21) * t118;
t114 = sin(qJ(5));
t117 = cos(qJ(5));
t101 = t114 * t115 - t117 * t118;
t102 = t114 * t118 + t115 * t117;
t140 = pkin(8) + t106;
t94 = t140 * t115;
t95 = t140 * t118;
t126 = t102 * MDP(24) - t101 * MDP(25) + (-t114 * t95 - t117 * t94) * MDP(27) - (-t114 * t94 + t117 * t95) * MDP(28);
t146 = t115 * MDP(17) + t118 * MDP(18) - t106 * t124 + t126;
t141 = cos(qJ(2));
t109 = -pkin(2) * t141 - pkin(1);
t145 = 0.2e1 * t109;
t144 = -2 * MDP(23);
t143 = 0.2e1 * MDP(28);
t113 = cos(pkin(9));
t116 = sin(qJ(2));
t98 = t112 * t116 - t113 * t141;
t142 = t98 * pkin(4);
t90 = t101 * MDP(27);
t139 = -t102 * MDP(28) - t90;
t99 = t112 * t141 + t113 * t116;
t138 = t115 * t99;
t104 = t147 * t141;
t127 = t147 * t116;
t88 = -t113 * t104 + t112 * t127;
t136 = t118 * t88;
t85 = t98 * pkin(3) - t99 * pkin(7) + t109;
t72 = t136 + (-pkin(8) * t99 + t85) * t115;
t137 = t117 * t72;
t135 = t118 * t99;
t134 = t115 * t118;
t78 = t102 * t99;
t76 = t78 * MDP(25);
t79 = t101 * t99;
t77 = t79 * MDP(24);
t133 = t99 * MDP(12);
t132 = MDP(22) * t102;
t131 = 0.2e1 * t141;
t130 = MDP(19) + MDP(26);
t129 = t98 * MDP(26) - t76 - t77;
t107 = -pkin(2) * t113 - pkin(3);
t128 = MDP(16) * t134;
t73 = -t115 * t88 + t118 * t85;
t71 = -pkin(8) * t135 + t142 + t73;
t68 = -t114 * t72 + t117 * t71;
t86 = -t104 * t112 - t113 * t127;
t125 = MDP(20) * t118 - MDP(21) * t115;
t123 = -MDP(11) - t125;
t122 = (MDP(27) * t117 - MDP(28) * t114) * pkin(4);
t121 = (MDP(17) * t118 - MDP(18) * t115) * t99;
t111 = t118 ^ 2;
t110 = t115 ^ 2;
t103 = -pkin(4) * t118 + t107;
t75 = pkin(4) * t138 + t86;
t74 = t115 * t85 + t136;
t69 = t114 * t71 + t137;
t1 = [MDP(1) + pkin(1) * MDP(9) * t131 + t133 * t145 + (t109 ^ 2 + t86 ^ 2 + t88 ^ 2) * MDP(14) + (MDP(15) * t111 - 0.2e1 * t128) * t99 ^ 2 + t130 * t98 ^ 2 - (-MDP(22) * t79 + t144 * t78) * t79 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t116 + MDP(5) * t131) * t116 + (MDP(11) * t145 + 0.2e1 * t121 - 0.2e1 * t76 - 0.2e1 * t77) * t98 + 0.2e1 * (t86 * t99 - t88 * t98) * MDP(13) + 0.2e1 * (t138 * t86 + t73 * t98) * MDP(20) + 0.2e1 * (t135 * t86 - t74 * t98) * MDP(21) + 0.2e1 * (t68 * t98 + t75 * t78) * MDP(27) + (-t69 * t98 - t75 * t79) * t143; t116 * MDP(6) + t141 * MDP(7) - t88 * MDP(12) - t79 * t132 + (t101 * t79 - t102 * t78) * MDP(23) + (t101 * t75 + t103 * t78) * MDP(27) + (t102 * t75 - t103 * t79) * MDP(28) + (-MDP(10) * t141 - t116 * MDP(9)) * pkin(6) + t123 * t86 + (MDP(15) * t134 + (-t110 + t111) * MDP(16) + t124 * t107) * t99 + t146 * t98 + ((-t112 * t98 - t113 * t99) * MDP(13) + (t112 * t88 - t113 * t86) * MDP(14)) * pkin(2); 0.2e1 * t128 + 0.2e1 * t103 * t90 + t110 * MDP(15) + MDP(8) + (t112 ^ 2 + t113 ^ 2) * MDP(14) * pkin(2) ^ 2 + (t101 * t144 + t103 * t143 + t132) * t102 - 0.2e1 * t125 * t107 + 0.2e1 * (t113 * MDP(11) - t112 * MDP(12)) * pkin(2); t133 + t109 * MDP(14) + (-t123 + t139) * t98; 0; MDP(14); t98 * MDP(19) + t73 * MDP(20) - t74 * MDP(21) + (t117 * t142 + t68) * MDP(27) + (-t137 + (-t71 - t142) * t114) * MDP(28) + t121 + t129; t146; t125 + t139; 0.2e1 * t122 + t130; t68 * MDP(27) - MDP(28) * t69 + t129; t126; t139; MDP(26) + t122; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
