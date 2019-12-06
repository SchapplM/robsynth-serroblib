% Calculate joint inertia matrix for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR7_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:36
% EndTime: 2019-12-05 16:37:38
% DurationCPUTime: 0.50s
% Computational Cost: add. (313->119), mult. (726->190), div. (0->0), fcn. (735->10), ass. (0->56)
t93 = cos(qJ(3));
t123 = pkin(7) * t93;
t90 = sin(qJ(3));
t80 = -pkin(3) * t93 - qJ(4) * t90 - pkin(2);
t85 = sin(pkin(10));
t87 = cos(pkin(10));
t72 = t87 * t123 + t85 * t80;
t68 = -pkin(8) * t93 + t72;
t73 = (pkin(4) * t85 - pkin(8) * t87 + pkin(7)) * t90;
t118 = t87 * t90;
t89 = sin(qJ(5));
t92 = cos(qJ(5));
t76 = t89 * t118 + t92 * t93;
t77 = t92 * t118 - t89 * t93;
t96 = -(-t68 * t89 + t73 * t92) * MDP(21) + (t68 * t92 + t73 * t89) * MDP(22) - t77 * MDP(18) + t76 * MDP(19);
t109 = qJ(4) * MDP(15) + MDP(14);
t125 = -2 * MDP(17);
t124 = pkin(7) * t85;
t86 = sin(pkin(5));
t94 = cos(qJ(2));
t119 = t86 * t94;
t91 = sin(qJ(2));
t120 = t86 * t91;
t88 = cos(pkin(5));
t75 = t93 * t120 + t88 * t90;
t64 = t87 * t119 + t75 * t85;
t122 = t64 * t85;
t121 = t80 * t87;
t117 = pkin(3) * MDP(15);
t116 = qJ(4) * t87;
t115 = MDP(12) * t87;
t114 = t87 * MDP(14);
t113 = t92 * MDP(16);
t112 = t93 * MDP(10);
t111 = t93 * MDP(12);
t107 = MDP(13) * t85 - t117;
t106 = MDP(13) * t87 + MDP(15) * pkin(7);
t104 = MDP(18) * t92 - MDP(19) * t89;
t66 = -t85 * t119 + t75 * t87;
t74 = t90 * t120 - t88 * t93;
t60 = -t66 * t89 + t74 * t92;
t61 = t66 * t92 + t74 * t89;
t103 = MDP(21) * t60 - MDP(22) * t61;
t79 = -pkin(4) * t87 - pkin(8) * t85 - pkin(3);
t101 = -(-t89 * t116 + t79 * t92) * MDP(21) + (t92 * t116 + t79 * t89) * MDP(22);
t100 = t76 * MDP(21) + t77 * MDP(22);
t99 = MDP(21) * t92 - MDP(22) * t89;
t98 = MDP(12) + t99;
t97 = -t87 * MDP(20) - t101;
t95 = qJ(4) ^ 2;
t84 = t90 ^ 2;
t83 = t87 ^ 2;
t82 = t85 ^ 2;
t71 = -t85 * t123 + t121;
t67 = -t121 + (pkin(4) + t124) * t93;
t1 = [MDP(1) + (t64 ^ 2 + t66 ^ 2 + t74 ^ 2) * MDP(15); (t93 * MDP(13) + t72 * MDP(15)) * t66 + (-t71 * MDP(15) + t100 + t111) * t64 + (t64 * t114 + (-MDP(14) * t66 + t103) * t85 + (MDP(12) * t85 + t106) * t74) * t90 + (-t91 * MDP(4) + (-t90 * MDP(11) + MDP(3) + t112) * t94) * t86; MDP(2) + 0.2e1 * pkin(2) * t112 + (t71 ^ 2 + t72 ^ 2) * MDP(15) + (pkin(7) ^ 2 * MDP(15) + t82 * MDP(20) + MDP(5)) * t84 + (t77 * MDP(16) + t76 * t125) * t77 + 0.2e1 * (t84 * t124 - t71 * t93) * MDP(12) + 0.2e1 * (pkin(7) * t84 * t87 + t72 * t93) * MDP(13) + 0.2e1 * t100 * t67 + 0.2e1 * (-pkin(2) * MDP(11) + t93 * MDP(6) - t71 * t114 + (-t72 * MDP(14) - t96) * t85) * t90; -t75 * MDP(11) + (t89 * t122 - t60 * t87) * MDP(21) + (t92 * t122 + t61 * t87) * MDP(22) + (-MDP(10) + t107 - t115) * t74 + t109 * (t66 * t87 + t122); (-pkin(7) * MDP(11) + MDP(8)) * t93 + (MDP(7) + t104 * t82 + (-MDP(10) - t117) * pkin(7)) * t90 + (-t90 * pkin(7) * MDP(12) + (-pkin(3) * t90 + qJ(4) * t93) * MDP(13) + t109 * t72 + t96) * t87 + (qJ(4) * t111 + t77 * t113 + (-t76 * t92 - t77 * t89) * MDP(17) + (qJ(4) * t76 + t67 * t89) * MDP(21) + (qJ(4) * t77 + t67 * t92) * MDP(22) + (-pkin(3) * MDP(12) + pkin(7) * MDP(13) + t97) * t90 - t109 * t71) * t85; MDP(9) + 0.2e1 * pkin(3) * t115 + (pkin(3) ^ 2 + t83 * t95) * MDP(15) + t83 * MDP(20) + (t95 * MDP(15) + (t89 * t125 + t113) * t92) * t82 + 0.2e1 * (-pkin(3) * MDP(13) - t104 * t87) * t85 + 0.2e1 * t101 * t87 + 0.2e1 * (t83 * MDP(14) + (t89 * MDP(21) + t92 * MDP(22) + MDP(14)) * t82) * qJ(4); t74 * MDP(15); (t98 * t85 + t106) * t90; -t98 * t87 + t107; MDP(15); t103; t85 * t90 * MDP(20) - t96; t104 * t85 + t97; t99; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
