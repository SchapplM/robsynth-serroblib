% Calculate joint inertia matrix for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP3_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:34
% EndTime: 2019-12-05 16:13:36
% DurationCPUTime: 0.39s
% Computational Cost: add. (263->101), mult. (547->143), div. (0->0), fcn. (480->6), ass. (0->40)
t101 = MDP(15) + MDP(19);
t124 = t101 * qJ(4);
t121 = pkin(6) * MDP(15);
t103 = MDP(12) + MDP(16);
t120 = MDP(14) + MDP(17);
t89 = cos(qJ(3));
t118 = pkin(6) * t89;
t87 = sin(qJ(3));
t84 = t87 ^ 2;
t88 = sin(qJ(2));
t116 = t84 * t88;
t75 = -t89 * pkin(3) - t87 * qJ(4) - pkin(2);
t86 = cos(pkin(8));
t115 = t86 * t75;
t114 = t88 * t89;
t85 = sin(pkin(8));
t64 = t86 * t118 + t85 * t75;
t112 = pkin(3) * MDP(15);
t62 = -t115 + (pkin(6) * t85 + pkin(4)) * t89;
t111 = t62 * MDP(19);
t66 = (pkin(4) * t85 - qJ(5) * t86 + pkin(6)) * t87;
t110 = t66 * MDP(19);
t73 = -t86 * pkin(4) - t85 * qJ(5) - pkin(3);
t109 = t73 * MDP(19);
t108 = t85 * MDP(13);
t107 = t85 * MDP(18);
t106 = t86 * MDP(12);
t105 = t86 * MDP(16);
t104 = t89 * MDP(10);
t102 = MDP(13) - MDP(18);
t97 = t103 * t85;
t94 = t85 * MDP(12) + t86 * MDP(13);
t93 = t85 * MDP(16) - t86 * MDP(18);
t92 = t102 * t85 - t103 * t86 + t109 - t112;
t90 = cos(qJ(2));
t71 = t86 * t114 - t90 * t85;
t69 = t85 * t114 + t90 * t86;
t63 = -t85 * t118 + t115;
t61 = -t89 * qJ(5) + t64;
t1 = [MDP(1) + t101 * (t84 * t88 ^ 2 + t69 ^ 2 + t71 ^ 2); -t88 * MDP(4) + (pkin(6) * t116 - t69 * t63 + t71 * t64) * MDP(15) + (t71 * t61 + t69 * t62) * MDP(19) + (MDP(3) + t104) * t90 + t102 * (t86 * t116 + t71 * t89) + (-t90 * MDP(11) + t88 * t110 + t120 * (t69 * t86 - t71 * t85)) * t87 + t103 * (t85 * t116 + t69 * t89); MDP(2) + (t63 ^ 2 + t64 ^ 2) * MDP(15) + (t61 ^ 2 + t62 ^ 2 + t66 ^ 2) * MDP(19) + 0.2e1 * (-t87 * MDP(11) + t104) * pkin(2) + 0.2e1 * ((-t63 * t86 - t64 * t85) * MDP(14) + (-t61 * t85 + t62 * t86) * MDP(17) + t93 * t66) * t87 + (MDP(5) + (0.2e1 * t94 + t121) * pkin(6)) * t84 + 0.2e1 * (-t63 * MDP(12) + t64 * MDP(13) + t62 * MDP(16) - t61 * MDP(18) + t87 * MDP(6)) * t89; (-t89 * MDP(11) + (-MDP(10) + t92) * t87) * t88 + (t120 + t124) * (t69 * t85 + t71 * t86); (-t63 * t85 + t64 * t86) * MDP(14) + (t61 * t86 + t62 * t85) * MDP(17) + (-pkin(6) * MDP(11) + MDP(8)) * t89 + (-t105 - t107 + t109) * t66 + (MDP(7) + t93 * t73 - t94 * pkin(3) + (-MDP(10) - t106 + t108 - t112) * pkin(6)) * t87 + (t89 * t97 + (-t63 * MDP(15) + t111) * t85 + (t64 * MDP(15) + t61 * MDP(19) + t102 * t89) * t86) * qJ(4); MDP(9) + (-0.2e1 * t105 - 0.2e1 * t107 + t109) * t73 + (0.2e1 * t106 - 0.2e1 * t108 + t112) * pkin(3) + (0.2e1 * t120 + t124) * (t85 ^ 2 + t86 ^ 2) * qJ(4); t101 * t88 * t87; t110 + (t102 * t86 + t121 + t97) * t87; t92; t101; t69 * MDP(19); t86 * t87 * MDP(17) + t89 * MDP(16) + t111; (MDP(19) * qJ(4) + MDP(17)) * t85; 0; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
