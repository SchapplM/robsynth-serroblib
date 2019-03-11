% Calculate joint inertia matrix for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPPRRR5_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:29:02
% EndTime: 2019-03-09 02:29:03
% DurationCPUTime: 0.30s
% Computational Cost: add. (292->92), mult. (489->119), div. (0->0), fcn. (477->6), ass. (0->45)
t91 = sin(qJ(6));
t94 = cos(qJ(6));
t101 = MDP(29) * t94 - MDP(30) * t91;
t125 = -MDP(22) - t101;
t113 = MDP(26) * t91 + MDP(27) * t94;
t124 = MDP(29) * t91 + MDP(30) * t94;
t116 = cos(qJ(5));
t92 = sin(qJ(5));
t93 = sin(qJ(4));
t95 = cos(qJ(4));
t72 = -t116 * t95 + t92 * t93;
t122 = t72 * MDP(23);
t121 = t93 * MDP(15) + t95 * MDP(16) + MDP(8);
t120 = t72 ^ 2;
t73 = t116 * t93 + t92 * t95;
t70 = t73 ^ 2;
t88 = -pkin(7) + qJ(2);
t117 = -pkin(8) + t88;
t75 = t117 * t93;
t76 = t117 * t95;
t60 = -t116 * t76 + t75 * t92;
t57 = t60 * t91;
t115 = t60 * t94;
t114 = t91 * t94;
t89 = pkin(1) + qJ(3);
t78 = pkin(4) * t93 + t89;
t107 = MDP(25) * t114;
t86 = t91 ^ 2;
t106 = MDP(24) * t86 + MDP(21) + 0.2e1 * t107;
t105 = pkin(5) * t72 - pkin(9) * t73;
t80 = pkin(4) * t92 + pkin(9);
t81 = -pkin(4) * t116 - pkin(5);
t104 = -t72 * t81 - t73 * t80;
t103 = MDP(15) * t95 - MDP(16) * t93;
t102 = -MDP(26) * t94 + MDP(27) * t91;
t99 = -t73 * MDP(23) + t125 * t72;
t61 = t116 * t75 + t76 * t92;
t87 = t94 ^ 2;
t98 = -t60 * MDP(22) - t61 * MDP(23) + (-MDP(20) + t113) * t73 + ((t86 - t87) * MDP(25) - MDP(24) * t114 - MDP(19)) * t72;
t97 = (MDP(22) * t116 - MDP(23) * t92) * pkin(4);
t96 = (qJ(2) ^ 2);
t55 = pkin(5) * t73 + pkin(9) * t72 + t78;
t54 = t55 * t91 + t61 * t94;
t53 = t55 * t94 - t61 * t91;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2 + t96) * MDP(6)) + (t89 ^ 2 + t96) * MDP(9) + t70 * MDP(28) + (MDP(10) * t95 - 0.2e1 * MDP(11) * t93) * t95 + (MDP(24) * t87 + MDP(17) - 0.2e1 * t107) * t120 + 0.2e1 * (t53 * t73 - t57 * t72) * MDP(29) + 0.2e1 * (-t115 * t72 - t54 * t73) * MDP(30) + 0.2e1 * (MDP(22) * t73 - t122) * t78 + (2 * (MDP(5) + MDP(7)) * qJ(2)) + 0.2e1 * (MDP(18) + t102) * t73 * t72 + 0.2e1 * t121 * t89; -(pkin(1) * MDP(6)) - t89 * MDP(9) + t125 * t73 + MDP(4) - t121 + t122; MDP(6) + MDP(9); qJ(2) * MDP(9) + MDP(7) + t124 * (-t70 - t120); 0; MDP(9); t95 * MDP(12) - t93 * MDP(13) + (t104 * t91 - t115) * MDP(29) + (t104 * t94 + t57) * MDP(30) + t103 * t88 + t98; 0; t103 + t99; -0.2e1 * t101 * t81 + MDP(14) + t106 + 0.2e1 * t97; (t105 * t91 - t115) * MDP(29) + (t105 * t94 + t57) * MDP(30) + t98; 0; t99; t106 + t97 + t101 * (pkin(5) - t81); 0.2e1 * pkin(5) * t101 + t106; t73 * MDP(28) + t53 * MDP(29) - t54 * MDP(30) + t102 * t72; -t101; -t124 * t73; -t124 * t80 + t113; -pkin(9) * t124 + t113; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
