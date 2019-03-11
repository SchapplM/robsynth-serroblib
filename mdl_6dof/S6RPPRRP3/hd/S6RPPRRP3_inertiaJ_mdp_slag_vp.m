% Calculate joint inertia matrix for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPRRP3_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:03:50
% EndTime: 2019-03-09 02:03:51
% DurationCPUTime: 0.34s
% Computational Cost: add. (368->117), mult. (610->166), div. (0->0), fcn. (478->6), ass. (0->50)
t87 = cos(qJ(4));
t121 = 0.2e1 * t87;
t109 = pkin(8) * MDP(25);
t120 = MDP(23) + t109;
t104 = MDP(21) - MDP(24);
t105 = MDP(20) + MDP(22);
t84 = sin(qJ(5));
t86 = cos(qJ(5));
t119 = -t104 * t86 - t105 * t84;
t98 = -pkin(5) * t86 - qJ(6) * t84;
t72 = -pkin(4) + t98;
t107 = t72 * MDP(25);
t118 = t104 * t84 - t105 * t86 - MDP(13) + t107;
t117 = 0.2e1 * pkin(5);
t82 = sin(pkin(9));
t74 = pkin(1) * t82 + qJ(3);
t116 = 0.2e1 * t74;
t83 = cos(pkin(9));
t76 = -pkin(1) * t83 - pkin(2);
t73 = -pkin(7) + t76;
t114 = t73 * t84;
t113 = t73 * t86;
t85 = sin(qJ(4));
t112 = t84 * t85;
t68 = pkin(4) * t85 - pkin(8) * t87 + t74;
t64 = t113 * t85 + t68 * t84;
t78 = t84 ^ 2;
t80 = t86 ^ 2;
t111 = t78 + t80;
t108 = MDP(25) * t87;
t106 = t86 * MDP(16);
t103 = t111 * MDP(23);
t102 = -MDP(25) * pkin(5) - MDP(22);
t99 = 0.2e1 * qJ(6) * MDP(24) + MDP(19);
t97 = -pkin(5) * t84 + qJ(6) * t86;
t61 = qJ(6) * t85 + t64;
t67 = t86 * t68;
t62 = -t67 + (-pkin(5) + t114) * t85;
t96 = t61 * t86 + t62 * t84;
t95 = t86 * MDP(17) - t84 * MDP(18);
t94 = t84 * MDP(17) + t86 * MDP(18);
t93 = (-t112 * t73 + t67) * MDP(20) - t64 * MDP(21);
t92 = t84 * MDP(22) - t86 * MDP(24);
t91 = t109 * t111 - MDP(14);
t90 = (MDP(25) * qJ(6) - t104) * t86 + (-MDP(20) + t102) * t84;
t81 = t87 ^ 2;
t79 = t85 ^ 2;
t77 = t80 * t87;
t65 = (-t73 - t97) * t87;
t1 = [MDP(1) + (t74 ^ 2 + t76 ^ 2) * MDP(7) + t79 * MDP(19) + (t61 ^ 2 + t62 ^ 2 + t65 ^ 2) * MDP(25) + (t82 ^ 2 + t83 ^ 2) * MDP(4) * pkin(1) ^ 2 + (MDP(13) * t116 + (-MDP(9) + t95) * t121) * t85 + 0.2e1 * t76 * MDP(5) + MDP(6) * t116 + 0.2e1 * (-t62 * MDP(22) + t61 * MDP(24) + t93) * t85 + (t80 * MDP(15) - 0.2e1 * t84 * t106 + MDP(8) + 0.2e1 * (-t84 * MDP(20) - t86 * MDP(21)) * t73) * t81 + (t74 * MDP(14) + (-t61 * t84 + t62 * t86) * MDP(23) + t92 * t65) * t121; (t65 * t85 + t87 * t96) * MDP(25); MDP(4) + MDP(7) + (t111 * t81 + t79) * MDP(25); MDP(5) + t76 * MDP(7) + (-t65 * t87 + t85 * t96) * MDP(25) + t119 * (t79 + t81); (-0.1e1 + t111) * t85 * t108; MDP(7) + (t111 * t79 + t81) * MDP(25); t77 * MDP(16) + (-MDP(22) * t86 - t84 * MDP(24) + t107) * t65 + (-t73 * MDP(14) + pkin(8) * t119 - MDP(11) + t94) * t85 + (MDP(10) + t73 * MDP(13) + t86 * t84 * MDP(15) - t78 * MDP(16) + (-pkin(4) * t84 + t113) * MDP(20) + (-pkin(4) * t86 - t114) * MDP(21) + t92 * t72) * t87 + t120 * t96; t77 * MDP(23) + (t78 * MDP(23) + t91) * t87 + t118 * t85; (t103 + t91) * t85 - t118 * t87; MDP(12) + t78 * MDP(15) + (pkin(8) ^ 2 * t111 + t72 ^ 2) * MDP(25) + 0.2e1 * pkin(8) * t103 + 0.2e1 * (MDP(20) * pkin(4) - MDP(22) * t72) * t86 + 0.2e1 * (-MDP(21) * pkin(4) - MDP(24) * t72 + t106) * t84; t67 * MDP(22) + t64 * MDP(24) + (-pkin(5) * t62 + qJ(6) * t61) * MDP(25) + ((t117 - t114) * MDP(22) + t99) * t85 + (MDP(23) * t98 + t95) * t87 + t93; t90 * t87; t90 * t85; MDP(23) * t97 + pkin(8) * t90 + t94; MDP(22) * t117 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(25) + t99; MDP(23) * t86 * t87 - MDP(22) * t85 + MDP(25) * t62; t84 * t108; MDP(25) * t112; t120 * t84; t102; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
