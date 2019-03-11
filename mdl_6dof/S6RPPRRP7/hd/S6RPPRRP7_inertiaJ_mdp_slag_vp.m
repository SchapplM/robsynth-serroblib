% Calculate joint inertia matrix for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRRP7_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:13:55
% EndTime: 2019-03-09 02:13:56
% DurationCPUTime: 0.29s
% Computational Cost: add. (490->108), mult. (824->154), div. (0->0), fcn. (846->6), ass. (0->56)
t119 = cos(qJ(4));
t88 = sin(pkin(9));
t89 = cos(pkin(9));
t92 = sin(qJ(4));
t71 = -t119 * t89 + t92 * t88;
t122 = t71 ^ 2;
t72 = t119 * t88 + t92 * t89;
t68 = t72 ^ 2;
t123 = -t68 - t122;
t90 = -pkin(1) - qJ(3);
t120 = -pkin(7) + t90;
t91 = sin(qJ(5));
t93 = cos(qJ(5));
t118 = t91 * t93;
t74 = t120 * t88;
t75 = t120 * t89;
t67 = t119 * t74 + t92 * t75;
t117 = t93 * t67;
t116 = -qJ(6) - pkin(8);
t78 = t88 ^ 2 + t89 ^ 2;
t86 = t91 ^ 2;
t87 = t93 ^ 2;
t115 = t86 + t87;
t114 = MDP(26) * pkin(5);
t113 = qJ(6) * t71;
t79 = t88 * pkin(3) + qJ(2);
t112 = MDP(21) * t91;
t65 = t72 * pkin(4) + t71 * pkin(8) + t79;
t61 = t93 * t65 - t91 * t67;
t59 = t72 * pkin(5) + t93 * t113 + t61;
t111 = t59 * MDP(26);
t110 = t71 * MDP(26);
t81 = -t93 * pkin(5) - pkin(4);
t109 = t81 * MDP(26);
t108 = t91 * MDP(24);
t107 = t93 * MDP(24);
t106 = MDP(19) * t118;
t105 = MDP(25) * pkin(5) - MDP(20);
t104 = MDP(23) + t114;
t60 = t117 + (t65 + t113) * t91;
t103 = t59 * t93 + t60 * t91;
t76 = t116 * t91;
t77 = t116 * t93;
t102 = t93 * t76 - t91 * t77;
t101 = -t76 * t91 - t77 * t93;
t100 = t88 * MDP(7) + t89 * MDP(8);
t99 = t115 * MDP(25) - MDP(17);
t98 = t61 * MDP(23) - (t91 * t65 + t117) * MDP(24);
t97 = t93 * MDP(23) - t108;
t96 = t91 * MDP(23) + t107;
t66 = -t119 * t75 + t92 * t74;
t95 = MDP(16) + t97;
t94 = qJ(2) ^ 2;
t70 = t78 * t90;
t63 = -t91 * t71 * pkin(5) + t66;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2) + t94) * MDP(6) + (t78 * t90 ^ 2 + t94) * MDP(10) + t68 * MDP(22) + (t59 ^ 2 + t60 ^ 2 + t63 ^ 2) * MDP(26) + (t87 * MDP(18) + MDP(11) - 0.2e1 * t106) * t122 - 0.2e1 * t79 * t71 * MDP(17) - 0.2e1 * t70 * MDP(9) + 0.2e1 * (t103 * MDP(25) - t96 * t66) * t71 + 0.2e1 * (MDP(5) + t100) * qJ(2) + 0.2e1 * ((-MDP(20) * t93 + MDP(12) + t112) * t71 + t79 * MDP(16) + t98) * t72; t63 * t110 + t70 * MDP(10) - pkin(1) * MDP(6) - t78 * MDP(9) + MDP(4) + (t60 * t72 * MDP(26) + t123 * MDP(24)) * t93 + (t123 * MDP(23) - t72 * t111) * t91; MDP(6) + t78 * MDP(10) + (t115 * t68 + t122) * MDP(26); qJ(2) * MDP(10) + t103 * MDP(26) + t99 * t71 + t95 * t72 + t100; 0; t115 * MDP(26) + MDP(10); -t67 * MDP(17) + (-t59 * t91 + t60 * t93) * MDP(25) + (t59 * t76 - t60 * t77 + t63 * t81) * MDP(26) - t95 * t66 + (t91 * MDP(20) + t93 * MDP(21) - t96 * pkin(8) - MDP(14)) * t72 + (-MDP(13) - MDP(18) * t118 + (t86 - t87) * MDP(19) + t102 * MDP(25) + t96 * pkin(4)) * t71; (-t95 + t109) * t71 + (t101 * MDP(26) + t99) * t72; t102 * MDP(26); MDP(15) + t86 * MDP(18) + 0.2e1 * t106 + 0.2e1 * t101 * MDP(25) + (t76 ^ 2 + t77 ^ 2 + t81 ^ 2) * MDP(26) + 0.2e1 * t97 * pkin(4); pkin(5) * t111 + t72 * MDP(22) + (t105 * t93 + t112) * t71 + t98; (-t104 * t91 - t107) * t72; t104 * t93 - t108; t76 * t114 + (-MDP(24) * pkin(8) + MDP(21)) * t93 + (-MDP(23) * pkin(8) - t105) * t91; MDP(26) * pkin(5) ^ 2 + MDP(22); t63 * MDP(26); t110; 0; t109; 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
