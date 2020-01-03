% Calculate joint inertia matrix for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRPRP3_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:15
% EndTime: 2019-12-31 19:51:16
% DurationCPUTime: 0.30s
% Computational Cost: add. (384->89), mult. (649->110), div. (0->0), fcn. (604->6), ass. (0->56)
t98 = sin(pkin(8));
t99 = cos(pkin(8));
t123 = t98 ^ 2 + t99 ^ 2;
t124 = t123 * qJ(3);
t100 = sin(qJ(4));
t121 = t100 * t98;
t131 = cos(qJ(4));
t81 = -t131 * t99 + t121;
t114 = t131 * t98;
t82 = t100 * t99 + t114;
t139 = t81 * MDP(16) + t82 * MDP(17);
t138 = t81 * MDP(18) - t82 * MDP(20);
t137 = t99 * MDP(7) - t98 * MDP(8);
t136 = 2 * MDP(9);
t135 = 2 * MDP(19);
t134 = t99 * pkin(3);
t101 = sin(qJ(2));
t89 = t101 * pkin(1) + qJ(3);
t132 = -pkin(7) - t89;
t102 = cos(qJ(2));
t130 = t102 * pkin(1);
t95 = t99 * pkin(7);
t79 = t99 * t89 + t95;
t63 = t100 * t79 - t132 * t114;
t64 = t132 * t121 + t131 * t79;
t129 = t63 * t82 - t64 * t81;
t113 = (-pkin(7) - qJ(3)) * t98;
t84 = t99 * qJ(3) + t95;
t68 = t100 * t84 - t131 * t113;
t69 = t100 * t113 + t131 * t84;
t128 = t68 * t82 - t69 * t81;
t90 = -pkin(2) - t134;
t65 = t81 * pkin(4) - t82 * qJ(5) + t90;
t62 = t65 - t130;
t127 = t62 + t65;
t91 = -pkin(2) - t130;
t83 = t91 - t134;
t126 = t83 + t90;
t125 = t123 * t89;
t122 = pkin(2) * MDP(10);
t119 = t62 * MDP(21);
t118 = t65 * MDP(21);
t117 = t91 * MDP(10);
t116 = (-pkin(4) * t82 - t81 * qJ(5)) * MDP(19) + t82 * MDP(13) - t81 * MDP(14);
t115 = MDP(4) + (MDP(11) * t82 - 0.2e1 * MDP(12) * t81) * t82;
t112 = t123 * MDP(10);
t111 = -MDP(21) * pkin(4) - MDP(18);
t110 = -MDP(16) + t111;
t109 = 0.2e1 * t137;
t108 = MDP(21) * qJ(5) - MDP(17) + MDP(20);
t107 = 0.2e1 * t138;
t106 = (t102 * MDP(5) - t101 * MDP(6)) * pkin(1);
t105 = -t137 + t138 + t139;
t104 = 0.2e1 * t139;
t75 = t82 * MDP(19);
t1 = [MDP(1) + (t63 ^ 2 + t64 ^ 2) * MDP(21) + t89 ^ 2 * t112 + (-t109 + t117) * t91 + t83 * t104 + (t107 + t119) * t62 + 0.2e1 * t106 + t125 * t136 + t129 * t135 + t115; (t124 + t125) * MDP(9) + (-t91 * pkin(2) + t89 * t124) * MDP(10) + (t128 + t129) * MDP(19) + (t62 * t65 + t63 * t68 + t64 * t69) * MDP(21) + t106 + (t126 * MDP(17) - t127 * MDP(20)) * t82 + (t126 * MDP(16) + t127 * MDP(18)) * t81 + t115 + t137 * (pkin(2) - t91); (t68 ^ 2 + t69 ^ 2) * MDP(21) + t90 * t104 + (t107 + t118) * t65 + qJ(3) ^ 2 * t112 + (t109 + t122) * pkin(2) + t124 * t136 + t128 * t135 + t115; t105 + t117 + t119; t105 + t118 - t122; MDP(10) + MDP(21); t108 * t64 + t110 * t63 + t116; t108 * t69 + t110 * t68 + t116; 0; MDP(15) + 0.2e1 * pkin(4) * MDP(18) + 0.2e1 * qJ(5) * MDP(20) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(21); t63 * MDP(21) + t75; t68 * MDP(21) + t75; 0; t111; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
