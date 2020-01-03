% Calculate joint inertia matrix for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRRRP2_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:46
% EndTime: 2020-01-03 12:11:47
% DurationCPUTime: 0.18s
% Computational Cost: add. (339->84), mult. (581->111), div. (0->0), fcn. (549->6), ass. (0->50)
t100 = cos(qJ(4));
t99 = sin(qJ(2));
t89 = t99 * pkin(1) + pkin(7);
t98 = sin(qJ(3));
t79 = (-pkin(8) - t89) * t98;
t101 = cos(qJ(3));
t96 = t101 * pkin(8);
t80 = t101 * t89 + t96;
t97 = sin(qJ(4));
t108 = -t100 * t80 - t97 * t79;
t112 = t100 * t79 - t97 * t80;
t132 = t112 * MDP(19) + t108 * MDP(20);
t85 = (-pkin(7) - pkin(8)) * t98;
t86 = t101 * pkin(7) + t96;
t107 = -t100 * t86 - t97 * t85;
t111 = t100 * t85 - t97 * t86;
t131 = t111 * MDP(19) + t107 * MDP(20);
t105 = t101 * MDP(12) - t98 * MDP(13);
t81 = -t100 * t101 + t97 * t98;
t82 = t100 * t98 + t97 * t101;
t130 = t81 * MDP(19) + t82 * MDP(20);
t129 = 2 * MDP(21);
t128 = pkin(3) * t97;
t127 = t81 * pkin(4);
t102 = cos(qJ(2));
t125 = t102 * pkin(1);
t119 = t82 * qJ(5);
t61 = t112 - t119;
t78 = t81 * qJ(5);
t62 = -t108 - t78;
t124 = -t61 * t82 - t62 * t81;
t63 = t111 - t119;
t64 = -t107 - t78;
t123 = -t63 * t82 - t64 * t81;
t122 = t82 * MDP(16) - t81 * MDP(17);
t120 = MDP(22) * pkin(4);
t115 = t100 * MDP(19);
t113 = -pkin(4) * t82 * MDP(21) + t122;
t92 = -t101 * pkin(3) - pkin(2);
t90 = t100 * pkin(3) + pkin(4);
t110 = (-t81 * t128 - t90 * t82) * MDP(21) + t101 * MDP(10) + t98 * MDP(9) + t122;
t109 = MDP(4) + (MDP(7) * t98 + 0.2e1 * MDP(8) * t101) * t98 + (MDP(14) * t82 - 0.2e1 * MDP(15) * t81) * t82;
t106 = -t98 * MDP(12) - t101 * MDP(13);
t84 = t92 - t125;
t104 = (t102 * MDP(5) - t99 * MDP(6)) * pkin(1);
t103 = 0.2e1 * t130;
t91 = -pkin(2) - t125;
t72 = t92 + t127;
t71 = t84 + t127;
t1 = [MDP(1) + t124 * t129 + (t61 ^ 2 + t62 ^ 2 + t71 ^ 2) * MDP(22) + t84 * t103 + t109 - 0.2e1 * t105 * t91 + 0.2e1 * t104; (t123 + t124) * MDP(21) + (t61 * t63 + t62 * t64 + t71 * t72) * MDP(22) + t104 + t109 + t105 * (pkin(2) - t91) + t130 * (t84 + t92); t123 * t129 + (t63 ^ 2 + t64 ^ 2 + t72 ^ 2) * MDP(22) + t92 * t103 + 0.2e1 * t105 * pkin(2) + t109; (t62 * t128 + t61 * t90) * MDP(22) + t106 * t89 + t110 + t132; (t64 * t128 + t63 * t90) * MDP(22) + t106 * pkin(7) + t110 + t131; t90 ^ 2 * MDP(22) + MDP(11) + MDP(18) + (0.2e1 * t115 + (MDP(22) * t128 - 0.2e1 * MDP(20)) * t97) * pkin(3); t61 * t120 + t113 + t132; t63 * t120 + t113 + t131; t90 * t120 + MDP(18) + (-t97 * MDP(20) + t115) * pkin(3); MDP(22) * pkin(4) ^ 2 + MDP(18); t71 * MDP(22); t72 * MDP(22); 0; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
