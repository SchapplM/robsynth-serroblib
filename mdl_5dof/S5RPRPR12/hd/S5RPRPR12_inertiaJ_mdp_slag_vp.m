% Calculate joint inertia matrix for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRPR12_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:16
% EndTime: 2019-12-31 18:30:17
% DurationCPUTime: 0.29s
% Computational Cost: add. (473->99), mult. (917->148), div. (0->0), fcn. (992->8), ass. (0->58)
t91 = sin(pkin(9));
t93 = cos(pkin(9));
t120 = t91 ^ 2 + t93 ^ 2;
t130 = MDP(17) * t120;
t85 = -t93 * pkin(4) - pkin(3);
t129 = 0.2e1 * t85;
t94 = cos(pkin(8));
t86 = -t94 * pkin(2) - pkin(1);
t128 = 0.2e1 * t86;
t127 = -2 * MDP(20);
t126 = cos(qJ(3));
t125 = pkin(1) * MDP(7);
t92 = sin(pkin(8));
t96 = sin(qJ(3));
t78 = t126 * t92 + t96 * t94;
t124 = t78 * t91;
t123 = t78 * t93;
t122 = pkin(6) + qJ(2);
t121 = pkin(7) + qJ(4);
t76 = -t126 * t94 + t96 * t92;
t69 = t76 * pkin(3) - t78 * qJ(4) + t86;
t80 = t122 * t92;
t82 = t122 * t94;
t74 = t126 * t82 - t96 * t80;
t63 = t91 * t69 + t93 * t74;
t118 = pkin(3) * MDP(18);
t117 = t92 * MDP(5);
t116 = t94 * MDP(4);
t95 = sin(qJ(5));
t97 = cos(qJ(5));
t77 = t97 * t91 + t95 * t93;
t65 = t77 * t78;
t115 = t65 * MDP(22);
t75 = t95 * t91 - t97 * t93;
t66 = t75 * t78;
t114 = t66 * MDP(19);
t113 = t66 * MDP(21);
t112 = t75 * MDP(24);
t111 = t76 * MDP(23);
t110 = t91 * MDP(16);
t109 = t93 * MDP(15);
t62 = t93 * t69 - t91 * t74;
t108 = t120 * MDP(18);
t107 = t62 * t93 + t63 * t91;
t106 = -t62 * t91 + t63 * t93;
t105 = MDP(15) * t91 + MDP(16) * t93;
t60 = t76 * pkin(4) - pkin(7) * t123 + t62;
t61 = -pkin(7) * t124 + t63;
t104 = (t97 * t60 - t95 * t61) * MDP(24) - (t95 * t60 + t97 * t61) * MDP(25);
t103 = t65 * MDP(24) - t66 * MDP(25);
t102 = -t77 * MDP(25) - t112;
t72 = t126 * t80 + t96 * t82;
t101 = t102 + t109 - t110;
t79 = t121 * t91;
t81 = t121 * t93;
t100 = t77 * MDP(21) - t75 * MDP(22) + (-t97 * t79 - t95 * t81) * MDP(24) - (-t95 * t79 + t97 * t81) * MDP(25);
t64 = pkin(4) * t124 + t72;
t1 = [MDP(1) + (t62 ^ 2 + t63 ^ 2 + t72 ^ 2) * MDP(18) + (MDP(14) * t128 + MDP(8) * t78) * t78 - (t65 * t127 - t114) * t66 + (0.2e1 * t116 - 0.2e1 * t117 + t125) * pkin(1) + (MDP(13) * t128 - 0.2e1 * t78 * MDP(9) + t111 - 0.2e1 * t113 - 0.2e1 * t115) * t76 + 0.2e1 * t103 * t64 + 0.2e1 * (t62 * MDP(15) - t63 * MDP(16) + t104) * t76 + 0.2e1 * (-MDP(17) * t107 + t105 * t72) * t78 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t92 ^ 2 + t94 ^ 2) * qJ(2); -t116 + t117 - t125 + t107 * MDP(18) + (MDP(14) - t130) * t78 + (MDP(13) + t101) * t76; MDP(7) + t108; t78 * MDP(10) - t72 * MDP(13) - t74 * MDP(14) + (-pkin(3) * t124 - t72 * t93) * MDP(15) + (-pkin(3) * t123 + t72 * t91) * MDP(16) + t106 * MDP(17) + (-t72 * pkin(3) + qJ(4) * t106) * MDP(18) - t77 * t114 + (-t77 * t65 + t66 * t75) * MDP(20) + (t64 * t75 + t85 * t65) * MDP(24) + (t64 * t77 - t85 * t66) * MDP(25) + (-qJ(4) * t105 - MDP(11) + t100) * t76; 0; t112 * t129 + MDP(12) + (0.2e1 * t109 - 0.2e1 * t110 + t118) * pkin(3) + (MDP(19) * t77 + MDP(25) * t129 + t127 * t75) * t77 + (t108 * qJ(4) + 0.2e1 * t130) * qJ(4); t72 * MDP(18) + t105 * t78 + t103; 0; -t101 - t118; MDP(18); t104 + t111 - t113 - t115; t102; t100; 0; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
