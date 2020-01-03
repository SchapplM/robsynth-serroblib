% Calculate joint inertia matrix for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR8_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:27
% EndTime: 2019-12-31 21:20:28
% DurationCPUTime: 0.32s
% Computational Cost: add. (348->108), mult. (631->146), div. (0->0), fcn. (597->6), ass. (0->54)
t130 = MDP(16) - MDP(19);
t129 = -MDP(17) + MDP(20);
t102 = cos(qJ(2));
t128 = 0.2e1 * t102;
t127 = -2 * MDP(19);
t126 = 2 * MDP(20);
t125 = pkin(3) + pkin(8);
t124 = -pkin(7) - pkin(6);
t101 = cos(qJ(3));
t99 = sin(qJ(2));
t84 = t124 * t99;
t85 = t124 * t102;
t98 = sin(qJ(3));
t74 = -t101 * t85 + t84 * t98;
t82 = -t101 * t102 + t98 * t99;
t66 = -pkin(4) * t82 + t74;
t97 = sin(qJ(5));
t63 = t66 * t97;
t88 = pkin(2) * t98 + qJ(4);
t123 = t82 * t88;
t122 = (MDP(21) * pkin(3));
t121 = qJ(4) * t82;
t100 = cos(qJ(5));
t120 = t100 * t97;
t64 = t66 * t100;
t83 = t101 * t99 + t102 * t98;
t118 = t83 * MDP(25);
t91 = -pkin(2) * t101 - pkin(3);
t117 = t91 * MDP(21);
t116 = t97 * MDP(27);
t94 = t100 * MDP(24);
t115 = t100 * MDP(28);
t114 = 0.2e1 * t83;
t112 = MDP(23) * t120;
t96 = t100 ^ 2;
t113 = t96 * MDP(22) + MDP(15) - 0.2e1 * t112;
t92 = -pkin(2) * t102 - pkin(1);
t111 = -t97 * MDP(25) + t94;
t87 = -pkin(8) + t91;
t110 = -t83 * t87 + t123;
t73 = -t101 * t84 - t85 * t98;
t109 = t125 * t83 + t121;
t108 = MDP(24) * t97 + MDP(25) * t100;
t107 = MDP(27) * t100 - MDP(28) * t97;
t106 = -qJ(4) * t83 + t92;
t105 = t126 + 0.2e1 * t115 + 0.2e1 * t116;
t95 = t97 ^ 2;
t104 = t64 * MDP(28) + (t94 + MDP(13)) * t83 + t129 * t74 - t130 * t73 + ((-t95 + t96) * MDP(23) + MDP(22) * t120 - MDP(14)) * t82;
t68 = pkin(3) * t82 + t106;
t65 = pkin(4) * t83 + t73;
t62 = t125 * t82 + t106;
t61 = t100 * t62 + t65 * t97;
t60 = t100 * t65 - t62 * t97;
t1 = [MDP(1) + pkin(1) * MDP(9) * t128 + (t68 ^ 2 + t73 ^ 2 + t74 ^ 2) * MDP(21) + (t92 * MDP(17) - t68 * MDP(20)) * t114 + (MDP(11) + MDP(26)) * t83 ^ 2 + (t95 * MDP(22) + 0.2e1 * t112) * t82 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t99 + MDP(5) * t128) * t99 + (0.2e1 * t92 * MDP(16) + t68 * t127 + (-MDP(12) + t108) * t114) * t82 + 0.2e1 * (t73 * t83 - t74 * t82) * MDP(18) + 0.2e1 * (t60 * t83 - t82 * t64) * MDP(27) + 0.2e1 * (-t61 * t83 + t82 * t63) * MDP(28); t104 + (-t110 * t100 + t63) * MDP(27) + t99 * MDP(6) + t102 * MDP(7) + (t83 * t91 - t123) * MDP(18) + (t73 * t91 + t74 * t88) * MDP(21) + (-t102 * MDP(10) - t99 * MDP(9)) * pkin(6) + (t110 * MDP(28) - t118) * t97; MDP(8) + ((2 * MDP(19)) + t117) * t91 + 0.2e1 * (MDP(16) * t101 - MDP(17) * t98) * pkin(2) + (t88 * MDP(21) + t105) * t88 + t113; t104 + (-t109 * t100 + t63) * MDP(27) + (-pkin(3) * t83 - t121) * MDP(18) + (-pkin(3) * t73 + qJ(4) * t74) * MDP(21) + (t109 * MDP(28) - t118) * t97; pkin(3) * t127 + qJ(4) * t126 + (-pkin(3) * t91 + qJ(4) * t88) * MDP(21) + (t130 * t101 + t129 * t98) * pkin(2) + t113 + (t116 + t115) * (qJ(4) + t88); (t127 + t122) * pkin(3) + (qJ(4) * MDP(21) + t105) * qJ(4) + t113; t73 * MDP(21) + (MDP(18) + t107) * t83; MDP(19) + t117; MDP(19) - t122; MDP(21); MDP(26) * t83 + t60 * MDP(27) - t61 * MDP(28) + t108 * t82; t107 * t87 + t111; -t107 * t125 + t111; t107; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
