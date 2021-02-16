% Calculate joint inertia matrix for
% S5RRPRR2
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
%   see S5RRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR2_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:23:36
% EndTime: 2021-01-15 21:23:37
% DurationCPUTime: 0.22s
% Computational Cost: add. (530->90), mult. (974->130), div. (0->0), fcn. (1103->8), ass. (0->51)
t92 = sin(pkin(9));
t93 = cos(pkin(9));
t96 = sin(qJ(2));
t99 = cos(qJ(2));
t83 = t92 * t96 - t93 * t99;
t84 = t92 * t99 + t93 * t96;
t95 = sin(qJ(4));
t98 = cos(qJ(4));
t73 = t98 * t83 + t95 * t84;
t91 = -t99 * pkin(2) - pkin(1);
t78 = t83 * pkin(3) + t91;
t121 = 0.2e1 * t73 * pkin(4) + 0.2e1 * t78;
t120 = 0.2e1 * t78;
t119 = 0.2e1 * t99;
t118 = pkin(2) * t92;
t90 = t93 * pkin(2) + pkin(3);
t81 = t118 * t98 + t95 * t90;
t97 = cos(qJ(5));
t117 = t97 * t81;
t116 = -qJ(3) - pkin(6);
t74 = -t95 * t83 + t98 * t84;
t94 = sin(qJ(5));
t62 = t97 * t73 + t94 * t74;
t115 = t62 * MDP(27);
t80 = -t118 * t95 + t98 * t90;
t79 = pkin(4) + t80;
t66 = t97 * t79 - t94 * t81;
t114 = t66 * MDP(27);
t113 = (-t94 * t79 - t117) * MDP(28);
t112 = t73 * MDP(20);
t111 = t80 * MDP(20);
t110 = t81 * MDP(21);
t109 = t83 * MDP(11);
t108 = t84 * MDP(12);
t107 = t91 * MDP(14);
t106 = MDP(19) + MDP(26);
t86 = t116 * t96;
t87 = t116 * t99;
t75 = t93 * t86 + t92 * t87;
t68 = -t84 * pkin(7) + t75;
t76 = t92 * t86 - t93 * t87;
t69 = -t83 * pkin(7) + t76;
t105 = t98 * t68 - t95 * t69;
t56 = -t74 * pkin(8) + t105;
t103 = -t95 * t68 - t98 * t69;
t57 = -t73 * pkin(8) - t103;
t63 = -t94 * t73 + t97 * t74;
t104 = t63 * MDP(24) - t62 * MDP(25) + (t97 * t56 - t94 * t57) * MDP(27) + (-t94 * t56 - t97 * t57) * MDP(28);
t102 = (t97 * MDP(27) - t94 * MDP(28)) * pkin(4);
t101 = t74 * MDP(17) - t73 * MDP(18) + t105 * MDP(20) + t103 * MDP(21) + t104;
t1 = [MDP(1) + pkin(1) * MDP(9) * t119 + 0.2e1 * (-t75 * t84 - t76 * t83) * MDP(13) + (t75 ^ 2 + t76 ^ 2) * MDP(14) + t112 * t120 + t115 * t121 + (t107 + 0.2e1 * t108 + 0.2e1 * t109) * t91 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t96 + MDP(5) * t119) * t96 + (MDP(15) * t74 - 0.2e1 * t73 * MDP(16) + MDP(21) * t120) * t74 + (MDP(22) * t63 - 0.2e1 * t62 * MDP(23) + MDP(28) * t121) * t63; t75 * MDP(11) - t76 * MDP(12) + t96 * MDP(6) + t99 * MDP(7) + (-t99 * MDP(10) - t96 * MDP(9)) * pkin(6) + ((-t83 * t92 - t84 * t93) * MDP(13) + (t75 * t93 + t76 * t92) * MDP(14)) * pkin(2) + t101; MDP(8) + 0.2e1 * t111 - 0.2e1 * t110 + 0.2e1 * t114 + 0.2e1 * t113 + t106 + (0.2e1 * t93 * MDP(11) - 0.2e1 * t92 * MDP(12) + (t92 ^ 2 + t93 ^ 2) * MDP(14) * pkin(2)) * pkin(2); t74 * MDP(21) + t63 * MDP(28) + t107 + t108 + t109 + t112 + t115; 0; MDP(14); t101; t111 - t110 + (t97 * pkin(4) + t66) * MDP(27) + (-t117 + (-pkin(4) - t79) * t94) * MDP(28) + t106; 0; 0.2e1 * t102 + t106; t104; MDP(26) + t113 + t114; 0; MDP(26) + t102; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
