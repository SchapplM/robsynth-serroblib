% Calculate joint inertia matrix for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP8_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:25
% EndTime: 2019-12-31 21:09:27
% DurationCPUTime: 0.48s
% Computational Cost: add. (387->147), mult. (685->188), div. (0->0), fcn. (532->4), ass. (0->47)
t92 = sin(qJ(2));
t123 = 0.2e1 * t92;
t91 = sin(qJ(3));
t111 = t91 * qJ(4);
t89 = pkin(3) + qJ(5);
t93 = cos(qJ(3));
t122 = -t89 * t93 - t111;
t121 = MDP(20) + MDP(23);
t120 = pkin(4) + pkin(7);
t94 = cos(qJ(2));
t119 = pkin(6) * t94;
t117 = t91 * t92;
t116 = t92 * t93;
t115 = pkin(3) * t117 + t92 * pkin(6);
t76 = -t94 * pkin(2) - t92 * pkin(7) - pkin(1);
t114 = t91 * t119 - t93 * t76;
t70 = t93 * t119 + t91 * t76;
t86 = t91 ^ 2;
t88 = t93 ^ 2;
t113 = t86 + t88;
t112 = MDP(21) * pkin(7);
t83 = qJ(4) * t93;
t110 = t94 * qJ(4);
t109 = t93 * MDP(12);
t108 = MDP(17) - MDP(20);
t107 = MDP(18) + MDP(22);
t85 = t94 * pkin(3);
t68 = t85 + t114;
t105 = -pkin(3) * MDP(21) + MDP(19);
t67 = t110 - t70;
t104 = -t93 * pkin(3) - t111;
t102 = t93 * MDP(13) - t91 * MDP(14);
t101 = -MDP(16) * t114 - t70 * MDP(17);
t100 = -t91 * MDP(23) - t93 * MDP(24);
t72 = -pkin(2) + t122;
t75 = -pkin(2) + t104;
t78 = t120 * t93;
t99 = pkin(2) * MDP(16) + t75 * MDP(19) + t78 * MDP(22) - t72 * MDP(24);
t77 = t120 * t91;
t98 = -pkin(2) * MDP(17) - t75 * MDP(20) + t77 * MDP(22) - t72 * MDP(23);
t97 = -t91 * MDP(13) - t93 * MDP(14) - t78 * MDP(23) + t77 * MDP(24);
t95 = qJ(4) ^ 2;
t71 = -t92 * t83 + t115;
t66 = (qJ(5) * t91 - t83) * t92 + t115;
t65 = -pkin(4) * t117 - t67;
t64 = pkin(4) * t116 + t94 * qJ(5) + t68;
t1 = [MDP(1) + (t67 ^ 2 + t68 ^ 2 + t71 ^ 2) * MDP(21) + (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) * MDP(25) + (t94 * MDP(15) + 0.2e1 * pkin(1) * MDP(9) + (MDP(5) - t102) * t123) * t94 + 0.2e1 * (-t68 * MDP(19) + t67 * MDP(20) - t65 * MDP(23) + t64 * MDP(24) - t101) * t94 + ((t67 * t91 + t68 * t93) * MDP(18) + (t64 * t93 - t65 * t91) * MDP(22) + (-t91 * MDP(19) - t93 * MDP(20)) * t71 + (-t93 * MDP(23) + t91 * MDP(24)) * t66) * t123 + (-0.2e1 * pkin(1) * MDP(10) + (t88 * MDP(11) - 0.2e1 * t91 * t109 + MDP(4) + 0.2e1 * (t91 * MDP(16) + t93 * MDP(17)) * pkin(6)) * t92) * t92; (t64 * t91 + t65 * t93) * MDP(22) + (t64 * t77 + t65 * t78) * MDP(25) + (t93 * MDP(19) - t91 * MDP(20) + MDP(21) * t75) * t71 + (MDP(25) * t72 + t100) * t66 + (-pkin(6) * MDP(10) + MDP(7) + (t108 * t93 + (MDP(16) - MDP(19)) * t91) * pkin(7) + t97) * t94 + (MDP(6) - pkin(6) * MDP(9) + (-t86 + t88) * MDP(12) + (-pkin(6) * MDP(16) + t98) * t93 + (t93 * MDP(11) + pkin(6) * MDP(17) - t99) * t91) * t92 + (MDP(18) + t112) * (-t67 * t93 + t68 * t91); MDP(8) + t86 * MDP(11) + (t113 * pkin(7) ^ 2 + t75 ^ 2) * MDP(21) + (t72 ^ 2 + t77 ^ 2 + t78 ^ 2) * MDP(25) + 0.2e1 * t113 * MDP(18) * pkin(7) + 0.2e1 * t99 * t93 + 0.2e1 * (t98 + t109) * t91; (0.2e1 * t85 + t114) * MDP(19) + (-t68 * pkin(3) - t67 * qJ(4)) * MDP(21) - t68 * MDP(24) + (t65 * qJ(4) - t64 * t89) * MDP(25) + (-MDP(15) + (-qJ(5) - t89) * MDP(24)) * t94 + (t104 * MDP(18) + t122 * MDP(22) + t100 * pkin(4) + t102) * t92 + t101 + t121 * (-0.2e1 * t110 + t70); (-pkin(3) * t91 + t83) * MDP(18) + (-t89 * t91 + t83) * MDP(22) + (t78 * qJ(4) - t77 * t89) * MDP(25) + ((qJ(4) * MDP(21) - t108) * t93 + (-MDP(16) + t105) * t91) * pkin(7) - t97; MDP(15) - 0.2e1 * pkin(3) * MDP(19) + (pkin(3) ^ 2 + t95) * MDP(21) + 0.2e1 * t89 * MDP(24) + (t89 ^ 2 + t95) * MDP(25) + 0.2e1 * t121 * qJ(4); t68 * MDP(21) + t64 * MDP(25) + (-MDP(19) + MDP(24)) * t94 + t107 * t116; t77 * MDP(25) + (t107 + t112) * t91; -t89 * MDP(25) - MDP(24) + t105; MDP(21) + MDP(25); -MDP(22) * t117 - t94 * MDP(23) + t65 * MDP(25); t93 * MDP(22) + t78 * MDP(25); MDP(25) * qJ(4) + MDP(23); 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
