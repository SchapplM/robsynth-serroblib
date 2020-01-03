% Calculate joint inertia matrix for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR13_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR13_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RPRRR13_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:15:29
% EndTime: 2019-12-31 19:15:30
% DurationCPUTime: 0.33s
% Computational Cost: add. (273->106), mult. (531->150), div. (0->0), fcn. (497->6), ass. (0->54)
t94 = sin(qJ(5));
t95 = sin(qJ(4));
t97 = cos(qJ(5));
t98 = cos(qJ(4));
t83 = t94 * t98 + t95 * t97;
t99 = cos(qJ(3));
t75 = t83 * t99;
t82 = t94 * t95 - t97 * t98;
t77 = t82 * t99;
t129 = -t77 * MDP(23) - t75 * MDP(24);
t124 = pkin(7) + pkin(8);
t85 = t124 * t95;
t86 = t124 * t98;
t106 = t83 * MDP(23) - t82 * MDP(24) + (-t85 * t97 - t86 * t94) * MDP(26) - (-t85 * t94 + t86 * t97) * MDP(27);
t127 = t95 * MDP(19) + t98 * MDP(20);
t128 = t95 * MDP(16) + t98 * MDP(17) - pkin(7) * t127 + t106;
t126 = -2 * MDP(22);
t125 = 0.2e1 * MDP(27);
t123 = pkin(8) * t99;
t96 = sin(qJ(3));
t122 = t96 * pkin(4);
t121 = (pkin(1) * MDP(6));
t120 = t95 * t98;
t100 = -pkin(1) - pkin(6);
t115 = t100 * t98;
t108 = t96 * t115;
t84 = pkin(3) * t96 - pkin(7) * t99 + qJ(2);
t63 = t108 + (t84 - t123) * t95;
t119 = t97 * t63;
t74 = t83 * t96;
t76 = t82 * t96;
t118 = -t74 * MDP(26) + t76 * MDP(27);
t116 = t100 * t95;
t114 = t82 * MDP(26);
t113 = t83 * MDP(21);
t110 = MDP(18) + MDP(25);
t109 = t96 * MDP(25) + t129;
t107 = MDP(15) * t120;
t80 = t98 * t84;
t62 = -t98 * t123 + t80 + (pkin(4) - t116) * t96;
t59 = t97 * t62 - t63 * t94;
t105 = t98 * MDP(16) - t95 * MDP(17);
t104 = t98 * MDP(19) - t95 * MDP(20);
t102 = (MDP(26) * t97 - MDP(27) * t94) * pkin(4);
t93 = t99 ^ 2;
t92 = t98 ^ 2;
t91 = t96 ^ 2;
t90 = t95 ^ 2;
t88 = -pkin(4) * t98 - pkin(3);
t81 = (pkin(4) * t95 - t100) * t99;
t69 = t95 * t84 + t108;
t68 = -t96 * t116 + t80;
t60 = t62 * t94 + t119;
t1 = [MDP(1) + t110 * t91 - (-t77 * MDP(21) + t75 * t126) * t77 + ((-2 * MDP(4) + t121) * pkin(1)) + (0.2e1 * t99 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + (t92 * MDP(14) + MDP(7) - 0.2e1 * t107) * t93 + 0.2e1 * (qJ(2) * MDP(12) + (-MDP(8) + t105) * t99 + t129) * t96 + 0.2e1 * (-t93 * t116 + t68 * t96) * MDP(19) + 0.2e1 * (-t93 * t115 - t69 * t96) * MDP(20) + 0.2e1 * (t59 * t96 + t75 * t81) * MDP(26) + (-t60 * t96 - t77 * t81) * t125; MDP(4) - t121 + (-t74 * t96 - t75 * t99) * MDP(26) + (t76 * t96 + t77 * t99) * MDP(27) + t127 * (-t91 - t93); MDP(6); -t77 * t113 + (-t75 * t83 + t77 * t82) * MDP(22) + (t75 * t88 + t81 * t82) * MDP(26) + (-t77 * t88 + t81 * t83) * MDP(27) + (-t100 * MDP(13) - MDP(10) + t128) * t96 + (MDP(9) + t100 * MDP(12) + MDP(14) * t120 + (-t90 + t92) * MDP(15) + (-pkin(3) * t95 + t115) * MDP(19) + (-pkin(3) * t98 - t116) * MDP(20)) * t99; -t96 * MDP(13) + (-t83 * MDP(27) + MDP(12) + t104 - t114) * t99; 0.2e1 * t107 + 0.2e1 * t88 * t114 + MDP(14) * t90 + MDP(11) + 0.2e1 * t104 * pkin(3) + (t88 * t125 + t82 * t126 + t113) * t83; t96 * MDP(18) + t68 * MDP(19) - t69 * MDP(20) + (t97 * t122 + t59) * MDP(26) + (-t119 + (-t62 - t122) * t94) * MDP(27) + t105 * t99 + t109; -t127 * t96 + t118; t128; 0.2e1 * t102 + t110; t59 * MDP(26) - t60 * MDP(27) + t109; t118; t106; MDP(25) + t102; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
