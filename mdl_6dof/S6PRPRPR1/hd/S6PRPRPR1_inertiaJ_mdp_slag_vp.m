% Calculate joint inertia matrix for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S6PRPRPR1_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:27:58
% EndTime: 2019-03-08 19:27:59
% DurationCPUTime: 0.31s
% Computational Cost: add. (397->103), mult. (855->174), div. (0->0), fcn. (971->12), ass. (0->53)
t96 = sin(qJ(6));
t99 = cos(qJ(6));
t110 = MDP(20) * t99 - MDP(21) * t96;
t120 = cos(qJ(4));
t90 = sin(pkin(12));
t93 = cos(pkin(12));
t97 = sin(qJ(4));
t79 = -t93 * t120 + t90 * t97;
t122 = t110 * t79;
t100 = cos(qJ(2));
t91 = sin(pkin(11));
t92 = sin(pkin(6));
t94 = cos(pkin(11));
t98 = sin(qJ(2));
t73 = (-t100 * t94 + t91 * t98) * t92;
t121 = t73 ^ 2;
t81 = t90 * t120 + t93 * t97;
t119 = t81 * t96;
t118 = t81 * t99;
t117 = t96 * t99;
t116 = MDP(14) * pkin(4);
t115 = MDP(19) * t79;
t87 = -t94 * pkin(2) - pkin(3);
t82 = -t120 * pkin(4) + t87;
t114 = t82 * MDP(14);
t113 = MDP(16) * t117;
t112 = pkin(2) * t91 + pkin(8);
t111 = t120 * MDP(11);
t109 = MDP(20) * t96 + MDP(21) * t99;
t108 = (-qJ(5) - t112) * t97;
t107 = t120 * t112;
t75 = (t100 * t91 + t94 * t98) * t92;
t95 = cos(pkin(6));
t106 = t95 * t120 - t75 * t97;
t105 = t97 * MDP(12) - t111;
t104 = (MDP(17) * t99 - MDP(18) * t96) * t81;
t103 = t96 * MDP(17) + t99 * MDP(18) - t109 * (pkin(4) * t90 + pkin(9));
t89 = t99 ^ 2;
t88 = t96 ^ 2;
t86 = -pkin(4) * t93 - pkin(5);
t78 = t81 ^ 2;
t77 = t120 * qJ(5) + t107;
t72 = t75 * t120 + t95 * t97;
t70 = t90 * t108 + t93 * t77;
t68 = -t93 * t108 + t77 * t90;
t67 = t79 * pkin(5) - t81 * pkin(9) + t82;
t66 = t90 * t106 + t93 * t72;
t64 = -t93 * t106 + t72 * t90;
t63 = t67 * t96 + t70 * t99;
t62 = t67 * t99 - t70 * t96;
t61 = t66 * t99 + t73 * t96;
t60 = -t66 * t96 + t73 * t99;
t1 = [MDP(1) + (t75 ^ 2 + t95 ^ 2 + t121) * MDP(5) + (t64 ^ 2 + t66 ^ 2 + t121) * MDP(14); (t64 * t81 - t66 * t79) * MDP(13) + (t64 * t68 + t66 * t70) * MDP(14) + (t64 * t119 + t60 * t79) * MDP(20) + (t64 * t118 - t61 * t79) * MDP(21) + (MDP(3) * t100 - MDP(4) * t98) * t92 + (t105 + t114) * t73 + (-t73 * t94 + t75 * t91) * MDP(5) * pkin(2); MDP(2) - 0.2e1 * t87 * t111 + (t68 ^ 2 + t70 ^ 2 + t82 ^ 2) * MDP(14) + (t89 * MDP(15) - 0.2e1 * t113) * t78 + (t91 ^ 2 + t94 ^ 2) * MDP(5) * pkin(2) ^ 2 + (0.2e1 * t87 * MDP(12) + MDP(6) * t97 + 0.2e1 * t120 * MDP(7)) * t97 + (0.2e1 * t104 + t115) * t79 + 0.2e1 * (t68 * t81 - t70 * t79) * MDP(13) + 0.2e1 * (t68 * t119 + t62 * t79) * MDP(20) + 0.2e1 * (t68 * t118 - t63 * t79) * MDP(21); t95 * MDP(5) + (t64 * t79 + t66 * t81) * MDP(14); (t68 * t79 + t70 * t81) * MDP(14); MDP(5) + (t79 ^ 2 + t78) * MDP(14); t106 * MDP(11) - t72 * MDP(12) - t110 * t64 + (-t64 * t93 + t66 * t90) * t116; t120 * MDP(9) - MDP(12) * t107 + (-t112 * MDP(11) + MDP(8)) * t97 - t110 * t68 + t103 * t79 + (MDP(15) * t117 + (-t88 + t89) * MDP(16) + t109 * t86) * t81 + ((-t79 * t90 - t81 * t93) * MDP(13) + (-t68 * t93 + t70 * t90) * MDP(14)) * pkin(4); -t122 + (-t79 * t93 + t81 * t90) * t116 - t105; 0.2e1 * t113 + t88 * MDP(15) + MDP(10) - 0.2e1 * t110 * t86 + (t90 ^ 2 + t93 ^ 2) * MDP(14) * pkin(4) ^ 2; t73 * MDP(14); t114 + t122; 0; 0; MDP(14); MDP(20) * t60 - MDP(21) * t61; t62 * MDP(20) - t63 * MDP(21) + t104 + t115; -t109 * t81; t103; t110; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
