% Calculate joint inertia matrix for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPPR2_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:06:32
% EndTime: 2019-03-08 21:06:33
% DurationCPUTime: 0.35s
% Computational Cost: add. (431->119), mult. (873->179), div. (0->0), fcn. (956->10), ass. (0->61)
t140 = MDP(12) + MDP(14);
t124 = MDP(13) + MDP(17);
t100 = cos(pkin(11));
t102 = sin(qJ(3));
t105 = cos(qJ(3));
t98 = sin(pkin(11));
t86 = t100 * t102 + t98 * t105;
t95 = -t105 * pkin(3) - pkin(2);
t112 = -t86 * qJ(5) + t95;
t85 = -t100 * t105 + t98 * t102;
t77 = t85 * pkin(4) + t112;
t139 = -0.2e1 * t77;
t138 = -qJ(4) - pkin(8);
t137 = MDP(13) * pkin(3);
t101 = sin(qJ(6));
t136 = t101 * t85;
t103 = sin(qJ(2));
t99 = sin(pkin(6));
t135 = t103 * t99;
t104 = cos(qJ(6));
t134 = t104 * t85;
t106 = cos(qJ(2));
t133 = t106 * t99;
t132 = cos(pkin(6));
t131 = t85 * MDP(15);
t130 = t86 * MDP(22);
t90 = t98 * pkin(3) + qJ(5);
t129 = t90 * MDP(17);
t94 = -t100 * pkin(3) - pkin(4);
t128 = t94 * MDP(17);
t127 = MDP(10) * t105;
t126 = t101 * MDP(23);
t125 = t104 * MDP(24);
t121 = t138 * t102;
t88 = t138 * t105;
t78 = -t100 * t121 - t98 * t88;
t80 = -t100 * t88 + t98 * t121;
t123 = t78 ^ 2 + t80 ^ 2;
t122 = t104 * t101 * MDP(19);
t120 = MDP(15) + t128;
t89 = -pkin(9) + t94;
t115 = t90 * t85 - t86 * t89;
t114 = t104 * MDP(23) - t101 * MDP(24);
t113 = -t125 - t126;
t111 = -MDP(16) + t113;
t110 = (MDP(20) * t101 + MDP(21) * t104) * t85;
t109 = -t102 * t135 + t132 * t105;
t108 = t95 * MDP(13) + t77 * MDP(17) - t131;
t97 = t104 ^ 2;
t96 = t101 ^ 2;
t83 = t132 * t102 + t105 * t135;
t75 = t100 * t83 + t109 * t98;
t73 = -t100 * t109 + t98 * t83;
t72 = -t85 * pkin(5) + t80;
t71 = t86 * pkin(5) + t78;
t70 = (pkin(4) + pkin(9)) * t85 + t112;
t69 = -t101 * t73 + t104 * t133;
t68 = t101 * t133 + t104 * t73;
t67 = t101 * t71 + t104 * t70;
t66 = -t101 * t70 + t104 * t71;
t1 = [MDP(1) + t124 * (t99 ^ 2 * t106 ^ 2 + t73 ^ 2 + t75 ^ 2); (-t75 * t134 + t68 * t86) * MDP(23) + (t75 * t136 + t69 * t86) * MDP(24) + (-t103 * MDP(4) + (-MDP(11) * t102 + t86 * MDP(16) + MDP(3) - t108 + t127) * t106) * t99 + t124 * (t73 * t78 + t75 * t80) + t140 * (t73 * t86 - t75 * t85); MDP(2) + 0.2e1 * pkin(2) * t127 + (t95 ^ 2 + t123) * MDP(13) + t131 * t139 + (t77 ^ 2 + t123) * MDP(17) + (t96 * MDP(18) + 0.2e1 * t122) * t85 ^ 2 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t102 + 0.2e1 * t105 * MDP(6)) * t102 + (MDP(16) * t139 + 0.2e1 * t110 + t130) * t86 + 0.2e1 * (-t72 * t134 + t66 * t86) * MDP(23) + 0.2e1 * (t72 * t136 - t67 * t86) * MDP(24) + 0.2e1 * t140 * (t78 * t86 - t80 * t85); t109 * MDP(10) - t83 * MDP(11) + (-t100 * t137 + t120) * t73 + (t98 * t137 - t111 + t129) * t75; t102 * MDP(7) + t105 * MDP(8) + t94 * t86 * MDP(14) + t78 * MDP(15) + t80 * MDP(16) + (t78 * t94 + t80 * t90) * MDP(17) + (-t90 * MDP(14) + (-t96 + t97) * MDP(19)) * t85 + (-t102 * MDP(10) - t105 * MDP(11)) * pkin(8) + (t86 * MDP(20) - t115 * MDP(23) + t72 * MDP(24)) * t104 + (MDP(18) * t134 - t86 * MDP(21) + t72 * MDP(23) + t115 * MDP(24)) * t101 + ((-t100 * t86 - t85 * t98) * MDP(12) + (-t100 * t78 + t80 * t98) * MDP(13)) * pkin(3); -0.2e1 * t122 + t97 * MDP(18) + MDP(9) + (0.2e1 * MDP(15) + t128) * t94 + (t100 ^ 2 + t98 ^ 2) * MDP(13) * pkin(3) ^ 2 + (0.2e1 * MDP(16) + 0.2e1 * t125 + 0.2e1 * t126 + t129) * t90; -t124 * t133; t111 * t86 + t108; 0; t124; t73 * MDP(17); t78 * MDP(17) + (MDP(14) + t114) * t86; t120; 0; MDP(17); t68 * MDP(23) + t69 * MDP(24); t66 * MDP(23) - t67 * MDP(24) + t110 + t130; (t89 * MDP(23) + MDP(20)) * t104 + (-t89 * MDP(24) - MDP(21)) * t101; t113; t114; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
