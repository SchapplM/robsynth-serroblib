% Calculate Gravitation load on the joints for
% S6RRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:39:27
% EndTime: 2019-03-09 18:39:29
% DurationCPUTime: 0.58s
% Computational Cost: add. (438->99), mult. (648->168), div. (0->0), fcn. (752->12), ass. (0->48)
t138 = MDP(10) - MDP(18);
t102 = cos(qJ(6));
t105 = cos(qJ(1));
t96 = sin(pkin(6));
t121 = t105 * t96;
t100 = sin(qJ(2));
t101 = sin(qJ(1));
t104 = cos(qJ(2));
t120 = cos(pkin(6));
t115 = t105 * t120;
t85 = t100 * t115 + t101 * t104;
t95 = qJ(3) + pkin(12) + qJ(5);
t92 = sin(t95);
t93 = cos(t95);
t76 = -t92 * t121 + t85 * t93;
t84 = t100 * t101 - t104 * t115;
t98 = sin(qJ(6));
t137 = -t102 * t84 + t76 * t98;
t136 = t102 * t76 + t84 * t98;
t135 = g(1) * t105 + g(2) * t101;
t132 = g(3) * t96;
t128 = t93 * t98;
t127 = t100 * t96;
t126 = t101 * t96;
t124 = t102 * t93;
t103 = cos(qJ(3));
t123 = t103 * t96;
t122 = t104 * t98;
t119 = t102 * t104;
t99 = sin(qJ(3));
t117 = t103 * t85 - t99 * t121;
t116 = t101 * t120;
t112 = t93 * t121 + t85 * t92;
t87 = -t100 * t116 + t104 * t105;
t78 = t93 * t126 - t87 * t92;
t79 = t92 * t126 + t87 * t93;
t83 = t120 * t92 + t93 * t127;
t114 = (g(1) * t79 + g(2) * t76 + g(3) * t83) * MDP(26) + (-MDP(32) * t102 + MDP(33) * t98 - MDP(25)) * (g(1) * t78 - g(2) * t112 + g(3) * (t120 * t93 - t92 * t127));
t80 = t101 * t123 - t87 * t99;
t111 = t103 * t121 + t85 * t99;
t86 = t105 * t100 + t104 * t116;
t107 = -g(1) * t86 - g(2) * t84 + t104 * t132;
t97 = -qJ(4) - pkin(9);
t94 = pkin(3) * t103 + pkin(2);
t81 = t103 * t87 + t99 * t126;
t74 = t102 * t79 + t86 * t98;
t73 = t102 * t86 - t79 * t98;
t1 = [(g(1) * t101 - g(2) * t105) * MDP(2) + t135 * MDP(3) + (g(1) * t85 - g(2) * t87) * MDP(9) + (g(1) * t117 - g(2) * t81) * MDP(16) + (-g(1) * t111 - g(2) * t80) * MDP(17) + (-g(1) * (-t101 * pkin(1) + t84 * t97 - t85 * t94) - g(2) * (pkin(1) * t105 - t86 * t97 + t87 * t94) - t135 * t96 * (pkin(3) * t99 + pkin(8))) * MDP(19) + (g(1) * t76 - g(2) * t79) * MDP(25) + (-g(1) * t112 - g(2) * t78) * MDP(26) + (g(1) * t136 - g(2) * t74) * MDP(32) + (-g(1) * t137 - g(2) * t73) * MDP(33) - t138 * (g(1) * t84 - g(2) * t86); (-g(1) * (-t86 * t94 - t87 * t97) - g(2) * (-t84 * t94 - t85 * t97) - (-t100 * t97 + t104 * t94) * t132) * MDP(19) + (-g(1) * (-t86 * t124 + t87 * t98) - g(2) * (-t84 * t124 + t85 * t98) - (t100 * t98 + t93 * t119) * t132) * MDP(32) + (-g(1) * (t102 * t87 + t86 * t128) - g(2) * (t102 * t85 + t84 * t128) - (t100 * t102 - t93 * t122) * t132) * MDP(33) + t138 * (g(1) * t87 + g(2) * t85 + g(3) * t127) + (-MDP(16) * t103 + t99 * MDP(17) - t93 * MDP(25) + MDP(26) * t92 - MDP(9)) * t107; (g(1) * t81 + g(2) * t117 - g(3) * (-t100 * t123 - t120 * t99)) * MDP(17) + t114 + (pkin(3) * MDP(19) + MDP(16)) * (-g(3) * (t120 * t103 - t99 * t127) + g(2) * t111 - g(1) * t80); t107 * MDP(19); t114; (-g(1) * t73 + g(2) * t137 - g(3) * (-t96 * t119 - t83 * t98)) * MDP(32) + (g(1) * t74 + g(2) * t136 - g(3) * (-t102 * t83 + t96 * t122)) * MDP(33);];
taug  = t1;
