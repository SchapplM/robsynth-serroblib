% Calculate Gravitation load on the joints for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:41:07
% EndTime: 2019-03-10 04:41:09
% DurationCPUTime: 0.67s
% Computational Cost: add. (552->124), mult. (970->215), div. (0->0), fcn. (1198->14), ass. (0->54)
t103 = cos(qJ(4));
t100 = sin(qJ(3));
t104 = cos(qJ(3));
t131 = cos(qJ(1));
t98 = sin(pkin(6));
t112 = t98 * t131;
t101 = sin(qJ(2));
t102 = sin(qJ(1));
t105 = cos(qJ(2));
t115 = cos(pkin(6));
t109 = t115 * t131;
t87 = t101 * t109 + t102 * t105;
t80 = -t100 * t112 + t87 * t104;
t86 = t102 * t101 - t105 * t109;
t99 = sin(qJ(4));
t138 = -t86 * t103 + t80 * t99;
t97 = qJ(4) + qJ(5);
t94 = sin(t97);
t95 = cos(t97);
t137 = t80 * t94 - t86 * t95;
t136 = t80 * t95 + t86 * t94;
t96 = qJ(6) + t97;
t92 = sin(t96);
t93 = cos(t96);
t135 = t80 * t92 - t86 * t93;
t134 = t80 * t93 + t86 * t92;
t133 = t80 * t103 + t86 * t99;
t132 = g(3) * t98;
t117 = t105 * t98;
t111 = t102 * t115;
t89 = -t101 * t111 + t131 * t105;
t83 = t102 * t98 * t100 + t89 * t104;
t88 = t131 * t101 + t105 * t111;
t73 = -t83 * t92 + t88 * t93;
t74 = t83 * t93 + t88 * t92;
t119 = t104 * t98;
t85 = t115 * t100 + t101 * t119;
t125 = (-g(1) * t73 + g(2) * t135 - g(3) * (-t93 * t117 - t85 * t92)) * MDP(37) + (g(1) * t74 + g(2) * t134 - g(3) * (t92 * t117 - t85 * t93)) * MDP(38);
t124 = t101 * t98;
t123 = t104 * t92;
t122 = t104 * t93;
t121 = t104 * t94;
t120 = t104 * t95;
t118 = t104 * t99;
t114 = t103 * t104;
t113 = t104 * t105;
t75 = -t83 * t94 + t88 * t95;
t76 = t83 * t95 + t88 * t94;
t110 = (-g(1) * t75 + g(2) * t137 - g(3) * (-t95 * t117 - t85 * t94)) * MDP(30) + (g(1) * t76 + g(2) * t136 - g(3) * (t94 * t117 - t85 * t95)) * MDP(31) + t125;
t107 = t87 * t100 + t104 * t112;
t82 = -t89 * t100 + t102 * t119;
t78 = t83 * t103 + t88 * t99;
t77 = t88 * t103 - t83 * t99;
t1 = [(g(1) * t102 - g(2) * t131) * MDP(2) + (g(1) * t131 + g(2) * t102) * MDP(3) + (g(1) * t87 - g(2) * t89) * MDP(9) + (-g(1) * t86 + g(2) * t88) * MDP(10) + (g(1) * t80 - g(2) * t83) * MDP(16) + (-g(1) * t107 - g(2) * t82) * MDP(17) + (g(1) * t133 - g(2) * t78) * MDP(23) + (-g(1) * t138 - g(2) * t77) * MDP(24) + (g(1) * t136 - g(2) * t76) * MDP(30) + (-g(1) * t137 - g(2) * t75) * MDP(31) + (g(1) * t134 - g(2) * t74) * MDP(37) + (-g(1) * t135 - g(2) * t73) * MDP(38); (g(1) * t89 + g(2) * t87 + g(3) * t124) * MDP(10) + (-g(1) * (-t88 * t114 + t89 * t99) - g(2) * (-t86 * t114 + t87 * t99) - (t101 * t99 + t103 * t113) * t132) * MDP(23) + (-g(1) * (t89 * t103 + t88 * t118) - g(2) * (t87 * t103 + t86 * t118) - (t101 * t103 - t99 * t113) * t132) * MDP(24) + (-g(1) * (-t88 * t120 + t89 * t94) - g(2) * (-t86 * t120 + t87 * t94) - (t101 * t94 + t95 * t113) * t132) * MDP(30) + (-g(1) * (t88 * t121 + t89 * t95) - g(2) * (t86 * t121 + t87 * t95) - (t101 * t95 - t94 * t113) * t132) * MDP(31) + (-g(1) * (-t88 * t122 + t89 * t92) - g(2) * (-t86 * t122 + t87 * t92) - (t101 * t92 + t93 * t113) * t132) * MDP(37) + (-g(1) * (t88 * t123 + t89 * t93) - g(2) * (t86 * t123 + t87 * t93) - (t101 * t93 - t92 * t113) * t132) * MDP(38) + (-t104 * MDP(16) + t100 * MDP(17) - MDP(9)) * (-g(1) * t88 - g(2) * t86 + g(3) * t117); (g(1) * t83 + g(2) * t80 + g(3) * t85) * MDP(17) + (-MDP(23) * t103 + MDP(24) * t99 - MDP(30) * t95 + MDP(31) * t94 - MDP(37) * t93 + MDP(38) * t92 - MDP(16)) * (g(1) * t82 - g(2) * t107 + g(3) * (-t100 * t124 + t115 * t104)); (-g(1) * t77 + g(2) * t138 - g(3) * (-t103 * t117 - t85 * t99)) * MDP(23) + (g(1) * t78 + g(2) * t133 - g(3) * (-t85 * t103 + t99 * t117)) * MDP(24) + t110; t110; t125;];
taug  = t1;
