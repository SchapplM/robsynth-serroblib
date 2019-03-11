% Calculate Gravitation load on the joints for
% S6RRRRRR5
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
%   see S6RRRRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:04:48
% EndTime: 2019-03-10 04:04:50
% DurationCPUTime: 0.50s
% Computational Cost: add. (554->100), mult. (748->171), div. (0->0), fcn. (890->14), ass. (0->51)
t100 = sin(qJ(6));
t104 = cos(qJ(6));
t107 = cos(qJ(1));
t99 = sin(pkin(6));
t122 = t107 * t99;
t102 = sin(qJ(2));
t103 = sin(qJ(1));
t106 = cos(qJ(2));
t121 = cos(pkin(6));
t115 = t107 * t121;
t86 = t102 * t115 + t103 * t106;
t98 = qJ(3) + qJ(4);
t97 = qJ(5) + t98;
t93 = sin(t97);
t94 = cos(t97);
t75 = -t122 * t93 + t86 * t94;
t85 = t102 * t103 - t106 * t115;
t132 = t100 * t75 - t104 * t85;
t131 = t100 * t85 + t104 * t75;
t130 = g(3) * t99;
t128 = t100 * t94;
t127 = t102 * t99;
t126 = t103 * t99;
t124 = t104 * t94;
t105 = cos(qJ(3));
t123 = t105 * t99;
t120 = t100 * t106;
t119 = t104 * t106;
t95 = sin(t98);
t96 = cos(t98);
t118 = -t122 * t95 + t86 * t96;
t101 = sin(qJ(3));
t117 = -t101 * t122 + t105 * t86;
t116 = t103 * t121;
t112 = t122 * t94 + t86 * t93;
t88 = -t102 * t116 + t106 * t107;
t77 = t126 * t94 - t88 * t93;
t78 = t126 * t93 + t88 * t94;
t84 = t121 * t93 + t127 * t94;
t114 = (g(1) * t78 + g(2) * t75 + g(3) * t84) * MDP(31) + (-MDP(37) * t104 + MDP(38) * t100 - MDP(30)) * (g(1) * t77 - g(2) * t112 + g(3) * (t121 * t94 - t127 * t93));
t111 = t122 * t96 + t86 * t95;
t79 = t126 * t96 - t88 * t95;
t80 = t126 * t95 + t88 * t96;
t113 = (-g(1) * t79 + g(2) * t111 - g(3) * (t121 * t96 - t127 * t95)) * MDP(23) + (g(1) * t80 + g(2) * t118 - g(3) * (-t121 * t95 - t127 * t96)) * MDP(24) + t114;
t110 = t101 * t86 + t105 * t122;
t87 = t102 * t107 + t106 * t116;
t82 = t101 * t126 + t105 * t88;
t81 = -t101 * t88 + t103 * t123;
t73 = t100 * t87 + t104 * t78;
t72 = -t100 * t78 + t104 * t87;
t1 = [(g(1) * t103 - g(2) * t107) * MDP(2) + (g(1) * t107 + g(2) * t103) * MDP(3) + (g(1) * t86 - g(2) * t88) * MDP(9) + (-g(1) * t85 + g(2) * t87) * MDP(10) + (g(1) * t117 - g(2) * t82) * MDP(16) + (-g(1) * t110 - g(2) * t81) * MDP(17) + (g(1) * t118 - g(2) * t80) * MDP(23) + (-g(1) * t111 - g(2) * t79) * MDP(24) + (g(1) * t75 - g(2) * t78) * MDP(30) + (-g(1) * t112 - g(2) * t77) * MDP(31) + (g(1) * t131 - g(2) * t73) * MDP(37) + (-g(1) * t132 - g(2) * t72) * MDP(38); (g(1) * t88 + g(2) * t86 + g(3) * t127) * MDP(10) + (-g(1) * (t100 * t88 - t124 * t87) - g(2) * (t100 * t86 - t124 * t85) - (t100 * t102 + t119 * t94) * t130) * MDP(37) + (-g(1) * (t104 * t88 + t128 * t87) - g(2) * (t104 * t86 + t128 * t85) - (t102 * t104 - t120 * t94) * t130) * MDP(38) + (-MDP(16) * t105 + MDP(17) * t101 - MDP(23) * t96 + MDP(24) * t95 - MDP(30) * t94 + MDP(31) * t93 - MDP(9)) * (-g(1) * t87 - g(2) * t85 + t106 * t130); (-g(1) * t81 + g(2) * t110 - g(3) * (-t101 * t127 + t105 * t121)) * MDP(16) + (g(1) * t82 + g(2) * t117 - g(3) * (-t101 * t121 - t102 * t123)) * MDP(17) + t113; t113; t114; (-g(1) * t72 + g(2) * t132 - g(3) * (-t100 * t84 - t119 * t99)) * MDP(37) + (g(1) * t73 + g(2) * t131 - g(3) * (-t104 * t84 + t120 * t99)) * MDP(38);];
taug  = t1;
