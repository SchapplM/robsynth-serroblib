% Calculate Gravitation load on the joints for
% S6RRRRRR6
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
%   see S6RRRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:21:36
% EndTime: 2019-03-10 04:21:37
% DurationCPUTime: 0.59s
% Computational Cost: add. (536->112), mult. (852->195), div. (0->0), fcn. (1032->14), ass. (0->55)
t105 = qJ(5) + qJ(6);
t101 = sin(t105);
t103 = cos(t105);
t106 = qJ(3) + qJ(4);
t102 = sin(t106);
t104 = cos(t106);
t107 = sin(pkin(6));
t142 = cos(qJ(1));
t123 = t107 * t142;
t110 = sin(qJ(2));
t111 = sin(qJ(1));
t114 = cos(qJ(2));
t135 = cos(pkin(6));
t120 = t135 * t142;
t95 = t110 * t120 + t111 * t114;
t86 = -t102 * t123 + t95 * t104;
t94 = t111 * t110 - t114 * t120;
t146 = t86 * t101 - t94 * t103;
t145 = t94 * t101 + t86 * t103;
t108 = sin(qJ(5));
t112 = cos(qJ(5));
t144 = t86 * t108 - t94 * t112;
t143 = t94 * t108 + t86 * t112;
t141 = g(3) * t107;
t126 = t107 * t114;
t128 = t107 * t111;
t121 = t111 * t135;
t97 = -t110 * t121 + t142 * t114;
t89 = t102 * t128 + t97 * t104;
t96 = t142 * t110 + t114 * t121;
t81 = -t89 * t101 + t96 * t103;
t82 = t96 * t101 + t89 * t103;
t129 = t107 * t110;
t93 = t135 * t102 + t104 * t129;
t140 = (-g(1) * t81 + g(2) * t146 - g(3) * (-t93 * t101 - t103 * t126)) * MDP(37) + (g(1) * t82 + g(2) * t145 - g(3) * (t101 * t126 - t93 * t103)) * MDP(38);
t134 = t101 * t104;
t133 = t103 * t104;
t132 = t104 * t108;
t131 = t104 * t112;
t130 = t104 * t114;
t113 = cos(qJ(3));
t127 = t107 * t113;
t125 = t108 * t114;
t124 = t112 * t114;
t109 = sin(qJ(3));
t122 = -t109 * t123 + t95 * t113;
t117 = t95 * t102 + t104 * t123;
t88 = -t97 * t102 + t104 * t128;
t119 = (g(1) * t89 + g(2) * t86 + g(3) * t93) * MDP(24) + (-t112 * MDP(30) + t108 * MDP(31) - t103 * MDP(37) + t101 * MDP(38) - MDP(23)) * (g(1) * t88 - g(2) * t117 + g(3) * (-t102 * t129 + t135 * t104));
t116 = t95 * t109 + t113 * t123;
t91 = t109 * t128 + t97 * t113;
t90 = -t97 * t109 + t111 * t127;
t84 = t96 * t108 + t89 * t112;
t83 = -t89 * t108 + t96 * t112;
t1 = [(g(1) * t111 - g(2) * t142) * MDP(2) + (g(1) * t142 + g(2) * t111) * MDP(3) + (g(1) * t95 - g(2) * t97) * MDP(9) + (-g(1) * t94 + g(2) * t96) * MDP(10) + (g(1) * t122 - g(2) * t91) * MDP(16) + (-g(1) * t116 - g(2) * t90) * MDP(17) + (g(1) * t86 - g(2) * t89) * MDP(23) + (-g(1) * t117 - g(2) * t88) * MDP(24) + (g(1) * t143 - g(2) * t84) * MDP(30) + (-g(1) * t144 - g(2) * t83) * MDP(31) + (g(1) * t145 - g(2) * t82) * MDP(37) + (-g(1) * t146 - g(2) * t81) * MDP(38); (g(1) * t97 + g(2) * t95 + g(3) * t129) * MDP(10) + (-g(1) * (t97 * t108 - t96 * t131) - g(2) * (t95 * t108 - t94 * t131) - (t104 * t124 + t108 * t110) * t141) * MDP(30) + (-g(1) * (t97 * t112 + t96 * t132) - g(2) * (t95 * t112 + t94 * t132) - (-t104 * t125 + t110 * t112) * t141) * MDP(31) + (-g(1) * (t97 * t101 - t96 * t133) - g(2) * (t95 * t101 - t94 * t133) - (t101 * t110 + t103 * t130) * t141) * MDP(37) + (-g(1) * (t97 * t103 + t96 * t134) - g(2) * (t95 * t103 + t94 * t134) - (-t101 * t130 + t103 * t110) * t141) * MDP(38) + (-t113 * MDP(16) + t109 * MDP(17) - t104 * MDP(23) + t102 * MDP(24) - MDP(9)) * (-g(1) * t96 - g(2) * t94 + g(3) * t126); (-g(1) * t90 + g(2) * t116 - g(3) * (-t109 * t129 + t135 * t113)) * MDP(16) + (g(1) * t91 + g(2) * t122 - g(3) * (-t135 * t109 - t110 * t127)) * MDP(17) + t119; t119; (-g(1) * t83 + g(2) * t144 - g(3) * (-t107 * t124 - t93 * t108)) * MDP(30) + (g(1) * t84 + g(2) * t143 - g(3) * (t107 * t125 - t93 * t112)) * MDP(31) + t140; t140;];
taug  = t1;
