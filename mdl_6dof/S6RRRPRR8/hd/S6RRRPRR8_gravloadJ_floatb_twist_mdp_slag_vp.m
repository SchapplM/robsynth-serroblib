% Calculate Gravitation load on the joints for
% S6RRRPRR8
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
%   see S6RRRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:52:23
% EndTime: 2019-03-09 18:52:26
% DurationCPUTime: 0.79s
% Computational Cost: add. (438->120), mult. (746->208), div. (0->0), fcn. (888->14), ass. (0->57)
t156 = MDP(10) - MDP(18);
t106 = qJ(5) + qJ(6);
t103 = sin(t106);
t104 = cos(t106);
t105 = qJ(3) + pkin(12);
t101 = sin(t105);
t102 = cos(t105);
t107 = sin(pkin(6));
t116 = cos(qJ(1));
t130 = t107 * t116;
t111 = sin(qJ(2));
t112 = sin(qJ(1));
t115 = cos(qJ(2));
t140 = cos(pkin(6));
t124 = t116 * t140;
t94 = t111 * t124 + t112 * t115;
t86 = t101 * t130 - t94 * t102;
t93 = t111 * t112 - t115 * t124;
t155 = t103 * t86 + t104 * t93;
t154 = -t103 * t93 + t104 * t86;
t109 = sin(qJ(5));
t113 = cos(qJ(5));
t153 = t109 * t86 + t113 * t93;
t152 = -t109 * t93 + t113 * t86;
t151 = g(1) * t116 + g(2) * t112;
t146 = g(3) * t107;
t131 = t107 * t115;
t133 = t107 * t112;
t125 = t112 * t140;
t96 = -t111 * t125 + t115 * t116;
t88 = t101 * t133 + t102 * t96;
t95 = t116 * t111 + t115 * t125;
t80 = -t103 * t88 + t104 * t95;
t81 = t103 * t95 + t104 * t88;
t134 = t107 * t111;
t92 = t140 * t101 + t102 * t134;
t145 = (-g(1) * t80 - g(2) * t155 - g(3) * (-t103 * t92 - t104 * t131)) * MDP(32) + (g(1) * t81 - g(2) * t154 - g(3) * (t103 * t131 - t104 * t92)) * MDP(33);
t139 = t102 * t103;
t138 = t102 * t104;
t137 = t102 * t109;
t136 = t102 * t113;
t135 = t102 * t115;
t114 = cos(qJ(3));
t132 = t107 * t114;
t129 = t109 * t115;
t128 = t113 * t115;
t110 = sin(qJ(3));
t126 = -t110 * t130 + t114 * t94;
t89 = -t110 * t96 + t112 * t132;
t122 = t94 * t110 + t114 * t130;
t118 = -g(1) * t95 - g(2) * t93 + g(3) * t131;
t108 = -qJ(4) - pkin(9);
t100 = pkin(3) * t114 + pkin(2);
t90 = t110 * t133 + t114 * t96;
t83 = t109 * t95 + t113 * t88;
t82 = -t109 * t88 + t113 * t95;
t1 = [(g(1) * t112 - g(2) * t116) * MDP(2) + t151 * MDP(3) + (g(1) * t94 - g(2) * t96) * MDP(9) + (g(1) * t126 - g(2) * t90) * MDP(16) + (-g(1) * t122 - g(2) * t89) * MDP(17) + (-g(1) * (-t112 * pkin(1) - t94 * t100 + t93 * t108) - g(2) * (pkin(1) * t116 + t96 * t100 - t95 * t108) - t151 * t107 * (pkin(3) * t110 + pkin(8))) * MDP(19) + (-g(1) * t152 - g(2) * t83) * MDP(25) + (g(1) * t153 - g(2) * t82) * MDP(26) + (-g(1) * t154 - g(2) * t81) * MDP(32) + (g(1) * t155 - g(2) * t80) * MDP(33) - t156 * (g(1) * t93 - g(2) * t95); (-g(1) * (-t100 * t95 - t108 * t96) - g(2) * (-t100 * t93 - t108 * t94) - (t100 * t115 - t108 * t111) * t146) * MDP(19) + (-g(1) * (t109 * t96 - t95 * t136) - g(2) * (t109 * t94 - t93 * t136) - (t102 * t128 + t109 * t111) * t146) * MDP(25) + (-g(1) * (t113 * t96 + t95 * t137) - g(2) * (t113 * t94 + t93 * t137) - (-t102 * t129 + t111 * t113) * t146) * MDP(26) + (-g(1) * (t103 * t96 - t95 * t138) - g(2) * (t103 * t94 - t93 * t138) - (t103 * t111 + t104 * t135) * t146) * MDP(32) + (-g(1) * (t104 * t96 + t95 * t139) - g(2) * (t104 * t94 + t93 * t139) - (-t103 * t135 + t104 * t111) * t146) * MDP(33) + t156 * (g(1) * t96 + g(2) * t94 + g(3) * t134) + (-t114 * MDP(16) + t110 * MDP(17) - MDP(9)) * t118; (g(1) * t90 + g(2) * t126 - g(3) * (-t140 * t110 - t111 * t132)) * MDP(17) + (-MDP(25) * t113 + MDP(26) * t109 - MDP(32) * t104 + MDP(33) * t103) * (g(1) * (-t101 * t96 + t102 * t133) + g(2) * (-t94 * t101 - t102 * t130) + g(3) * (-t101 * t134 + t140 * t102)) + (pkin(3) * MDP(19) + MDP(16)) * (-g(3) * (-t110 * t134 + t140 * t114) + g(2) * t122 - g(1) * t89); t118 * MDP(19); (-g(1) * t82 - g(2) * t153 - g(3) * (-t107 * t128 - t109 * t92)) * MDP(25) + (g(1) * t83 - g(2) * t152 - g(3) * (t107 * t129 - t113 * t92)) * MDP(26) + t145; t145;];
taug  = t1;
