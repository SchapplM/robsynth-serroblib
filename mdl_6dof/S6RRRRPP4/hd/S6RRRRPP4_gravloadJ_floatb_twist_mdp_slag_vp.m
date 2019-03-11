% Calculate Gravitation load on the joints for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRRPP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:02:36
% EndTime: 2019-03-09 21:02:37
% DurationCPUTime: 0.54s
% Computational Cost: add. (475->115), mult. (527->160), div. (0->0), fcn. (508->10), ass. (0->58)
t135 = qJ(3) + qJ(4);
t128 = cos(t135);
t139 = cos(qJ(3));
t117 = t139 * pkin(3) + pkin(4) * t128;
t115 = pkin(2) + t117;
t140 = cos(qJ(2));
t110 = t140 * t115;
t134 = -qJ(5) - pkin(9) - pkin(8);
t137 = sin(qJ(2));
t161 = t134 * t137;
t171 = t110 - t161;
t170 = MDP(10) - MDP(25) - MDP(28);
t138 = sin(qJ(1));
t141 = cos(qJ(1));
t150 = g(1) * t141 + g(2) * t138;
t101 = -g(3) * t140 + t150 * t137;
t127 = sin(t135);
t169 = pkin(4) * t127;
t126 = pkin(10) + t135;
t123 = sin(t126);
t168 = pkin(5) * t123;
t165 = g(3) * t137;
t124 = cos(t126);
t163 = t124 * t137;
t162 = t128 * t141;
t160 = t138 * t128;
t159 = t138 * t140;
t158 = t140 * t141;
t136 = sin(qJ(3));
t116 = pkin(3) * t136 + t169;
t114 = t141 * t116;
t157 = t141 * t123;
t156 = -t140 * t114 + t138 * t117;
t155 = t127 * t158;
t93 = t123 * t159 + t124 * t141;
t94 = t124 * t159 - t157;
t154 = -t93 * pkin(5) + t94 * qJ(6);
t95 = -t138 * t124 + t140 * t157;
t96 = t138 * t123 + t124 * t158;
t153 = -t95 * pkin(5) + t96 * qJ(6);
t100 = t138 * t127 + t128 * t158;
t145 = g(1) * t95 + g(2) * t93 + t123 * t165;
t97 = t127 * t159 + t162;
t98 = t127 * t141 - t128 * t159;
t99 = -t155 + t160;
t152 = t145 * MDP(27) + (-g(1) * t96 - g(2) * t94 - g(3) * t163) * MDP(29) + (-g(1) * t99 + g(2) * t97 + t127 * t165) * MDP(23) + (g(1) * t100 - g(2) * t98 + t128 * t165) * MDP(24);
t151 = -t116 * t159 - t117 * t141;
t148 = pkin(5) * t124 + qJ(6) * t123;
t144 = t115 * t158 + (pkin(1) - t161) * t141 + (pkin(7) + t116) * t138;
t143 = t141 * pkin(7) + t114 + (-pkin(1) - t171) * t138;
t102 = t150 * t140 + t165;
t122 = pkin(4) * t160;
t118 = qJ(6) * t163;
t109 = t138 * t136 + t139 * t158;
t108 = -t136 * t158 + t138 * t139;
t107 = t136 * t141 - t139 * t159;
t106 = t136 * t159 + t139 * t141;
t1 = [t150 * MDP(3) + (-g(1) * t107 - g(2) * t109) * MDP(16) + (-g(1) * t106 - g(2) * t108) * MDP(17) + (-g(1) * t98 - g(2) * t100) * MDP(23) + (-g(1) * t97 - g(2) * t99) * MDP(24) + (-g(1) * t143 - g(2) * t144) * MDP(26) + (g(1) * t94 - g(2) * t96) * MDP(27) + (g(1) * t93 - g(2) * t95) * MDP(29) + (-g(1) * (-pkin(5) * t94 - qJ(6) * t93 + t143) - g(2) * (t96 * pkin(5) + t95 * qJ(6) + t144)) * MDP(30) + (t140 * MDP(9) - t170 * t137 + MDP(2)) * (g(1) * t138 - g(2) * t141); (-g(3) * t171 + t150 * (t115 * t137 + t134 * t140)) * MDP(26) + (-g(3) * t110 + (-g(3) * t148 + t150 * t134) * t140 + (g(3) * t134 + t150 * (t115 + t148)) * t137) * MDP(30) + t170 * t102 + (t139 * MDP(16) - t136 * MDP(17) + t128 * MDP(23) - t127 * MDP(24) + t124 * MDP(27) + t123 * MDP(29) + MDP(9)) * t101; (-g(1) * t108 + g(2) * t106 + t136 * t165) * MDP(16) + (g(1) * t109 - g(2) * t107 + t139 * t165) * MDP(17) + (-g(1) * t156 - g(2) * t151 + t116 * t165) * MDP(26) + (-g(1) * (t153 + t156) - g(2) * (t151 + t154) - g(3) * (t118 + (-t116 - t168) * t137)) * MDP(30) + t152; (-g(1) * t122 + (g(2) * t162 + t102 * t127) * pkin(4)) * MDP(26) + (-g(1) * (-pkin(4) * t155 + t122 + t153) - g(2) * (-t97 * pkin(4) + t154) - g(3) * (t118 + (-t168 - t169) * t137)) * MDP(30) + t152; (-MDP(26) - MDP(30)) * t101; -t145 * MDP(30);];
taug  = t1;
