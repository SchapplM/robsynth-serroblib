% Calculate Gravitation load on the joints for
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:51:41
% EndTime: 2019-03-09 20:51:42
% DurationCPUTime: 0.54s
% Computational Cost: add. (482->102), mult. (644->132), div. (0->0), fcn. (617->8), ass. (0->59)
t182 = MDP(17) - MDP(26) + MDP(31);
t179 = MDP(23) + MDP(25) + MDP(29);
t178 = MDP(24) - MDP(27) - MDP(30);
t132 = qJ(2) + qJ(3);
t129 = sin(t132);
t135 = sin(qJ(1));
t138 = cos(qJ(1));
t154 = g(1) * t138 + g(2) * t135;
t181 = t154 * t129;
t130 = cos(t132);
t180 = -t130 * pkin(3) - t129 * pkin(9);
t177 = -pkin(4) - pkin(5);
t134 = sin(qJ(2));
t176 = pkin(2) * t134;
t175 = g(1) * t135;
t139 = -pkin(8) - pkin(7);
t172 = g(2) * t139;
t133 = sin(qJ(4));
t171 = qJ(5) * t133;
t170 = t129 * t133;
t136 = cos(qJ(4));
t169 = t129 * t136;
t168 = t129 * t138;
t167 = t130 * t135;
t166 = t130 * t136;
t165 = t130 * t138;
t164 = t133 * t135;
t163 = t135 * t136;
t162 = t136 * t138;
t161 = t138 * t133;
t160 = t135 * t176;
t159 = t138 * t176;
t158 = -pkin(3) - t171;
t110 = t130 * t164 + t162;
t111 = t130 * t163 - t161;
t157 = -t110 * pkin(4) + qJ(5) * t111;
t112 = t130 * t161 - t163;
t113 = t130 * t162 + t164;
t156 = -t112 * pkin(4) + qJ(5) * t113;
t155 = pkin(4) * t166 + t130 * t171 - t180;
t117 = pkin(9) * t167;
t152 = -qJ(6) * t167 + t117;
t121 = pkin(9) * t165;
t151 = -qJ(6) * t165 + t121;
t137 = cos(qJ(2));
t131 = t137 * pkin(2);
t150 = t131 + t155;
t128 = t131 + pkin(1);
t148 = -t128 + t180;
t147 = -t111 * pkin(4) - t110 * qJ(5) - t138 * t139;
t146 = pkin(3) * t165 + t113 * pkin(4) + pkin(9) * t168 + t112 * qJ(5) + t138 * t128;
t144 = g(1) * t112 + g(2) * t110 + g(3) * t170;
t103 = -g(3) * t130 + t181;
t142 = t182 * (g(3) * t129 + t154 * t130) + (-t178 * t133 + t179 * t136 + MDP(16)) * t103;
t141 = (pkin(4) * t136 - t158) * t181;
t140 = (g(3) * qJ(6) + t154 * (-t177 * t136 - t158)) * t129;
t118 = pkin(5) * t166;
t114 = qJ(5) * t169;
t1 = [t154 * MDP(3) + (-g(1) * t147 - g(2) * t146 + (-g(1) * t148 + t172) * t135) * MDP(28) + (-g(1) * (-t111 * pkin(5) + t147) - g(2) * (t113 * pkin(5) - qJ(6) * t168 + t146) + (-g(1) * (qJ(6) * t129 + t148) + t172) * t135) * MDP(32) + t179 * (g(1) * t111 - g(2) * t113) - t178 * (g(1) * t110 - g(2) * t112) - t182 * (-g(2) * t168 + t129 * t175) + (-t134 * MDP(10) + t130 * MDP(16) + t137 * MDP(9) + MDP(2)) * (-g(2) * t138 + t175); t142 + (-g(1) * (t151 - t159) - g(2) * (t152 - t160) - g(3) * (t118 + t150) + t140) * MDP(32) + (-g(1) * (t121 - t159) - g(2) * (t117 - t160) - g(3) * t150 + t141) * MDP(28) + (-g(3) * t137 + t134 * t154) * MDP(9) + (g(3) * t134 + t137 * t154) * MDP(10); (-g(1) * t121 - g(2) * t117 - g(3) * t155 + t141) * MDP(28) + (-g(1) * t151 - g(2) * t152 - g(3) * (t118 + t155) + t140) * MDP(32) + t142; (-g(1) * t156 - g(2) * t157 - g(3) * (-pkin(4) * t170 + t114)) * MDP(28) + (-g(1) * (-pkin(5) * t112 + t156) - g(2) * (-pkin(5) * t110 + t157) - g(3) * (t177 * t170 + t114)) * MDP(32) + t179 * t144 + t178 * (g(1) * t113 + g(2) * t111 + g(3) * t169); -(MDP(28) + MDP(32)) * t144; t103 * MDP(32);];
taug  = t1;
