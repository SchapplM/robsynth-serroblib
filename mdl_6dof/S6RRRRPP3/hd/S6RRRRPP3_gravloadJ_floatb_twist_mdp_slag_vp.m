% Calculate Gravitation load on the joints for
% S6RRRRPP3
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
%   see S6RRRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:56:30
% EndTime: 2019-03-09 20:56:32
% DurationCPUTime: 0.62s
% Computational Cost: add. (484->102), mult. (649->128), div. (0->0), fcn. (624->8), ass. (0->58)
t185 = MDP(24) - MDP(27) - MDP(30);
t183 = MDP(23) - MDP(26) + MDP(31);
t182 = MDP(17) - MDP(25) - MDP(29);
t134 = qJ(2) + qJ(3);
t131 = sin(t134);
t137 = sin(qJ(1));
t140 = cos(qJ(1));
t155 = g(1) * t140 + g(2) * t137;
t178 = t131 * t155;
t132 = cos(t134);
t181 = -t132 * pkin(3) - t131 * pkin(9);
t136 = sin(qJ(2));
t177 = pkin(2) * t136;
t141 = -pkin(8) - pkin(7);
t174 = g(2) * t141;
t126 = t131 * pkin(5);
t173 = -pkin(4) - qJ(6);
t135 = sin(qJ(4));
t172 = qJ(5) * t135;
t171 = t131 * t135;
t138 = cos(qJ(4));
t170 = t131 * t138;
t169 = t131 * t140;
t168 = t132 * t137;
t167 = t132 * t138;
t166 = t132 * t140;
t165 = t135 * t137;
t164 = t137 * t138;
t163 = t138 * t140;
t162 = t140 * t135;
t161 = -pkin(3) - t172;
t110 = t132 * t165 + t163;
t111 = t132 * t164 - t162;
t160 = -t110 * pkin(4) + qJ(5) * t111;
t112 = t132 * t162 - t164;
t113 = t132 * t163 + t165;
t159 = -t112 * pkin(4) + qJ(5) * t113;
t158 = pkin(4) * t167 + t132 * t172 - t181;
t118 = pkin(9) * t168;
t157 = -t137 * t177 + t118;
t122 = pkin(9) * t166;
t156 = -t140 * t177 + t122;
t139 = cos(qJ(2));
t133 = t139 * pkin(2);
t130 = t133 + pkin(1);
t152 = -t130 + t181;
t151 = qJ(6) * t167 + t126 + t158;
t150 = -t111 * pkin(4) - t110 * qJ(5) - t140 * t141;
t149 = pkin(3) * t166 + t113 * pkin(4) + pkin(9) * t169 + t112 * qJ(5) + t140 * t130;
t147 = g(1) * t112 + g(2) * t110 + g(3) * t171;
t146 = g(1) * t113 + g(2) * t111 + g(3) * t170;
t144 = t182 * (g(3) * t131 + t132 * t155) + (-t135 * t185 + t183 * t138 + MDP(16)) * (-g(3) * t132 + t178);
t143 = (pkin(4) * t138 - t161) * t178;
t142 = (-t138 * t173 - t161) * t178;
t123 = pkin(5) * t166;
t119 = pkin(5) * t168;
t114 = qJ(5) * t170;
t1 = [t155 * MDP(3) + (-g(1) * t150 - g(2) * t149 + (-g(1) * t152 + t174) * t137) * MDP(28) + (-g(1) * (-t111 * qJ(6) + t150) - g(2) * (pkin(5) * t169 + t113 * qJ(6) + t149) + (-g(1) * (t152 - t126) + t174) * t137) * MDP(32) + t183 * (g(1) * t111 - g(2) * t113) - t185 * (g(1) * t110 - g(2) * t112) + (-t136 * MDP(10) + t132 * MDP(16) + t139 * MDP(9) - t182 * t131 + MDP(2)) * (g(1) * t137 - g(2) * t140); t144 + (-g(3) * t139 + t155 * t136) * MDP(9) + (g(3) * t136 + t155 * t139) * MDP(10) + (-g(1) * (t123 + t156) - g(2) * (t119 + t157) - g(3) * (t133 + t151) + t142) * MDP(32) + (-g(1) * t156 - g(2) * t157 - g(3) * (t133 + t158) + t143) * MDP(28); (-g(1) * t122 - g(2) * t118 - g(3) * t158 + t143) * MDP(28) + (-g(1) * (t122 + t123) - g(2) * (t118 + t119) - g(3) * t151 + t142) * MDP(32) + t144; (-g(1) * t159 - g(2) * t160 - g(3) * (-pkin(4) * t171 + t114)) * MDP(28) + (-g(1) * (-qJ(6) * t112 + t159) - g(2) * (-qJ(6) * t110 + t160) - g(3) * (t171 * t173 + t114)) * MDP(32) + t183 * t147 + t185 * t146; -(MDP(28) + MDP(32)) * t147; -t146 * MDP(32);];
taug  = t1;
