% Calculate Gravitation load on the joints for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:44:15
% EndTime: 2019-03-09 15:44:18
% DurationCPUTime: 1.04s
% Computational Cost: add. (507->150), mult. (871->242), div. (0->0), fcn. (1011->14), ass. (0->70)
t143 = sin(qJ(2));
t144 = sin(qJ(1));
t146 = cos(qJ(2));
t147 = cos(qJ(1));
t180 = cos(pkin(6));
t158 = t147 * t180;
t113 = t143 * t144 - t146 * t158;
t136 = pkin(12) + qJ(6);
t131 = sin(t136);
t133 = cos(t136);
t114 = t143 * t158 + t144 * t146;
t137 = qJ(3) + pkin(11);
t132 = sin(t137);
t134 = cos(t137);
t139 = sin(pkin(6));
t167 = t139 * t147;
t96 = t114 * t134 - t132 * t167;
t184 = -t113 * t133 + t131 * t96;
t183 = t113 * t131 + t133 * t96;
t182 = MDP(10) - MDP(18);
t181 = g(3) * t139;
t159 = t144 * t180;
t116 = -t143 * t159 + t146 * t147;
t142 = sin(qJ(3));
t177 = t116 * t142;
t176 = t131 * t134;
t175 = t133 * t134;
t138 = sin(pkin(12));
t174 = t134 * t138;
t140 = cos(pkin(12));
t173 = t134 * t140;
t172 = t134 * t146;
t171 = t139 * t143;
t170 = t139 * t144;
t145 = cos(qJ(3));
t169 = t139 * t145;
t168 = t139 * t146;
t141 = -qJ(4) - pkin(9);
t166 = t141 * t143;
t130 = pkin(3) * t145 + pkin(2);
t165 = -t113 * t130 - t114 * t141;
t115 = t147 * t143 + t146 * t159;
t164 = -t115 * t130 - t116 * t141;
t163 = t142 * t170;
t162 = t142 * t171;
t161 = t144 * t169;
t160 = t145 * t167;
t124 = t142 * t167;
t157 = t180 * t145;
t156 = t114 * t145 - t124;
t154 = pkin(4) * t134 + qJ(5) * t132;
t153 = t147 * pkin(1) + pkin(3) * t163 + pkin(8) * t170 - t115 * t141 + t116 * t130;
t95 = t114 * t132 + t134 * t167;
t152 = t114 * t142 + t160;
t151 = -pkin(1) * t144 + pkin(3) * t124 + pkin(8) * t167 + t113 * t141 - t114 * t130;
t107 = t132 * t171 - t134 * t180;
t99 = t116 * t132 - t134 * t170;
t150 = g(1) * t99 + g(2) * t95 + g(3) * t107;
t149 = -g(1) * t115 - g(2) * t113 + g(3) * t168;
t148 = g(1) * t116 + g(2) * t114 + g(3) * t171;
t129 = pkin(3) * t157;
t121 = pkin(3) * t161;
t117 = t130 * t168;
t108 = t132 * t180 + t134 * t171;
t102 = t116 * t145 + t163;
t101 = t161 - t177;
t100 = t116 * t134 + t132 * t170;
t93 = t100 * t133 + t115 * t131;
t92 = -t100 * t131 + t115 * t133;
t1 = [(g(1) * t144 - g(2) * t147) * MDP(2) + (g(1) * t147 + g(2) * t144) * MDP(3) + (g(1) * t114 - g(2) * t116) * MDP(9) + (g(1) * t156 - g(2) * t102) * MDP(16) + (-g(1) * t152 - g(2) * t101) * MDP(17) + (-g(1) * t151 - g(2) * t153) * MDP(19) + (-g(1) * (-t113 * t138 - t140 * t96) - g(2) * (t100 * t140 + t115 * t138)) * MDP(20) + (-g(1) * (-t113 * t140 + t138 * t96) - g(2) * (-t100 * t138 + t115 * t140)) * MDP(21) + (g(1) * t95 - g(2) * t99) * MDP(22) + (-g(1) * (-pkin(4) * t96 - qJ(5) * t95 + t151) - g(2) * (pkin(4) * t100 + qJ(5) * t99 + t153)) * MDP(23) + (g(1) * t183 - g(2) * t93) * MDP(29) + (-g(1) * t184 - g(2) * t92) * MDP(30) - t182 * (g(1) * t113 - g(2) * t115); (-g(1) * t164 - g(2) * t165 - g(3) * (-t139 * t166 + t117)) * MDP(19) + (-g(1) * (-t115 * t173 + t116 * t138) - g(2) * (-t113 * t173 + t114 * t138) - (t138 * t143 + t140 * t172) * t181) * MDP(20) + (-g(1) * (t115 * t174 + t116 * t140) - g(2) * (t113 * t174 + t114 * t140) - (-t138 * t172 + t140 * t143) * t181) * MDP(21) + (-g(1) * (-t115 * t154 + t164) - g(2) * (-t113 * t154 + t165) - g(3) * t117 - (t146 * t154 - t166) * t181) * MDP(23) + (-g(1) * (-t115 * t175 + t116 * t131) - g(2) * (-t113 * t175 + t114 * t131) - (t131 * t143 + t133 * t172) * t181) * MDP(29) + (-g(1) * (t115 * t176 + t116 * t133) - g(2) * (t113 * t176 + t114 * t133) - (-t131 * t172 + t133 * t143) * t181) * MDP(30) + t182 * t148 + (-t145 * MDP(16) + t142 * MDP(17) - t132 * MDP(22) - MDP(9)) * t149; (-g(1) * t101 + g(2) * t152 - g(3) * (t157 - t162)) * MDP(16) + (g(1) * t102 + g(2) * t156 - g(3) * (-t142 * t180 - t143 * t169)) * MDP(17) + (-g(1) * t121 - g(3) * t129 + (g(2) * t160 + t142 * t148) * pkin(3)) * MDP(19) + (-g(1) * t100 - g(2) * t96 - g(3) * t108) * MDP(22) + (-g(1) * (-pkin(3) * t177 - pkin(4) * t99 + qJ(5) * t100 + t121) - g(2) * (-pkin(3) * t152 - t95 * pkin(4) + t96 * qJ(5)) - g(3) * (-pkin(3) * t162 - pkin(4) * t107 + qJ(5) * t108 + t129)) * MDP(23) + (MDP(20) * t140 - MDP(21) * t138 + MDP(29) * t133 - MDP(30) * t131) * t150; (MDP(19) + MDP(23)) * t149; -t150 * MDP(23); (-g(1) * t92 + g(2) * t184 - g(3) * (-t108 * t131 - t133 * t168)) * MDP(29) + (g(1) * t93 + g(2) * t183 - g(3) * (-t108 * t133 + t131 * t168)) * MDP(30);];
taug  = t1;
