% Calculate Gravitation load on the joints for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:06:46
% EndTime: 2019-03-10 05:06:49
% DurationCPUTime: 1.10s
% Computational Cost: add. (802->152), mult. (1914->270), div. (0->0), fcn. (2454->16), ass. (0->72)
t148 = sin(qJ(2));
t151 = cos(qJ(2));
t152 = cos(qJ(1));
t177 = cos(pkin(6));
t162 = t152 * t177;
t178 = sin(qJ(1));
t132 = t178 * t148 - t151 * t162;
t133 = t148 * t162 + t178 * t151;
t147 = sin(qJ(3));
t176 = cos(pkin(7));
t163 = t147 * t176;
t143 = sin(pkin(7));
t144 = sin(pkin(6));
t171 = t143 * t144;
t165 = t152 * t171;
t179 = cos(qJ(3));
t112 = -t132 * t163 + t133 * t179 - t147 * t165;
t164 = t144 * t176;
t125 = t132 * t143 - t152 * t164;
t142 = qJ(4) + qJ(5);
t140 = sin(t142);
t141 = cos(t142);
t101 = t112 * t141 + t125 * t140;
t157 = t176 * t179;
t111 = t132 * t157 + t133 * t147 + t179 * t165;
t145 = sin(qJ(6));
t149 = cos(qJ(6));
t187 = t101 * t145 - t111 * t149;
t186 = t101 * t149 + t111 * t145;
t183 = t112 * t140 - t125 * t141;
t146 = sin(qJ(4));
t150 = cos(qJ(4));
t182 = t112 * t150 + t125 * t146;
t181 = t112 * t146 - t125 * t150;
t158 = t177 * t178;
t154 = t152 * t148 + t151 * t158;
t180 = t154 * t176 - t178 * t171;
t175 = t140 * t143;
t174 = t141 * t143;
t173 = t141 * t145;
t172 = t141 * t149;
t170 = t143 * t146;
t169 = t143 * t150;
t168 = t144 * t148;
t167 = t144 * t151;
t166 = t143 * t168;
t161 = t177 * t143;
t134 = -t148 * t158 + t152 * t151;
t116 = t134 * t179 - t180 * t147;
t126 = t154 * t143 + t178 * t164;
t103 = -t116 * t140 + t126 * t141;
t104 = t116 * t141 + t126 * t140;
t123 = t147 * t161 + (t179 * t148 + t151 * t163) * t144;
t131 = -t143 * t167 + t177 * t176;
t110 = t123 * t141 + t131 * t140;
t160 = (g(1) * t104 + g(2) * t101 + g(3) * t110) * MDP(31) + (-t149 * MDP(37) + t145 * MDP(38) - MDP(30)) * (g(1) * t103 - g(2) * t183 + g(3) * (-t123 * t140 + t131 * t141));
t130 = (-t148 * t163 + t179 * t151) * t144;
t129 = (t147 * t151 + t148 * t157) * t144;
t122 = t147 * t168 - t157 * t167 - t179 * t161;
t121 = t130 * t141 + t140 * t166;
t120 = -t134 * t163 - t154 * t179;
t119 = t134 * t157 - t154 * t147;
t118 = -t132 * t179 - t133 * t163;
t117 = -t132 * t147 + t133 * t157;
t115 = t134 * t147 + t180 * t179;
t108 = t120 * t141 + t134 * t175;
t107 = t118 * t141 + t133 * t175;
t106 = t116 * t150 + t126 * t146;
t105 = -t116 * t146 + t126 * t150;
t99 = t104 * t149 + t115 * t145;
t98 = -t104 * t145 + t115 * t149;
t1 = [(g(1) * t178 - g(2) * t152) * MDP(2) + (g(1) * t152 + g(2) * t178) * MDP(3) + (g(1) * t133 - g(2) * t134) * MDP(9) + (-g(1) * t132 + g(2) * t154) * MDP(10) + (g(1) * t112 - g(2) * t116) * MDP(16) + (-g(1) * t111 + g(2) * t115) * MDP(17) + (g(1) * t182 - g(2) * t106) * MDP(23) + (-g(1) * t181 - g(2) * t105) * MDP(24) + (g(1) * t101 - g(2) * t104) * MDP(30) + (-g(1) * t183 - g(2) * t103) * MDP(31) + (g(1) * t186 - g(2) * t99) * MDP(37) + (-g(1) * t187 - g(2) * t98) * MDP(38); (g(1) * t154 + g(2) * t132 - g(3) * t167) * MDP(9) + (g(1) * t134 + g(2) * t133 + g(3) * t168) * MDP(10) + (-g(1) * t120 - g(2) * t118 - g(3) * t130) * MDP(16) + (g(1) * t119 + g(2) * t117 + g(3) * t129) * MDP(17) + (-g(1) * (t120 * t150 + t134 * t170) - g(2) * (t118 * t150 + t133 * t170) - g(3) * (t130 * t150 + t146 * t166)) * MDP(23) + (-g(1) * (-t120 * t146 + t134 * t169) - g(2) * (-t118 * t146 + t133 * t169) - g(3) * (-t130 * t146 + t150 * t166)) * MDP(24) + (-g(1) * t108 - g(2) * t107 - g(3) * t121) * MDP(30) + (-g(1) * (-t120 * t140 + t134 * t174) - g(2) * (-t118 * t140 + t133 * t174) - g(3) * (-t130 * t140 + t141 * t166)) * MDP(31) + (-g(1) * (t108 * t149 + t119 * t145) - g(2) * (t107 * t149 + t117 * t145) - g(3) * (t121 * t149 + t129 * t145)) * MDP(37) + (-g(1) * (-t108 * t145 + t119 * t149) - g(2) * (-t107 * t145 + t117 * t149) - g(3) * (-t121 * t145 + t129 * t149)) * MDP(38); (g(1) * t116 + g(2) * t112 + g(3) * t123) * MDP(17) + (-g(1) * (-t115 * t172 + t116 * t145) - g(2) * (-t111 * t172 + t112 * t145) - g(3) * (-t122 * t172 + t123 * t145)) * MDP(37) + (-g(1) * (t115 * t173 + t116 * t149) - g(2) * (t111 * t173 + t112 * t149) - g(3) * (t122 * t173 + t123 * t149)) * MDP(38) + (MDP(23) * t150 - MDP(24) * t146 + t141 * MDP(30) - MDP(31) * t140 + MDP(16)) * (g(1) * t115 + g(2) * t111 + g(3) * t122); (-g(1) * t105 + g(2) * t181 - g(3) * (-t123 * t146 + t131 * t150)) * MDP(23) + (g(1) * t106 + g(2) * t182 - g(3) * (-t123 * t150 - t131 * t146)) * MDP(24) + t160; t160; (-g(1) * t98 + g(2) * t187 - g(3) * (-t110 * t145 + t122 * t149)) * MDP(37) + (g(1) * t99 + g(2) * t186 - g(3) * (-t110 * t149 - t122 * t145)) * MDP(38);];
taug  = t1;
