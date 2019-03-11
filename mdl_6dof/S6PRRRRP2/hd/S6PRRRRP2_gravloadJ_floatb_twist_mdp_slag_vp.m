% Calculate Gravitation load on the joints for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:04:21
% EndTime: 2019-03-09 00:04:22
% DurationCPUTime: 0.60s
% Computational Cost: add. (615->94), mult. (1046->151), div. (0->0), fcn. (1252->12), ass. (0->56)
t147 = sin(pkin(11));
t151 = sin(qJ(2));
t154 = cos(qJ(2));
t178 = cos(pkin(11));
t179 = cos(pkin(6));
t166 = t179 * t178;
t135 = t147 * t154 + t151 * t166;
t168 = t147 * t179;
t137 = -t151 * t168 + t178 * t154;
t146 = qJ(3) + qJ(4);
t144 = sin(t146);
t145 = cos(t146);
t148 = sin(pkin(6));
t167 = t148 * t178;
t173 = t148 * t151;
t174 = t147 * t148;
t192 = g(3) * (-t144 * t173 + t179 * t145) + g(2) * (-t135 * t144 - t145 * t167) + g(1) * (-t137 * t144 + t145 * t174);
t188 = MDP(18) - MDP(27);
t187 = MDP(24) + MDP(26);
t186 = MDP(25) - MDP(28);
t118 = t135 * t145 - t144 * t167;
t120 = t137 * t145 + t144 * t174;
t129 = t179 * t144 + t145 * t173;
t184 = g(1) * t120 + g(2) * t118 + g(3) * t129;
t134 = t147 * t151 - t154 * t166;
t136 = t178 * t151 + t154 * t168;
t183 = g(1) * t136 + g(2) * t134;
t149 = sin(qJ(5));
t176 = t145 * t149;
t152 = cos(qJ(5));
t175 = t145 * t152;
t153 = cos(qJ(3));
t172 = t148 * t153;
t171 = t148 * t154;
t170 = t152 * t154;
t169 = t149 * t171;
t165 = t153 * pkin(3) + pkin(4) * t145 + pkin(10) * t144 + pkin(2);
t164 = t184 * t188 + (t149 * t186 - t152 * t187 - MDP(17)) * t192;
t103 = t118 * t149 - t134 * t152;
t105 = t120 * t149 - t136 * t152;
t121 = t129 * t149 + t148 * t170;
t163 = g(1) * t105 + g(2) * t103 + g(3) * t121;
t157 = -pkin(10) * t184 - t192 * (pkin(5) * t152 + qJ(6) * t149 + pkin(4));
t150 = sin(qJ(3));
t156 = -g(1) * (-t137 * t150 + t147 * t172) - g(2) * (-t135 * t150 - t153 * t167) - g(3) * (-t150 * t173 + t179 * t153);
t155 = -pkin(9) - pkin(8);
t126 = (t145 * t170 + t149 * t151) * t148;
t125 = t145 * t169 - t152 * t173;
t122 = t129 * t152 - t169;
t110 = -t136 * t175 + t137 * t149;
t109 = -t136 * t176 - t137 * t152;
t108 = -t134 * t175 + t135 * t149;
t107 = -t134 * t176 - t135 * t152;
t106 = t120 * t152 + t136 * t149;
t104 = t118 * t152 + t134 * t149;
t1 = [(-MDP(1) - MDP(29)) * g(3); (g(1) * t137 + g(2) * t135 + g(3) * t173) * MDP(4) + (-g(1) * (t110 * pkin(5) + t109 * qJ(6) - t137 * t155) - g(2) * (t108 * pkin(5) + t107 * qJ(6) - t135 * t155) + t183 * t165 + (-t126 * pkin(5) - t125 * qJ(6) - (-t151 * t155 + t165 * t154) * t148) * g(3)) * MDP(29) + t187 * (-g(1) * t110 - g(2) * t108 - g(3) * t126) + t186 * (g(1) * t109 + g(2) * t107 + g(3) * t125) + (-MDP(10) * t153 + MDP(11) * t150 - t145 * MDP(17) + t188 * t144 - MDP(3)) * (g(3) * t171 - t183); t156 * MDP(10) + (-g(1) * (-t137 * t153 - t150 * t174) - g(2) * (-t135 * t153 + t150 * t167) - g(3) * (-t179 * t150 - t151 * t172)) * MDP(11) + (t156 * pkin(3) + t157) * MDP(29) + t164; t157 * MDP(29) + t164; (-g(1) * (-t105 * pkin(5) + t106 * qJ(6)) - g(2) * (-t103 * pkin(5) + t104 * qJ(6)) - g(3) * (-t121 * pkin(5) + t122 * qJ(6))) * MDP(29) + t187 * t163 + t186 * (g(1) * t106 + g(2) * t104 + g(3) * t122); -t163 * MDP(29);];
taug  = t1;
