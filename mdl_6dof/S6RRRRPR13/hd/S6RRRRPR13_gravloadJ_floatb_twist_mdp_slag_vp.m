% Calculate Gravitation load on the joints for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR13_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR13_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:04:05
% EndTime: 2019-03-10 00:04:09
% DurationCPUTime: 1.21s
% Computational Cost: add. (608->137), mult. (1567->218), div. (0->0), fcn. (1959->12), ass. (0->55)
t150 = sin(qJ(2));
t154 = cos(qJ(2));
t185 = cos(pkin(6));
t194 = cos(qJ(1));
t172 = t185 * t194;
t193 = sin(qJ(1));
t136 = t150 * t172 + t154 * t193;
t149 = sin(qJ(3));
t153 = cos(qJ(3));
t146 = sin(pkin(6));
t175 = t146 * t194;
t121 = t136 * t153 - t149 * t175;
t135 = t150 * t193 - t154 * t172;
t148 = sin(qJ(4));
t152 = cos(qJ(4));
t107 = t121 * t148 - t135 * t152;
t108 = t121 * t152 + t135 * t148;
t147 = sin(qJ(6));
t151 = cos(qJ(6));
t203 = t107 * t151 - t108 * t147;
t202 = t107 * t147 + t108 * t151;
t199 = MDP(17) - MDP(26);
t177 = MDP(23) + MDP(25);
t198 = MDP(24) - MDP(27);
t171 = t185 * t193;
t138 = -t150 * t171 + t154 * t194;
t174 = t146 * t193;
t125 = t138 * t153 + t149 * t174;
t137 = t150 * t194 + t154 * t171;
t111 = t125 * t148 - t137 * t152;
t112 = t125 * t152 + t137 * t148;
t101 = t111 * t151 - t112 * t147;
t102 = t111 * t147 + t112 * t151;
t182 = t146 * t150;
t134 = t149 * t185 + t153 * t182;
t178 = t152 * t154;
t118 = t134 * t148 + t146 * t178;
t181 = t146 * t154;
t176 = t148 * t181;
t119 = t134 * t152 - t176;
t201 = -(g(2) * t203 + g(3) * (t118 * t151 - t119 * t147) + g(1) * t101) * MDP(34) + (g(2) * t202 + g(3) * (t118 * t147 + t119 * t151) + g(1) * t102) * MDP(35);
t200 = g(1) * t137 + g(2) * t135;
t180 = t148 * t153;
t179 = t152 * t153;
t173 = -t136 * t149 - t153 * t175;
t162 = pkin(3) * t153 + pkin(10) * t149 + pkin(2);
t160 = g(1) * t111 + g(2) * t107 + g(3) * t118;
t127 = (t148 * t150 + t153 * t178) * t146;
t126 = -t152 * t182 + t153 * t176;
t124 = t138 * t149 - t153 * t174;
t117 = -t137 * t179 + t138 * t148;
t116 = -t137 * t180 - t138 * t152;
t115 = -t135 * t179 + t136 * t148;
t114 = -t135 * t180 - t136 * t152;
t1 = [(g(1) * t193 - g(2) * t194) * MDP(2) + (g(1) * t194 + g(2) * t193) * MDP(3) + (g(1) * t136 - g(2) * t138) * MDP(9) + (-g(1) * t135 + g(2) * t137) * MDP(10) + (g(1) * t121 - g(2) * t125) * MDP(16) + (-g(1) * (-pkin(1) * t193 - t136 * pkin(2) - pkin(3) * t121 - pkin(4) * t108 + pkin(8) * t175 - t135 * pkin(9) + pkin(10) * t173 - qJ(5) * t107) - g(2) * (pkin(1) * t194 + t138 * pkin(2) + t125 * pkin(3) + t112 * pkin(4) + pkin(8) * t174 + t137 * pkin(9) + t124 * pkin(10) + t111 * qJ(5))) * MDP(28) + (g(1) * t202 - g(2) * t102) * MDP(34) + (g(1) * t203 - g(2) * t101) * MDP(35) + t177 * (g(1) * t108 - g(2) * t112) + t198 * (-g(1) * t107 + g(2) * t111) + t199 * (g(1) * t173 + g(2) * t124); (g(1) * t138 + g(2) * t136 + g(3) * t182) * MDP(10) + (-g(1) * (pkin(4) * t117 + pkin(9) * t138 + qJ(5) * t116) - g(2) * (pkin(4) * t115 + pkin(9) * t136 + qJ(5) * t114) + t200 * t162 + (-t127 * pkin(4) - t126 * qJ(5) - (pkin(9) * t150 + t154 * t162) * t146) * g(3)) * MDP(28) + (-g(1) * (t116 * t147 + t117 * t151) - g(2) * (t114 * t147 + t115 * t151) - g(3) * (t126 * t147 + t127 * t151)) * MDP(34) + (-g(1) * (t116 * t151 - t117 * t147) - g(2) * (t114 * t151 - t115 * t147) - g(3) * (t126 * t151 - t127 * t147)) * MDP(35) + t198 * (g(1) * t116 + g(2) * t114 + g(3) * t126) + t177 * (-g(1) * t117 - g(2) * t115 - g(3) * t127) + (-t153 * MDP(16) + t199 * t149 - MDP(9)) * (g(3) * t181 - t200); (-MDP(28) * pkin(10) + t199) * (g(1) * t125 + g(2) * t121 + g(3) * t134) + (-(pkin(4) * t152 + qJ(5) * t148 + pkin(3)) * MDP(28) + (t147 * t152 - t148 * t151) * MDP(35) - (t147 * t148 + t151 * t152) * MDP(34) + t198 * t148 - t177 * t152 - MDP(16)) * (g(3) * (-t149 * t182 + t153 * t185) + g(2) * t173 - g(1) * t124); (-g(1) * (-pkin(4) * t111 + qJ(5) * t112) - g(2) * (-pkin(4) * t107 + qJ(5) * t108) - g(3) * (-pkin(4) * t118 + qJ(5) * t119)) * MDP(28) + t177 * t160 + t198 * (g(1) * t112 + g(2) * t108 + g(3) * t119) - t201; -t160 * MDP(28); t201;];
taug  = t1;
