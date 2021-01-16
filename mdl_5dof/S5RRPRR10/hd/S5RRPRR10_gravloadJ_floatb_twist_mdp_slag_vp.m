% Calculate Gravitation load on the joints for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:02:46
% EndTime: 2021-01-15 22:02:51
% DurationCPUTime: 0.68s
% Computational Cost: add. (303->115), mult. (470->192), div. (0->0), fcn. (533->16), ass. (0->81)
t141 = cos(pkin(5));
t139 = qJ(2) + pkin(10);
t137 = sin(t139);
t149 = cos(qJ(1));
t171 = t149 * t137;
t138 = cos(t139);
t145 = sin(qJ(1));
t179 = t145 * t138;
t108 = t141 * t171 + t179;
t143 = sin(qJ(4));
t147 = cos(qJ(4));
t140 = sin(pkin(5));
t184 = t140 * t143;
t188 = g(3) * (-t137 * t184 + t141 * t147);
t127 = t149 * t138;
t180 = t145 * t137;
t193 = t141 * t180 - t127;
t195 = -(g(1) * t193 - g(2) * t108) * t143 - t188;
t136 = pkin(5) - t139;
t194 = sin(t136) / 0.2e1;
t135 = pkin(5) + t139;
t128 = sin(t135);
t190 = -t128 / 0.2e1 + t194;
t187 = g(3) * t140;
t186 = t138 * t140;
t142 = sin(qJ(5));
t185 = t138 * t142;
t183 = t140 * t147;
t182 = t142 * t145;
t181 = t142 * t149;
t144 = sin(qJ(2));
t178 = t145 * t144;
t177 = t145 * t147;
t148 = cos(qJ(2));
t176 = t145 * t148;
t146 = cos(qJ(5));
t175 = t146 * t145;
t174 = t146 * t147;
t173 = t146 * t149;
t172 = t147 * t149;
t170 = t149 * t144;
t169 = t149 * t148;
t167 = t145 * t184;
t165 = t145 * t174;
t164 = t142 * t177;
t126 = t149 * t184;
t163 = t142 * t172;
t162 = t146 * t172;
t161 = t108 * t147 - t126;
t159 = g(1) * t145 - g(2) * t149;
t130 = cos(t135);
t131 = cos(t136);
t125 = t130 + t131;
t158 = t180 - t149 * t125 / 0.2e1;
t157 = t171 + t145 * t125 / 0.2e1;
t156 = t159 * t140;
t155 = -t147 * t193 + t167;
t153 = t141 * t169 - t178;
t121 = -t141 * t176 - t170;
t110 = t141 * t165 - t181;
t118 = t141 * t182 + t162;
t151 = -t110 * t137 + t118 * t138 + t146 * t167;
t111 = t141 * t163 - t175;
t116 = -t141 * t173 - t164;
t150 = -t111 * t137 + t116 * t138 + t142 * t126;
t134 = t148 * pkin(2) + pkin(1);
t123 = t141 * t144 * pkin(2) - t140 * (pkin(7) + qJ(3));
t122 = -t141 * t178 + t169;
t120 = -t141 * t170 - t176;
t117 = t141 * t175 - t163;
t115 = t141 * t181 - t165;
t114 = t137 * t183 + t141 * t143;
t112 = t141 * t162 + t182;
t109 = t141 * t164 + t173;
t104 = t145 * t190 + t127;
t103 = t149 * t190 - t179;
t102 = t140 * t177 + t143 * t193;
t101 = t108 * t143 + t140 * t172;
t100 = t109 * t137 + t117 * t138 - t142 * t167;
t99 = -t112 * t137 + t115 * t138 + t146 * t126;
t1 = [t159 * MDP(2) + (-g(1) * t120 - g(2) * t122) * MDP(9) + (g(1) * t153 - g(2) * t121) * MDP(10) + (-g(1) * t103 - g(2) * t104) * MDP(11) + (-g(1) * t158 + g(2) * t157) * MDP(12) + (-g(1) * (-t123 * t149 - t145 * t134) - g(2) * (-t145 * t123 + t149 * t134)) * MDP(14) + (g(1) * t161 - g(2) * t155) * MDP(20) + (-g(1) * t101 - g(2) * t102) * MDP(21) + (-g(1) * t99 - g(2) * t151) * MDP(27) + (g(1) * t150 - g(2) * t100) * MDP(28) + (-t140 * MDP(13) + MDP(3)) * (g(1) * t149 + g(2) * t145); (g(1) * t122 - g(2) * t120 + t144 * t187) * MDP(10) + (g(1) * t157 + g(2) * t158 - g(3) * (t194 + t128 / 0.2e1)) * MDP(11) + (g(1) * t104 - g(2) * t103 - g(3) * (-t131 / 0.2e1 + t130 / 0.2e1)) * MDP(12) + (-g(1) * (-t110 * t138 - t118 * t137) - g(2) * (t112 * t138 + t115 * t137) - (t137 * t142 + t138 * t174) * t187) * MDP(27) + (-g(1) * (t109 * t138 - t117 * t137) - g(2) * (-t111 * t138 - t116 * t137) - (t137 * t146 - t147 * t185) * t187) * MDP(28) + (-t147 * MDP(20) + MDP(21) * t143) * (-g(1) * (t141 * t179 + t171) + g(2) * (t141 * t127 - t180) + g(3) * t186) + (pkin(2) * MDP(14) + MDP(9)) * (-g(1) * t121 - g(2) * t153 - t148 * t187); (-g(3) * t141 - t156) * MDP(14); (-g(1) * t102 + g(2) * t101 - t188) * MDP(20) + (g(1) * t155 + g(2) * t161 + g(3) * t114) * MDP(21) + (-t147 * t156 + t195) * t146 * MDP(27) + (t159 * t183 - t195) * t142 * MDP(28); (-g(1) * t100 - g(2) * t150 - g(3) * (-t114 * t142 - t146 * t186)) * MDP(27) + (g(1) * t151 - g(2) * t99 - g(3) * (-t114 * t146 + t140 * t185)) * MDP(28);];
taug = t1;
