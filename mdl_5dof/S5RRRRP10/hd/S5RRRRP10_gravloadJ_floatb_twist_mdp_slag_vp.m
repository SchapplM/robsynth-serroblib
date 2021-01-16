% Calculate Gravitation load on the joints for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:35:55
% EndTime: 2021-01-16 00:35:59
% DurationCPUTime: 0.87s
% Computational Cost: add. (366->114), mult. (827->196), div. (0->0), fcn. (946->10), ass. (0->65)
t141 = cos(pkin(5));
t144 = sin(qJ(3));
t145 = sin(qJ(2));
t174 = t144 * t145;
t140 = sin(pkin(5));
t148 = cos(qJ(3));
t178 = t140 * t148;
t119 = t141 * t174 + t178;
t150 = cos(qJ(1));
t146 = sin(qJ(1));
t149 = cos(qJ(2));
t166 = t150 * t149;
t163 = g(1) * (-t119 * t146 + t144 * t166) + g(3) * (t140 * t174 - t141 * t148);
t170 = t146 * t149;
t189 = t163 + g(2) * (t119 * t150 + t144 * t170);
t188 = MDP(17) - MDP(27);
t165 = MDP(23) + MDP(25);
t164 = MDP(24) + MDP(26);
t180 = g(3) * t140;
t179 = t140 * t145;
t177 = t141 * t144;
t143 = sin(qJ(4));
t176 = t143 * t149;
t175 = t144 * t140;
t173 = t145 * t146;
t172 = t145 * t148;
t171 = t145 * t150;
t147 = cos(qJ(4));
t169 = t147 * t149;
t168 = t148 * t149;
t130 = t143 * t145 + t147 * t168;
t167 = t150 * t130;
t139 = pkin(4) * t147 + pkin(3);
t142 = qJ(5) + pkin(9);
t123 = -t139 * t144 + t142 * t148;
t122 = t141 * t172 - t175;
t110 = t122 * t143 + t141 * t169;
t129 = t143 * t168 - t145 * t147;
t162 = t110 * t150 + t146 * t129;
t160 = t139 * t148 + t142 * t144;
t118 = pkin(2) + t160;
t138 = pkin(4) * t143 + pkin(8);
t161 = t118 * t149 + t138 * t145;
t109 = -t118 * t145 + t138 * t149;
t126 = -t141 * t166 + t173;
t159 = g(2) * t126 - t149 * t180;
t156 = -t122 * t147 + t141 * t176;
t124 = t141 * t173 - t166;
t114 = -t124 * t148 + t146 * t175;
t125 = -t141 * t171 - t170;
t113 = -t125 * t144 + t150 * t178;
t127 = t141 * t170 + t171;
t153 = t143 * t172 + t169;
t121 = t140 * t172 + t177;
t117 = t129 * pkin(4);
t116 = t127 * t148;
t115 = t129 * t141;
t112 = -t122 * t150 - t146 * t168;
t108 = pkin(1) + t161;
t107 = (t153 * t141 - t143 * t175) * pkin(4);
t106 = t161 * t141;
t105 = t123 * t141 * t145 - t140 * t160;
t104 = -t156 * t146 - t167;
t103 = -t140 * (-pkin(7) + t123) + t109 * t141;
t1 = [(g(1) * t146 - g(2) * t150) * MDP(2) + (g(1) * t150 + g(2) * t146) * MDP(3) + (-g(1) * t125 + g(2) * t124) * MDP(9) + (-g(1) * t126 + g(2) * t127) * MDP(10) + (-g(1) * t112 - g(2) * t114) * MDP(16) + (-g(1) * (t103 * t150 - t108 * t146) - g(2) * (t103 * t146 + t108 * t150)) * MDP(28) + t164 * (-g(1) * t162 - g(2) * (t110 * t146 - t150 * t129)) + t165 * (-g(1) * (t112 * t147 - t126 * t143) + g(2) * t104) - t188 * (g(2) * (t124 * t144 + t146 * t178) + g(1) * t113); (-g(1) * t124 - g(2) * t125 + g(3) * t179) * MDP(10) + (g(1) * t116 + t159 * t148) * MDP(16) + (-g(1) * (-t106 * t146 + t109 * t150) - g(2) * (t106 * t150 + t109 * t146) - t161 * t180) * MDP(28) + t165 * (-g(1) * (-t116 * t147 - t124 * t143) - g(2) * (t146 * (-t147 * t172 + t176) + t141 * t167) - t130 * t180) + t164 * (-g(1) * (t115 * t146 + t150 * t153) - g(2) * (-t115 * t150 + t146 * t153) + t129 * t180) + (t144 * t188 - MDP(9)) * (-g(1) * t127 - t159); (g(2) * t113 + t163) * MDP(16) + (-g(1) * (-t105 * t146 + t123 * t166) - g(2) * (t105 * t150 + t123 * t170) - g(3) * (t123 * t179 + t160 * t141)) * MDP(28) - t188 * (g(2) * (t125 * t148 + t150 * t175) - g(3) * t121 - g(1) * t114) + (-t164 * t143 + t165 * t147) * t189; (-g(1) * (t107 * t146 - t117 * t150) - g(2) * (-t107 * t150 - t117 * t146) - g(3) * (-t153 * t140 - t143 * t177) * pkin(4)) * MDP(28) + t165 * (-g(1) * (-(-t122 * t146 + t148 * t166) * t143 + t127 * t147) + g(2) * t162 - g(3) * (-t121 * t143 - t140 * t169)) + t164 * (-g(1) * t104 - g(2) * (-t146 * t130 + t156 * t150) - g(3) * (-t121 * t147 + t140 * t176)); -t189 * MDP(28);];
taug = t1;
