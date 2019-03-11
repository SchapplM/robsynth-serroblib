% Calculate Gravitation load on the joints for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP14_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP14_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:09:52
% EndTime: 2019-03-09 13:09:54
% DurationCPUTime: 0.81s
% Computational Cost: add. (455->124), mult. (1135->185), div. (0->0), fcn. (1359->10), ass. (0->56)
t208 = MDP(9) - MDP(12);
t207 = -MDP(10) + MDP(13);
t205 = MDP(21) - MDP(30);
t187 = MDP(27) + MDP(29);
t204 = MDP(28) - MDP(31);
t160 = sin(qJ(2));
t161 = sin(qJ(1));
t164 = cos(qJ(2));
t165 = cos(qJ(1));
t198 = cos(pkin(6));
t184 = t165 * t198;
t144 = t160 * t184 + t161 * t164;
t185 = t161 * t198;
t146 = -t160 * t185 + t164 * t165;
t206 = -g(1) * t146 - g(2) * t144;
t158 = sin(qJ(5));
t162 = cos(qJ(5));
t143 = t160 * t161 - t164 * t184;
t159 = sin(qJ(4));
t163 = cos(qJ(4));
t157 = sin(pkin(6));
t192 = t157 * t165;
t174 = -t143 * t159 + t163 * t192;
t109 = t144 * t162 + t158 * t174;
t110 = -t144 * t158 + t162 * t174;
t195 = t157 * t160;
t194 = t157 * t161;
t193 = t157 * t164;
t191 = t158 * t159;
t190 = t159 * t162;
t189 = t160 * t162;
t188 = pkin(2) * t193 + qJ(3) * t195;
t186 = t158 * t195;
t183 = pkin(4) * t159 - pkin(10) * t163;
t145 = t165 * t160 + t164 * t185;
t175 = t165 * pkin(1) + t146 * pkin(2) + pkin(8) * t194 + qJ(3) * t145;
t173 = t143 * t163 + t159 * t192;
t122 = t145 * t159 + t163 * t194;
t107 = t122 * t158 - t146 * t162;
t142 = -t159 * t193 + t163 * t198;
t119 = t142 * t158 - t157 * t189;
t172 = g(1) * t107 - g(2) * t109 + g(3) * t119;
t167 = -t161 * pkin(1) - t144 * pkin(2) + pkin(8) * t192 - t143 * qJ(3);
t114 = -g(1) * t145 - g(2) * t143 + g(3) * t193;
t139 = t145 * pkin(2);
t137 = t143 * pkin(2);
t128 = (t158 * t164 + t159 * t189) * t157;
t127 = t159 * t186 - t162 * t193;
t121 = -t145 * t163 + t159 * t194;
t120 = t142 * t162 + t186;
t118 = -t145 * t158 + t146 * t190;
t117 = t145 * t162 + t146 * t191;
t116 = -t143 * t158 + t144 * t190;
t115 = t143 * t162 + t144 * t191;
t108 = t122 * t162 + t146 * t158;
t1 = [(g(1) * t161 - g(2) * t165) * MDP(2) + (-g(1) * t167 - g(2) * t175) * MDP(14) + (-g(1) * t174 - g(2) * t122) * MDP(20) + (-g(1) * (pkin(3) * t192 + pkin(4) * t174 + t110 * pkin(5) - t144 * pkin(9) + pkin(10) * t173 + t109 * qJ(6) + t167) - g(2) * (pkin(3) * t194 + pkin(4) * t122 + pkin(5) * t108 + pkin(9) * t146 + pkin(10) * t121 + qJ(6) * t107 + t175)) * MDP(32) + t187 * (-g(1) * t110 - g(2) * t108) + t204 * (g(1) * t109 + g(2) * t107) + t205 * (g(1) * t173 + g(2) * t121) + t207 * (g(1) * t143 - g(2) * t145) + t208 * (g(1) * t144 - g(2) * t146) + (-t157 * MDP(11) + MDP(3)) * (g(1) * t165 + g(2) * t161); (-g(1) * (qJ(3) * t146 - t139) - g(2) * (qJ(3) * t144 - t137) - g(3) * t188) * MDP(14) + (-g(1) * (pkin(5) * t118 - pkin(9) * t145 + qJ(6) * t117 - t139) - g(2) * (pkin(5) * t116 - t143 * pkin(9) + qJ(6) * t115 - t137) + t206 * (qJ(3) + t183) + (-pkin(5) * t128 - t127 * qJ(6) - t188 - (pkin(9) * t164 + t160 * t183) * t157) * g(3)) * MDP(32) + t204 * (g(1) * t117 + g(2) * t115 + g(3) * t127) - t208 * t114 + t187 * (-g(1) * t118 - g(2) * t116 - g(3) * t128) + (-MDP(20) * t159 - t205 * t163 - t207) * (g(3) * t195 - t206); (MDP(14) + MDP(32)) * t114; (-pkin(10) * MDP(32) + t205) * (g(1) * t122 - g(2) * t174 + g(3) * t142) + (MDP(20) + MDP(32) * (pkin(5) * t162 + qJ(6) * t158 + pkin(4)) + t187 * t162 - t204 * t158) * (-g(3) * (-t159 * t198 - t163 * t193) - g(2) * t173 + g(1) * t121); (-g(1) * (-pkin(5) * t107 + qJ(6) * t108) - g(2) * (pkin(5) * t109 - qJ(6) * t110) - g(3) * (-pkin(5) * t119 + qJ(6) * t120)) * MDP(32) + t187 * t172 + t204 * (g(1) * t108 - g(2) * t110 + g(3) * t120); -t172 * MDP(32);];
taug  = t1;
