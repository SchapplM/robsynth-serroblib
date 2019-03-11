% Calculate Gravitation load on the joints for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:41:25
% EndTime: 2019-03-09 06:41:29
% DurationCPUTime: 1.17s
% Computational Cost: add. (682->119), mult. (1838->193), div. (0->0), fcn. (2337->14), ass. (0->58)
t187 = sin(pkin(12));
t192 = cos(pkin(6));
t173 = t192 * t187;
t190 = cos(pkin(12));
t198 = sin(qJ(1));
t200 = cos(qJ(1));
t134 = t200 * t173 + t198 * t190;
t148 = sin(qJ(3));
t199 = cos(qJ(3));
t175 = t192 * t190;
t161 = -t200 * t175 + t198 * t187;
t188 = sin(pkin(7));
t189 = sin(pkin(6));
t171 = t189 * t188;
t191 = cos(pkin(7));
t210 = t161 * t191 + t200 * t171;
t121 = -t134 * t199 + t210 * t148;
t147 = sin(qJ(4));
t150 = cos(qJ(4));
t172 = t191 * t189;
t202 = t161 * t188 - t200 * t172;
t113 = t121 * t150 - t147 * t202;
t118 = t134 * t148 + t210 * t199;
t146 = sin(qJ(5));
t149 = cos(qJ(5));
t213 = t113 * t146 + t118 * t149;
t186 = t118 * t146;
t212 = t113 * t149 - t186;
t112 = t121 * t147 + t150 * t202;
t203 = MDP(21) - MDP(29);
t157 = t198 * t175 + t200 * t187;
t206 = t157 * t191 - t198 * t171;
t204 = t190 * t172 + t192 * t188;
t135 = -t198 * t173 + t200 * t190;
t122 = t135 * t148 + t199 * t206;
t185 = t122 * t146;
t184 = t146 * t150;
t183 = t149 * t150;
t178 = t189 * t198;
t182 = t200 * pkin(1) + qJ(2) * t178;
t179 = t200 * t189;
t180 = -t198 * pkin(1) + qJ(2) * t179;
t170 = t189 * t187;
t123 = t135 * t199 - t148 * t206;
t151 = -t157 * t188 - t198 * t172;
t115 = t123 * t150 - t151 * t147;
t108 = -t115 * t146 + t122 * t149;
t114 = t123 * t147 + t151 * t150;
t128 = t148 * t204 + t199 * t170;
t156 = -t190 * t171 + t192 * t191;
t116 = t128 * t147 - t156 * t150;
t168 = g(1) * t114 - g(2) * t112 + g(3) * t116;
t145 = -qJ(6) - pkin(11);
t143 = pkin(5) * t149 + pkin(4);
t127 = t148 * t170 - t199 * t204;
t117 = t128 * t150 + t156 * t147;
t109 = t115 * t149 + t185;
t1 = [(g(1) * t198 - g(2) * t200) * MDP(2) + (g(1) * t200 + g(2) * t198) * MDP(3) + (g(1) * t134 - g(2) * t135) * MDP(4) + (-g(1) * t161 + g(2) * t157) * MDP(5) + (-g(1) * t179 - g(2) * t178) * MDP(6) + (-g(1) * t180 - g(2) * t182) * MDP(7) + (-g(1) * t121 - g(2) * t123) * MDP(13) + (-g(1) * t118 + g(2) * t122) * MDP(14) + (-g(1) * t113 - g(2) * t115) * MDP(20) + (-g(1) * t212 - g(2) * t109) * MDP(27) + (g(1) * t213 - g(2) * t108) * MDP(28) + (-g(1) * (-t134 * pkin(2) + t121 * pkin(3) - pkin(5) * t186 - pkin(10) * t118 - t112 * t145 + t113 * t143 + t180) - g(2) * (t135 * pkin(2) + t123 * pkin(3) + pkin(5) * t185 + t122 * pkin(10) - t114 * t145 + t115 * t143 + t182) + (g(1) * t202 + g(2) * t151) * pkin(9)) * MDP(30) + t203 * (g(1) * t112 + g(2) * t114); (MDP(30) + MDP(7)) * (-g(1) * t178 + g(2) * t179 - g(3) * t192); (-g(1) * (-t122 * t183 + t123 * t146) - g(2) * (-t118 * t183 - t121 * t146) - g(3) * (-t127 * t183 + t128 * t146)) * MDP(27) + (-g(1) * (t122 * t184 + t123 * t149) - g(2) * (t118 * t184 - t121 * t149) - g(3) * (t127 * t184 + t128 * t149)) * MDP(28) + (MDP(14) - (pkin(5) * t146 + pkin(10)) * MDP(30)) * (g(1) * t123 - g(2) * t121 + g(3) * t128) + (MDP(13) + t150 * MDP(20) + (t143 * t150 - t145 * t147 + pkin(3)) * MDP(30) - t203 * t147) * (g(1) * t122 + g(2) * t118 + g(3) * t127); (-g(1) * (-t114 * t143 - t115 * t145) - g(2) * (t112 * t143 + t113 * t145) - g(3) * (-t116 * t143 - t117 * t145)) * MDP(30) + t203 * (g(1) * t115 - g(2) * t113 + g(3) * t117) + (MDP(27) * t149 - MDP(28) * t146 + MDP(20)) * t168; (g(1) * t109 - g(2) * t212 - g(3) * (-t117 * t149 - t127 * t146)) * MDP(28) + (MDP(30) * pkin(5) + MDP(27)) * (-g(1) * t108 - g(2) * t213 - g(3) * (-t117 * t146 + t127 * t149)); -t168 * MDP(30);];
taug  = t1;
