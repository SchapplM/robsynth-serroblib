% Calculate Gravitation load on the joints for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:09:31
% EndTime: 2019-03-09 16:09:34
% DurationCPUTime: 1.03s
% Computational Cost: add. (399->119), mult. (984->187), div. (0->0), fcn. (1152->10), ass. (0->55)
t190 = -MDP(16) - MDP(18) + MDP(23);
t146 = sin(qJ(3));
t194 = -qJ(4) * t146 - pkin(2);
t147 = sin(qJ(2));
t148 = sin(qJ(1));
t151 = cos(qJ(2));
t185 = cos(pkin(6));
t188 = cos(qJ(1));
t160 = t185 * t188;
t129 = t147 * t160 + t148 * t151;
t150 = cos(qJ(3));
t144 = sin(pkin(6));
t166 = t144 * t188;
t111 = t129 * t146 + t150 * t166;
t128 = t147 * t148 - t151 * t160;
t145 = sin(qJ(6));
t149 = cos(qJ(6));
t193 = t111 * t145 + t128 * t149;
t192 = t111 * t149 - t128 * t145;
t191 = MDP(17) - MDP(20) - MDP(22);
t189 = MDP(19) - MDP(10) - MDP(24);
t187 = g(3) * t144;
t186 = pkin(9) - qJ(5);
t181 = t128 * t150;
t165 = t148 * t185;
t130 = t147 * t188 + t151 * t165;
t180 = t130 * t150;
t179 = t144 * t147;
t178 = t144 * t148;
t177 = t144 * t150;
t176 = t144 * t151;
t175 = t145 * t146;
t174 = t146 * t149;
t173 = t146 * t151;
t172 = t149 * t151;
t171 = t150 * t151;
t168 = -pkin(3) * t181 + t128 * t194;
t167 = -pkin(3) * t180 + t130 * t194;
t112 = t129 * t150 - t146 * t166;
t164 = -t111 * pkin(3) + qJ(4) * t112;
t131 = -t147 * t165 + t151 * t188;
t115 = t131 * t146 - t148 * t177;
t116 = t131 * t150 + t146 * t178;
t163 = -t115 * pkin(3) + qJ(4) * t116;
t126 = t146 * t179 - t150 * t185;
t127 = t146 * t185 + t147 * t177;
t162 = -t126 * pkin(3) + qJ(4) * t127;
t161 = pkin(2) * t176 + pkin(9) * t179 + (pkin(3) * t171 + qJ(4) * t173) * t144;
t156 = t188 * pkin(1) + t131 * pkin(2) + t116 * pkin(3) + pkin(8) * t178 + qJ(4) * t115;
t155 = g(1) * t115 + g(2) * t111 + g(3) * t126;
t153 = -g(1) * t130 - g(2) * t128 + g(3) * t176;
t152 = -pkin(1) * t148 - t129 * pkin(2) - pkin(3) * t112 + pkin(8) * t166 - qJ(4) * t111;
t101 = t115 * t149 - t130 * t145;
t100 = -t115 * t145 - t130 * t149;
t1 = [(g(1) * t148 - g(2) * t188) * MDP(2) + (g(1) * t188 + g(2) * t148) * MDP(3) + (g(1) * t129 - g(2) * t131) * MDP(9) + (-g(1) * (-pkin(9) * t128 + t152) - g(2) * (pkin(9) * t130 + t156)) * MDP(21) + (-g(1) * (-pkin(4) * t112 - t128 * t186 + t152) - g(2) * (pkin(4) * t116 + t130 * t186 + t156)) * MDP(25) + (g(1) * t192 - g(2) * t101) * MDP(31) + (-g(1) * t193 - g(2) * t100) * MDP(32) + t191 * (-g(1) * t111 + g(2) * t115) + t190 * (-g(1) * t112 + g(2) * t116) + t189 * (g(1) * t128 - g(2) * t130); (-g(1) * (pkin(9) * t131 + t167) - g(2) * (pkin(9) * t129 + t168) - g(3) * t161) * MDP(21) + (-g(1) * (-pkin(4) * t180 + t131 * t186 + t167) - g(2) * (-pkin(4) * t181 + t129 * t186 + t168) - g(3) * ((pkin(4) * t171 - qJ(5) * t147) * t144 + t161)) * MDP(25) + (-g(1) * (-t130 * t174 - t131 * t145) - g(2) * (-t128 * t174 - t129 * t145) - (-t145 * t147 + t146 * t172) * t187) * MDP(31) + (-g(1) * (t130 * t175 - t131 * t149) - g(2) * (t128 * t175 - t129 * t149) - (-t145 * t173 - t147 * t149) * t187) * MDP(32) - t189 * (g(1) * t131 + g(2) * t129 + g(3) * t179) + (t191 * t146 + t190 * t150 - MDP(9)) * t153; (-g(1) * t163 - g(2) * t164 - g(3) * t162) * MDP(21) + (-g(1) * (-pkin(4) * t115 + t163) - g(2) * (-pkin(4) * t111 + t164) - g(3) * (-pkin(4) * t126 + t162)) * MDP(25) + (-MDP(31) * t149 + MDP(32) * t145 + t191) * (g(1) * t116 + g(2) * t112 + g(3) * t127) - t190 * t155; -(MDP(21) + MDP(25)) * t155; -t153 * MDP(25); (-g(1) * t100 + g(2) * t193 - g(3) * (-t126 * t145 + t144 * t172)) * MDP(31) + (g(1) * t101 + g(2) * t192 - g(3) * (-t126 * t149 - t145 * t176)) * MDP(32);];
taug  = t1;
