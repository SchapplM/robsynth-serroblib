% Calculate Gravitation load on the joints for
% S6RPRRRP12
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
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:50:41
% EndTime: 2019-03-09 06:50:44
% DurationCPUTime: 1.32s
% Computational Cost: add. (1050->131), mult. (2891->208), div. (0->0), fcn. (3733->14), ass. (0->63)
t191 = cos(qJ(1));
t227 = sin(pkin(12));
t232 = cos(pkin(6));
t216 = t232 * t227;
t230 = cos(pkin(12));
t239 = sin(qJ(1));
t176 = t191 * t216 + t230 * t239;
t188 = sin(qJ(3));
t240 = cos(qJ(3));
t217 = t232 * t230;
t199 = -t191 * t217 + t239 * t227;
t228 = sin(pkin(7));
t229 = sin(pkin(6));
t214 = t229 * t228;
t231 = cos(pkin(7));
t250 = t191 * t214 + t199 * t231;
t161 = -t176 * t240 + t250 * t188;
t187 = sin(qJ(4));
t190 = cos(qJ(4));
t215 = t231 * t229;
t241 = -t191 * t215 + t199 * t228;
t146 = t161 * t190 - t187 * t241;
t158 = t176 * t188 + t250 * t240;
t186 = sin(qJ(5));
t189 = cos(qJ(5));
t133 = t146 * t186 + t158 * t189;
t134 = t146 * t189 - t158 * t186;
t145 = t161 * t187 + t190 * t241;
t244 = MDP(21) - MDP(30);
t243 = MDP(27) + MDP(29);
t242 = MDP(28) - MDP(31);
t195 = t191 * t227 + t217 * t239;
t246 = t195 * t231 - t239 * t214;
t245 = t230 * t215 + t228 * t232;
t226 = t186 * t190;
t225 = t189 * t190;
t221 = t229 * t239;
t224 = t191 * pkin(1) + qJ(2) * t221;
t223 = t191 * t229;
t222 = -pkin(1) * t239 + qJ(2) * t223;
t213 = t229 * t227;
t177 = t191 * t230 - t216 * t239;
t163 = t177 * t240 - t188 * t246;
t192 = -t195 * t228 - t215 * t239;
t148 = t163 * t190 - t187 * t192;
t162 = t177 * t188 + t240 * t246;
t135 = t148 * t186 - t162 * t189;
t169 = t188 * t245 + t240 * t213;
t175 = -t214 * t230 + t231 * t232;
t157 = t169 * t190 + t175 * t187;
t168 = t188 * t213 - t240 * t245;
t141 = t157 * t186 - t168 * t189;
t207 = g(1) * t135 - g(2) * t133 + g(3) * t141;
t150 = -t168 * t225 + t169 * t186;
t149 = -t168 * t226 - t169 * t189;
t147 = t163 * t187 + t190 * t192;
t142 = t157 * t189 + t168 * t186;
t140 = -t162 * t225 + t163 * t186;
t139 = -t162 * t226 - t163 * t189;
t138 = -t158 * t225 - t161 * t186;
t137 = -t158 * t226 + t161 * t189;
t136 = t148 * t189 + t162 * t186;
t1 = [(g(1) * t239 - g(2) * t191) * MDP(2) + (g(1) * t191 + g(2) * t239) * MDP(3) + (g(1) * t176 - g(2) * t177) * MDP(4) + (-g(1) * t199 + g(2) * t195) * MDP(5) + (-g(1) * t223 - g(2) * t221) * MDP(6) + (-g(1) * t222 - g(2) * t224) * MDP(7) + (-g(1) * t161 - g(2) * t163) * MDP(13) + (-g(1) * t158 + g(2) * t162) * MDP(14) + (-g(1) * t146 - g(2) * t148) * MDP(20) + (-g(1) * (-t176 * pkin(2) + t161 * pkin(3) + t146 * pkin(4) + t134 * pkin(5) - pkin(10) * t158 + t145 * pkin(11) + t133 * qJ(6) + t222) - g(2) * (t177 * pkin(2) + t163 * pkin(3) + t148 * pkin(4) + t136 * pkin(5) + t162 * pkin(10) + t147 * pkin(11) + t135 * qJ(6) + t224) + (g(1) * t241 + g(2) * t192) * pkin(9)) * MDP(32) + t242 * (g(1) * t133 + g(2) * t135) + t244 * (g(1) * t145 + g(2) * t147) + t243 * (-g(1) * t134 - g(2) * t136); (MDP(32) + MDP(7)) * (-g(1) * t221 + g(2) * t223 - g(3) * t232); (g(1) * t163 - g(2) * t161 + g(3) * t169) * MDP(14) + (-g(1) * (pkin(5) * t140 + pkin(10) * t163 + qJ(6) * t139) - g(2) * (pkin(5) * t138 - pkin(10) * t161 + qJ(6) * t137) - g(3) * (pkin(5) * t150 + pkin(10) * t169 + qJ(6) * t149)) * MDP(32) + t242 * (g(1) * t139 + g(2) * t137 + g(3) * t149) + t243 * (-g(1) * t140 - g(2) * t138 - g(3) * t150) + (MDP(13) + t190 * MDP(20) + (pkin(4) * t190 + pkin(11) * t187 + pkin(3)) * MDP(32) - t244 * t187) * (g(1) * t162 + g(2) * t158 + g(3) * t168); (-MDP(32) * pkin(11) + t244) * (g(1) * t148 - g(2) * t146 + g(3) * t157) + (MDP(20) + MDP(32) * (pkin(5) * t189 + qJ(6) * t186 + pkin(4)) + t243 * t189 - t242 * t186) * (-g(3) * (-t169 * t187 + t175 * t190) - g(2) * t145 + g(1) * t147); (-g(1) * (-pkin(5) * t135 + qJ(6) * t136) - g(2) * (pkin(5) * t133 - qJ(6) * t134) - g(3) * (-pkin(5) * t141 + qJ(6) * t142)) * MDP(32) + t243 * t207 + t242 * (g(1) * t136 - g(2) * t134 + g(3) * t142); -t207 * MDP(32);];
taug  = t1;
