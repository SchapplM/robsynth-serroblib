% Calculate Gravitation load on the joints for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:38:31
% EndTime: 2019-03-09 21:38:34
% DurationCPUTime: 0.98s
% Computational Cost: add. (716->147), mult. (1806->221), div. (0->0), fcn. (2209->10), ass. (0->68)
t225 = cos(qJ(3));
t275 = -pkin(3) * t225 - pkin(2);
t221 = sin(qJ(4));
t274 = qJ(5) * t221 + pkin(3);
t273 = MDP(24) - MDP(27) - MDP(30);
t223 = sin(qJ(2));
t226 = cos(qJ(2));
t267 = cos(pkin(6));
t271 = cos(qJ(1));
t239 = t267 * t271;
t270 = sin(qJ(1));
t203 = t223 * t239 + t270 * t226;
t222 = sin(qJ(3));
t220 = sin(pkin(6));
t246 = t220 * t271;
t176 = t203 * t225 - t222 * t246;
t202 = t270 * t223 - t226 * t239;
t224 = cos(qJ(4));
t151 = t176 * t221 - t202 * t224;
t152 = t176 * t224 + t202 * t221;
t272 = MDP(17) - MDP(26) + MDP(31);
t244 = MDP(23) + MDP(25) + MDP(29);
t268 = pkin(10) - qJ(6);
t240 = -t203 * t222 - t225 * t246;
t265 = t240 * t224;
t238 = t267 * t270;
t205 = -t223 * t238 + t271 * t226;
t245 = t220 * t270;
t179 = t205 * t222 - t225 * t245;
t264 = t179 * t224;
t258 = t220 * t223;
t200 = -t222 * t258 + t267 * t225;
t263 = t200 * t224;
t261 = t202 * t222;
t204 = t271 * t223 + t226 * t238;
t259 = t204 * t222;
t257 = t220 * t226;
t256 = t221 * t225;
t255 = t224 * t225;
t254 = t225 * t226;
t251 = t222 * t257;
t250 = t221 * t257;
t249 = pkin(4) * t265 + t274 * t240;
t248 = -pkin(4) * t264 - t274 * t179;
t247 = pkin(4) * t263 + t274 * t200;
t243 = -t151 * pkin(4) + qJ(5) * t152;
t180 = t205 * t225 + t222 * t245;
t155 = t180 * t221 - t204 * t224;
t156 = t180 * t224 + t204 * t221;
t242 = -t155 * pkin(4) + qJ(5) * t156;
t201 = t267 * t222 + t225 * t258;
t173 = t201 * t221 + t224 * t257;
t174 = t201 * t224 - t250;
t241 = -t173 * pkin(4) + qJ(5) * t174;
t236 = g(1) * t155 + g(2) * t151 + g(3) * t173;
t144 = g(1) * t179 - g(2) * t240 - g(3) * t200;
t182 = -t224 * t258 + t225 * t250;
t183 = (t221 * t223 + t224 * t254) * t220;
t232 = t220 * pkin(3) * t254 + pkin(2) * t257 + t183 * pkin(4) + pkin(9) * t258 + pkin(10) * t251 + t182 * qJ(5);
t160 = -t202 * t256 - t203 * t224;
t161 = -t202 * t255 + t203 * t221;
t230 = t161 * pkin(4) + pkin(9) * t203 - pkin(10) * t261 + t160 * qJ(5) + t275 * t202;
t162 = -t204 * t256 - t205 * t224;
t163 = -t204 * t255 + t205 * t221;
t229 = t163 * pkin(4) + pkin(9) * t205 - pkin(10) * t259 + qJ(5) * t162 + t275 * t204;
t228 = t271 * pkin(1) + t205 * pkin(2) + t180 * pkin(3) + t156 * pkin(4) + pkin(8) * t245 + pkin(9) * t204 + qJ(5) * t155;
t227 = -t270 * pkin(1) - t203 * pkin(2) - pkin(3) * t176 - pkin(4) * t152 + pkin(8) * t246 - t202 * pkin(9) - qJ(5) * t151;
t1 = [(g(1) * t270 - g(2) * t271) * MDP(2) + (g(1) * t271 + g(2) * t270) * MDP(3) + (g(1) * t203 - g(2) * t205) * MDP(9) + (-g(1) * t202 + g(2) * t204) * MDP(10) + (g(1) * t176 - g(2) * t180) * MDP(16) + (-g(1) * (pkin(10) * t240 + t227) - g(2) * (pkin(10) * t179 + t228)) * MDP(28) + (-g(1) * (-pkin(5) * t152 + t240 * t268 + t227) - g(2) * (pkin(5) * t156 + t268 * t179 + t228)) * MDP(32) + t272 * (g(1) * t240 + g(2) * t179) + t244 * (g(1) * t152 - g(2) * t156) + t273 * (-g(1) * t151 + g(2) * t155); (g(1) * t205 + g(2) * t203 + g(3) * t258) * MDP(10) + (-g(1) * t229 - g(2) * t230 - g(3) * t232) * MDP(28) + (-g(1) * (pkin(5) * t163 + qJ(6) * t259 + t229) - g(2) * (pkin(5) * t161 + qJ(6) * t261 + t230) - g(3) * (t183 * pkin(5) - qJ(6) * t251 + t232)) * MDP(32) + t244 * (-g(1) * t163 - g(2) * t161 - g(3) * t183) + t273 * (g(1) * t162 + g(2) * t160 + g(3) * t182) + (-MDP(16) * t225 + t222 * t272 - MDP(9)) * (-g(1) * t204 - g(2) * t202 + g(3) * t257); (-g(1) * (pkin(10) * t180 + t248) - g(2) * (pkin(10) * t176 + t249) - g(3) * (pkin(10) * t201 + t247)) * MDP(28) + (-g(1) * (-pkin(5) * t264 + t268 * t180 + t248) - g(2) * (pkin(5) * t265 + t268 * t176 + t249) - g(3) * (pkin(5) * t263 + t268 * t201 + t247)) * MDP(32) + t272 * (g(1) * t180 + g(2) * t176 + g(3) * t201) + (-t221 * t273 + t244 * t224 + MDP(16)) * t144; (-g(1) * t242 - g(2) * t243 - g(3) * t241) * MDP(28) + (-g(1) * (-pkin(5) * t155 + t242) - g(2) * (-pkin(5) * t151 + t243) - g(3) * (-pkin(5) * t173 + t241)) * MDP(32) + t244 * t236 + t273 * (g(1) * t156 + g(2) * t152 + g(3) * t174); -(MDP(28) + MDP(32)) * t236; t144 * MDP(32);];
taug  = t1;
