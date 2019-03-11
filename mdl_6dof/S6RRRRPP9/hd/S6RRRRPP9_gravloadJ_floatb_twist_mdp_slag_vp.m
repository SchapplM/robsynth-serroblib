% Calculate Gravitation load on the joints for
% S6RRRRPP9
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
%   see S6RRRRPP9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:49:49
% EndTime: 2019-03-09 21:49:52
% DurationCPUTime: 1.08s
% Computational Cost: add. (721->148), mult. (1819->221), div. (0->0), fcn. (2227->10), ass. (0->68)
t279 = MDP(24) - MDP(27) - MDP(30);
t278 = MDP(26) - MDP(23) - MDP(31);
t227 = cos(qJ(3));
t282 = -pkin(3) * t227 - pkin(2);
t223 = sin(qJ(4));
t281 = qJ(5) * t223 + pkin(3);
t280 = MDP(17) - MDP(25) - MDP(29);
t225 = sin(qJ(2));
t228 = cos(qJ(2));
t273 = cos(pkin(6));
t276 = cos(qJ(1));
t245 = t273 * t276;
t275 = sin(qJ(1));
t205 = t225 * t245 + t275 * t228;
t224 = sin(qJ(3));
t222 = sin(pkin(6));
t251 = t222 * t276;
t178 = t205 * t227 - t224 * t251;
t204 = t275 * t225 - t228 * t245;
t226 = cos(qJ(4));
t153 = t178 * t223 - t204 * t226;
t154 = t178 * t226 + t204 * t223;
t277 = pkin(5) + pkin(10);
t246 = -t205 * t224 - t227 * t251;
t271 = t246 * t226;
t244 = t273 * t275;
t207 = -t225 * t244 + t276 * t228;
t250 = t222 * t275;
t181 = t207 * t224 - t227 * t250;
t270 = t181 * t226;
t264 = t222 * t225;
t202 = -t224 * t264 + t273 * t227;
t269 = t202 * t226;
t267 = t204 * t224;
t206 = t276 * t225 + t228 * t244;
t265 = t206 * t224;
t263 = t222 * t228;
t262 = t223 * t227;
t261 = t226 * t227;
t260 = t227 * t228;
t256 = t224 * t263;
t255 = t223 * t263;
t254 = pkin(4) * t271 + t281 * t246;
t253 = -pkin(4) * t270 - t281 * t181;
t252 = pkin(4) * t269 + t281 * t202;
t249 = -t153 * pkin(4) + qJ(5) * t154;
t182 = t207 * t227 + t224 * t250;
t157 = t182 * t223 - t206 * t226;
t158 = t182 * t226 + t206 * t223;
t248 = -t157 * pkin(4) + qJ(5) * t158;
t203 = t273 * t224 + t227 * t264;
t175 = t203 * t223 + t226 * t263;
t176 = t203 * t226 - t255;
t247 = -t175 * pkin(4) + qJ(5) * t176;
t240 = g(1) * t157 + g(2) * t153 + g(3) * t175;
t239 = g(1) * t158 + g(2) * t154 + g(3) * t176;
t184 = -t226 * t264 + t227 * t255;
t185 = (t223 * t225 + t226 * t260) * t222;
t234 = t222 * pkin(3) * t260 + pkin(2) * t263 + t185 * pkin(4) + pkin(9) * t264 + pkin(10) * t256 + t184 * qJ(5);
t162 = -t204 * t262 - t205 * t226;
t163 = -t204 * t261 + t205 * t223;
t232 = t163 * pkin(4) + pkin(9) * t205 - pkin(10) * t267 + t162 * qJ(5) + t282 * t204;
t164 = -t206 * t262 - t207 * t226;
t165 = -t206 * t261 + t207 * t223;
t231 = t165 * pkin(4) + pkin(9) * t207 - pkin(10) * t265 + qJ(5) * t164 + t282 * t206;
t230 = t276 * pkin(1) + t207 * pkin(2) + t182 * pkin(3) + t158 * pkin(4) + pkin(8) * t250 + pkin(9) * t206 + qJ(5) * t157;
t229 = -t275 * pkin(1) - t205 * pkin(2) - pkin(3) * t178 - pkin(4) * t154 + pkin(8) * t251 - t204 * pkin(9) - qJ(5) * t153;
t1 = [(g(1) * t275 - g(2) * t276) * MDP(2) + (g(1) * t276 + g(2) * t275) * MDP(3) + (g(1) * t205 - g(2) * t207) * MDP(9) + (-g(1) * t204 + g(2) * t206) * MDP(10) + (g(1) * t178 - g(2) * t182) * MDP(16) + (-g(1) * (pkin(10) * t246 + t229) - g(2) * (pkin(10) * t181 + t230)) * MDP(28) + (-g(1) * (-qJ(6) * t154 + t246 * t277 + t229) - g(2) * (qJ(6) * t158 + t277 * t181 + t230)) * MDP(32) + t279 * (-g(1) * t153 + g(2) * t157) + t278 * (-g(1) * t154 + g(2) * t158) + t280 * (g(1) * t246 + g(2) * t181); (g(1) * t207 + g(2) * t205 + g(3) * t264) * MDP(10) + (-g(1) * t231 - g(2) * t232 - g(3) * t234) * MDP(28) + (-g(1) * (-pkin(5) * t265 + qJ(6) * t165 + t231) - g(2) * (-pkin(5) * t267 + qJ(6) * t163 + t232) - g(3) * (pkin(5) * t256 + t185 * qJ(6) + t234)) * MDP(32) + t279 * (g(1) * t164 + g(2) * t162 + g(3) * t184) + t278 * (g(1) * t165 + g(2) * t163 + g(3) * t185) + (-MDP(16) * t227 + t280 * t224 - MDP(9)) * (-g(1) * t206 - g(2) * t204 + g(3) * t263); (-g(1) * (pkin(10) * t182 + t253) - g(2) * (pkin(10) * t178 + t254) - g(3) * (pkin(10) * t203 + t252)) * MDP(28) + (-g(1) * (-qJ(6) * t270 + t277 * t182 + t253) - g(2) * (qJ(6) * t271 + t277 * t178 + t254) - g(3) * (qJ(6) * t269 + t277 * t203 + t252)) * MDP(32) + t280 * (g(1) * t182 + g(2) * t178 + g(3) * t203) + (-t279 * t223 - t278 * t226 + MDP(16)) * (g(1) * t181 - g(2) * t246 - g(3) * t202); (-g(1) * t248 - g(2) * t249 - g(3) * t247) * MDP(28) + (-g(1) * (-qJ(6) * t157 + t248) - g(2) * (-qJ(6) * t153 + t249) - g(3) * (-qJ(6) * t175 + t247)) * MDP(32) - t278 * t240 + t279 * t239; -(MDP(28) + MDP(32)) * t240; -t239 * MDP(32);];
taug  = t1;
