% Calculate joint inertia matrix for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR4_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:17:49
% EndTime: 2019-03-09 18:17:52
% DurationCPUTime: 1.00s
% Computational Cost: add. (1665->204), mult. (3172->269), div. (0->0), fcn. (3689->10), ass. (0->113)
t247 = sin(qJ(3));
t232 = pkin(2) * t247 + qJ(4);
t243 = sin(pkin(11));
t244 = cos(pkin(11));
t278 = t243 ^ 2 + t244 ^ 2;
t280 = t278 * t232;
t246 = sin(qJ(5));
t250 = cos(qJ(5));
t217 = t243 * t246 - t250 * t244;
t218 = t243 * t250 + t244 * t246;
t245 = sin(qJ(6));
t249 = cos(qJ(6));
t187 = t249 * t217 + t218 * t245;
t188 = -t217 * t245 + t218 * t249;
t283 = t188 * MDP(31) - t187 * MDP(32);
t269 = t218 * MDP(24) - t217 * MDP(25) + t283;
t248 = sin(qJ(2));
t293 = cos(qJ(3));
t294 = cos(qJ(2));
t220 = t247 * t248 - t293 * t294;
t221 = t247 * t294 + t293 * t248;
t236 = -t294 * pkin(2) - pkin(1);
t186 = t220 * pkin(3) - t221 * qJ(4) + t236;
t304 = pkin(8) + pkin(7);
t228 = t304 * t248;
t229 = t304 * t294;
t199 = -t247 * t228 + t293 * t229;
t162 = t244 * t186 - t199 * t243;
t240 = t244 * pkin(9);
t150 = pkin(4) * t220 - t221 * t240 + t162;
t163 = t243 * t186 + t244 * t199;
t284 = t221 * t243;
t157 = -pkin(9) * t284 + t163;
t138 = t250 * t150 - t157 * t246;
t139 = t150 * t246 + t157 * t250;
t175 = t218 * t221;
t176 = t217 * t221;
t306 = -t176 * MDP(24) - t175 * MDP(25) + t138 * MDP(27) - t139 * MDP(28);
t210 = (-pkin(9) - t232) * t243;
t211 = t232 * t244 + t240;
t179 = t250 * t210 - t211 * t246;
t291 = pkin(10) * t218;
t167 = t179 - t291;
t180 = t210 * t246 + t211 * t250;
t212 = t217 * pkin(10);
t168 = t180 - t212;
t146 = t167 * t249 - t168 * t245;
t147 = t167 * t245 + t168 * t249;
t303 = t146 * MDP(34) - t147 * MDP(35);
t224 = (-pkin(9) - qJ(4)) * t243;
t225 = qJ(4) * t244 + t240;
t194 = t250 * t224 - t225 * t246;
t172 = t194 - t291;
t195 = t224 * t246 + t225 * t250;
t173 = t195 - t212;
t153 = t172 * t249 - t173 * t245;
t154 = t172 * t245 + t173 * t249;
t302 = t153 * MDP(34) - t154 * MDP(35);
t158 = t249 * t175 - t176 * t245;
t159 = -t175 * t245 - t176 * t249;
t301 = t159 * MDP(31) - t158 * MDP(32);
t300 = t187 * MDP(34) + t188 * MDP(35);
t299 = t217 * MDP(27) + t218 * MDP(28);
t298 = t244 * MDP(18) - t243 * MDP(19);
t297 = 2 * MDP(20);
t296 = -2 * MDP(23);
t295 = -2 * MDP(30);
t292 = pkin(5) * t220;
t290 = t217 * pkin(5);
t289 = t244 * pkin(4);
t287 = pkin(3) * MDP(21);
t137 = -pkin(10) * t175 + t139;
t286 = t137 * t249;
t198 = t293 * t228 + t247 * t229;
t285 = t198 * t244;
t279 = t278 * qJ(4);
t277 = MDP(22) * t176;
t276 = MDP(29) * t159;
t235 = -t293 * pkin(2) - pkin(3);
t275 = t235 * MDP(21);
t273 = 0.2e1 * t294;
t272 = MDP(26) + MDP(33);
t271 = t220 * MDP(33) + t301;
t233 = -pkin(3) - t289;
t135 = pkin(10) * t176 + t138 + t292;
t132 = t249 * t135 - t137 * t245;
t270 = t278 * MDP(21);
t268 = -pkin(3) * t221 - qJ(4) * t220;
t267 = MDP(15) + (MDP(22) * t218 + t217 * t296) * t218 + (MDP(29) * t188 + t187 * t295) * t188;
t266 = -t162 * t243 + t163 * t244;
t265 = -t220 * t232 + t221 * t235;
t264 = 0.2e1 * t298;
t263 = t243 * MDP(18) + t244 * MDP(19);
t260 = t175 * MDP(27) - t176 * MDP(28);
t259 = t132 * MDP(34) - (t135 * t245 + t286) * MDP(35);
t258 = t158 * MDP(34) + t159 * MDP(35);
t171 = pkin(4) * t284 + t198;
t223 = t235 - t289;
t257 = (MDP(34) * t249 - MDP(35) * t245) * pkin(5);
t256 = 0.2e1 * t299;
t255 = 0.2e1 * t300;
t254 = -t298 + t299 + t300;
t253 = (t293 * MDP(16) - t247 * MDP(17)) * pkin(2);
t252 = (-t158 * t188 - t159 * t187) * MDP(30) + t266 * MDP(20) + t188 * t276 + (-t175 * t218 + t176 * t217) * MDP(23) - t218 * t277 + t221 * MDP(13) - t198 * MDP(16) - t199 * MDP(17) + (-MDP(14) + t269) * t220;
t201 = t233 + t290;
t200 = t223 + t290;
t191 = t198 * t243;
t165 = t171 * t218;
t164 = t171 * t217;
t160 = t175 * pkin(5) + t171;
t143 = t160 * t188;
t142 = t160 * t187;
t1 = [pkin(1) * MDP(9) * t273 + MDP(1) + (t162 ^ 2 + t163 ^ 2 + t198 ^ 2) * MDP(21) + (MDP(11) * t221 + 0.2e1 * t236 * MDP(17)) * t221 + t272 * t220 ^ 2 - (t175 * t296 - t277) * t176 + (t158 * t295 + t276) * t159 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t248 + MDP(5) * t273) * t248 + 0.2e1 * t260 * t171 + 0.2e1 * t258 * t160 + 0.2e1 * ((-t162 * t244 - t163 * t243) * MDP(20) + t263 * t198) * t221 + 0.2e1 * (-t221 * MDP(12) + t236 * MDP(16) + t162 * MDP(18) - t163 * MDP(19) + t259 + t301 + t306) * t220; t252 + t248 * MDP(6) + (t198 * t235 + t232 * t266) * MDP(21) + (t243 * t265 - t285) * MDP(18) + (t244 * t265 + t191) * MDP(19) + t294 * MDP(7) + (-t294 * MDP(10) - t248 * MDP(9)) * pkin(7) + (t146 * t220 + t158 * t200 + t142) * MDP(34) + (-t147 * t220 + t159 * t200 + t143) * MDP(35) + (t175 * t223 + t179 * t220 + t164) * MDP(27) + (-t176 * t223 - t180 * t220 + t165) * MDP(28); MDP(8) + t280 * t297 + t232 ^ 2 * t270 + (-t264 + t275) * t235 + t223 * t256 + t200 * t255 + 0.2e1 * t253 + t267; (t244 * t268 + t191) * MDP(19) + t252 + (-pkin(3) * t198 + qJ(4) * t266) * MDP(21) + (t243 * t268 - t285) * MDP(18) + (t175 * t233 + t194 * t220 + t164) * MDP(27) + (-t176 * t233 - t195 * t220 + t165) * MDP(28) + (t153 * t220 + t158 * t201 + t142) * MDP(34) + (-t154 * t220 + t159 * t201 + t143) * MDP(35); (t279 + t280) * MDP(20) + (-pkin(3) * t235 + qJ(4) * t280) * MDP(21) + t253 + t267 + t298 * (pkin(3) - t235) + t300 * (t200 + t201) + t299 * (t223 + t233); t279 * t297 + qJ(4) ^ 2 * t270 + t233 * t256 + t201 * t255 + (t264 + t287) * pkin(3) + t267; t198 * MDP(21) + t221 * t263 + t258 + t260; t254 + t275; t254 - t287; MDP(21); t220 * MDP(26) + (t249 * t292 + t132) * MDP(34) + (-t286 + (-t135 - t292) * t245) * MDP(35) + t271 + t306; t179 * MDP(27) - t180 * MDP(28) + t269 + t303; t194 * MDP(27) - t195 * MDP(28) + t269 + t302; 0; 0.2e1 * t257 + t272; t259 + t271; t283 + t303; t283 + t302; 0; MDP(33) + t257; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
