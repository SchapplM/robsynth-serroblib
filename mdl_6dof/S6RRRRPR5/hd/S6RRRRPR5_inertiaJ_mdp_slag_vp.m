% Calculate joint inertia matrix for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR5_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:16:20
% EndTime: 2019-03-09 22:16:24
% DurationCPUTime: 1.11s
% Computational Cost: add. (1165->207), mult. (2101->272), div. (0->0), fcn. (2223->8), ass. (0->113)
t222 = sin(qJ(3));
t205 = pkin(2) * t222 + pkin(9);
t221 = sin(qJ(4));
t218 = t221 ^ 2;
t225 = cos(qJ(4));
t219 = t225 ^ 2;
t275 = t218 + t219;
t277 = t275 * t205;
t220 = sin(qJ(6));
t224 = cos(qJ(6));
t186 = t220 * t221 + t224 * t225;
t189 = -t220 * t225 + t221 * t224;
t279 = t189 * MDP(31) - t186 * MDP(32);
t306 = t221 * MDP(20) + t225 * MDP(21) - t279;
t223 = sin(qJ(2));
t292 = cos(qJ(3));
t293 = cos(qJ(2));
t187 = t222 * t223 - t292 * t293;
t190 = t222 * t293 + t292 * t223;
t207 = -t293 * pkin(2) - pkin(1);
t155 = t187 * pkin(3) - t190 * pkin(9) + t207;
t302 = pkin(8) + pkin(7);
t197 = t302 * t223;
t199 = t302 * t293;
t167 = -t222 * t197 + t292 * t199;
t145 = t221 * t155 + t225 * t167;
t281 = -t225 * t155 + t221 * t167;
t285 = t190 * t225;
t295 = -pkin(4) - pkin(5);
t136 = -pkin(10) * t285 + t295 * t187 + t281;
t183 = t187 * qJ(5);
t141 = t183 + t145;
t138 = pkin(10) * t190 * t221 + t141;
t236 = -(t136 * t220 + t138 * t224) * MDP(35) + (t136 * t224 - t138 * t220) * MDP(34);
t304 = -t281 * MDP(23) - t145 * MDP(24) - t236;
t303 = 0.2e1 * t190;
t301 = t225 * pkin(4) + t221 * qJ(5);
t300 = MDP(34) * t186 + MDP(35) * t189;
t299 = MDP(24) * t221;
t298 = 0.2e1 * t207;
t297 = 2 * MDP(26);
t296 = -2 * MDP(30);
t294 = pkin(9) - pkin(10);
t291 = pkin(9) * t187;
t206 = -t292 * pkin(2) - pkin(3);
t290 = pkin(3) - t206;
t289 = -pkin(10) + t205;
t288 = pkin(4) * MDP(28);
t287 = qJ(5) * t225;
t286 = t187 * t205;
t282 = t221 * t225;
t181 = t206 - t301;
t254 = pkin(3) + t301;
t278 = -t181 + t254;
t276 = t275 * pkin(9);
t274 = MDP(28) * t254;
t273 = MDP(28) * t221;
t149 = t186 * t190;
t272 = MDP(29) * t149;
t270 = MDP(34) * (qJ(5) * t220 - t224 * t295);
t268 = MDP(35) * (t224 * qJ(5) + t220 * t295);
t265 = t141 * MDP(28);
t142 = -pkin(4) * t187 + t281;
t264 = t142 * MDP(28);
t165 = t292 * t197 + t199 * t222;
t248 = pkin(4) * t221 - t287;
t146 = t248 * t190 + t165;
t263 = t146 * MDP(27);
t148 = t189 * t190;
t262 = t148 * MDP(32);
t261 = t149 * MDP(31);
t260 = t181 * MDP(28);
t257 = 0.2e1 * t293;
t256 = MDP(22) + MDP(33);
t255 = 0.2e1 * pkin(4) * MDP(25);
t253 = t294 * t221;
t252 = MDP(19) * t282;
t251 = t289 * t221;
t250 = MDP(28) * t275;
t249 = -pkin(3) * t190 - t291;
t247 = t190 * t254 + t291;
t246 = -t248 * MDP(26) + t306;
t245 = t218 * MDP(18) + MDP(15) + 0.2e1 * t252 + (MDP(29) * t189 + t186 * t296) * t189;
t244 = t181 * t190 - t286;
t243 = t190 * t206 - t286;
t242 = t225 * MDP(20) - t221 * MDP(21);
t240 = MDP(23) * t225 - t299;
t239 = -t165 * MDP(23) - t146 * MDP(25);
t238 = -0.2e1 * t225 * MDP(25) - 0.2e1 * t221 * MDP(27);
t237 = -t261 - t262;
t184 = t289 * t225;
t152 = t184 * t220 - t224 * t251;
t153 = t224 * t184 + t220 * t251;
t235 = t152 * MDP(34) + t153 * MDP(35);
t198 = t294 * t225;
t164 = t198 * t220 - t224 * t253;
t166 = t224 * t198 + t220 * t253;
t234 = t164 * MDP(34) + t166 * MDP(35);
t233 = MDP(34) * t224 - MDP(35) * t220;
t232 = -MDP(25) - t233;
t231 = -MDP(33) - t268 - t270;
t230 = 0.2e1 * t300;
t229 = (t292 * MDP(16) - t222 * MDP(17)) * pkin(2);
t228 = (MDP(28) * qJ(5) - MDP(24) + MDP(27)) * t225 + (-MDP(23) - MDP(25) - t288) * t221;
t227 = (t141 * t225 + t142 * t221) * MDP(26) + (t148 * t189 - t149 * t186) * MDP(30) + t189 * t272 - t167 * MDP(17) + (t299 - MDP(16)) * t165 + ((-t218 + t219) * MDP(19) + MDP(18) * t282 + MDP(13)) * t190 + (-MDP(14) + t306) * t187;
t216 = t225 * pkin(5);
t210 = t221 * MDP(26);
t182 = t216 + t254;
t169 = t216 - t181;
t143 = (t295 * t221 + t287) * t190 - t165;
t140 = t143 * t189;
t139 = t143 * t186;
t1 = [(t141 ^ 2 + t142 ^ 2 + t146 ^ 2) * MDP(28) + pkin(1) * MDP(9) * t257 + MDP(1) + (-t148 * t296 + t272) * t149 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t223 + MDP(5) * t257) * t223 + (MDP(16) * t298 - 0.2e1 * t261 - 0.2e1 * t262 + t256 * t187 + (-MDP(12) + t242) * t303) * t187 + 0.2e1 * (-t148 * MDP(34) + t149 * MDP(35)) * t143 + 0.2e1 * (-t142 * MDP(25) + t141 * MDP(27) + t304) * t187 + ((-t141 * t221 + t142 * t225) * MDP(26) + (t221 * MDP(23) + t225 * MDP(24)) * t165 + (t221 * MDP(25) - t225 * MDP(27)) * t146) * t303 + (MDP(17) * t298 + (t219 * MDP(18) + MDP(11) - 0.2e1 * t252) * t190) * t190; t227 + t223 * MDP(6) + (-t148 * t169 + t152 * t187 + t139) * MDP(34) + (t149 * t169 + t153 * t187 + t140) * MDP(35) + t146 * t260 + t293 * MDP(7) + (t243 * MDP(24) - t244 * MDP(27) + t205 * t265 + t239) * t225 + (t243 * MDP(23) + t244 * MDP(25) + t205 * t264 - t263) * t221 + (-t293 * MDP(10) - t223 * MDP(9)) * pkin(7); MDP(8) + t277 * t297 + t205 ^ 2 * t250 + (t238 + t260) * t181 + t169 * t230 + t245 - 0.2e1 * t240 * t206 + 0.2e1 * t229; (t249 * MDP(23) - t247 * MDP(25) + pkin(9) * t264 - t263) * t221 + (t249 * MDP(24) + t247 * MDP(27) + pkin(9) * t265 + t239) * t225 + t227 + (-t182 * t148 + t164 * t187 + t139) * MDP(34) + (t182 * t149 + t166 * t187 + t140) * MDP(35) - t146 * t274; (t276 + t277) * MDP(26) + (pkin(9) * t277 - t181 * t254) * MDP(28) + t229 + (t290 * MDP(23) + t278 * MDP(25)) * t225 + (-t290 * MDP(24) + t278 * MDP(27)) * t221 + t245 + t300 * (t169 + t182); t276 * t297 + pkin(9) ^ 2 * t250 - (t238 - t274) * t254 + t182 * t230 + 0.2e1 * t240 * pkin(3) + t245; -t281 * MDP(25) + (0.2e1 * t183 + t145) * MDP(27) + (-pkin(4) * t142 + qJ(5) * t141) * MDP(28) + (MDP(22) - t231 + t255) * t187 + (-MDP(26) * t301 + t242) * t190 + t237 + t304; t228 * t205 + t235 + t246; t228 * pkin(9) + t234 + t246; t255 + 0.2e1 * qJ(5) * MDP(27) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(28) + 0.2e1 * t270 + 0.2e1 * t268 + t256; MDP(26) * t285 + t232 * t187 + t264; t205 * t273 + t210; pkin(9) * t273 + t210; t232 - t288; MDP(28); -t187 * MDP(33) + t236 - t237; -t235 + t279; -t234 + t279; t231; t233; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
