% Calculate joint inertia matrix for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR4_inertiaJ_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:51:29
% EndTime: 2019-03-10 03:51:34
% DurationCPUTime: 1.28s
% Computational Cost: add. (2167->229), mult. (4324->304), div. (0->0), fcn. (4840->10), ass. (0->122)
t214 = sin(qJ(3));
t219 = cos(qJ(3));
t268 = pkin(8) + pkin(9);
t195 = t268 * t214;
t196 = t268 * t219;
t213 = sin(qJ(4));
t218 = cos(qJ(4));
t170 = -t218 * t195 - t196 * t213;
t171 = -t195 * t213 + t196 * t218;
t190 = t213 * t214 - t218 * t219;
t191 = t213 * t219 + t214 * t218;
t157 = -pkin(10) * t191 + t170;
t158 = -pkin(10) * t190 + t171;
t212 = sin(qJ(5));
t217 = cos(qJ(5));
t138 = t217 * t157 - t158 * t212;
t139 = t157 * t212 + t158 * t217;
t163 = t217 * t190 + t191 * t212;
t164 = -t190 * t212 + t191 * t217;
t128 = -pkin(11) * t164 + t138;
t129 = -pkin(11) * t163 + t139;
t211 = sin(qJ(6));
t216 = cos(qJ(6));
t144 = t216 * t163 + t164 * t211;
t145 = -t163 * t211 + t164 * t216;
t235 = t145 * MDP(34) - t144 * MDP(35) + (t128 * t216 - t129 * t211) * MDP(37) - (t128 * t211 + t129 * t216) * MDP(38);
t224 = t164 * MDP(27) - t163 * MDP(28) + t138 * MDP(30) - t139 * MDP(31) + t235;
t223 = t191 * MDP(20) - t190 * MDP(21) + t170 * MDP(23) - t171 * MDP(24) + t224;
t289 = t223 - (MDP(16) * t214 + MDP(17) * t219) * pkin(8) + t214 * MDP(13) + t219 * MDP(14);
t255 = MDP(38) * t211;
t226 = (MDP(37) * t216 - t255) * pkin(5);
t241 = t217 * MDP(30);
t288 = (-MDP(31) * t212 + t241) * pkin(4);
t240 = t218 * MDP(23);
t287 = (-MDP(24) * t213 + t240) * pkin(3);
t215 = sin(qJ(2));
t179 = t191 * t215;
t180 = t190 * t215;
t152 = t217 * t179 - t180 * t212;
t153 = -t179 * t212 - t180 * t217;
t133 = t216 * t152 + t153 * t211;
t134 = -t152 * t211 + t153 * t216;
t259 = t134 * MDP(34) - t133 * MDP(35);
t286 = t153 * MDP(27) - t152 * MDP(28);
t285 = -t180 * MDP(20) - t179 * MDP(21);
t202 = pkin(3) * t218 + pkin(4);
t198 = t217 * t202;
t267 = pkin(3) * t213;
t183 = -t212 * t267 + t198;
t181 = pkin(5) + t183;
t185 = t202 * t212 + t217 * t267;
t260 = t216 * t185;
t246 = (t211 * t181 + t260) * MDP(38);
t175 = t216 * t181;
t155 = -t185 * t211 + t175;
t247 = t155 * MDP(37);
t284 = t247 - t246;
t220 = cos(qJ(2));
t194 = -pkin(2) * t220 - t215 * pkin(8) - pkin(1);
t189 = t219 * t194;
t262 = pkin(9) * t215;
t265 = pkin(7) * t214;
t165 = -t219 * t262 + t189 + (-pkin(3) - t265) * t220;
t263 = pkin(7) * t220;
t238 = t219 * t263;
t167 = t238 + (t194 - t262) * t214;
t146 = t218 * t165 - t213 * t167;
t147 = t213 * t165 + t218 * t167;
t135 = -pkin(4) * t220 + t180 * pkin(10) + t146;
t140 = -pkin(10) * t179 + t147;
t126 = t217 * t135 - t212 * t140;
t127 = t212 * t135 + t217 * t140;
t124 = -pkin(5) * t220 - t153 * pkin(11) + t126;
t125 = -pkin(11) * t152 + t127;
t117 = t216 * t124 - t211 * t125;
t118 = t211 * t124 + t216 * t125;
t281 = t117 * MDP(37) - t118 * MDP(38) + t259;
t282 = t126 * MDP(30) - t127 * MDP(31) + t281 + t286;
t283 = t146 * MDP(23) - t147 * MDP(24) + t282 + t285;
t277 = -2 * MDP(19);
t276 = 0.2e1 * MDP(23);
t275 = 0.2e1 * MDP(24);
t274 = -2 * MDP(26);
t273 = 0.2e1 * MDP(30);
t272 = 0.2e1 * MDP(31);
t271 = -2 * MDP(33);
t270 = 0.2e1 * MDP(37);
t269 = 0.2e1 * MDP(38);
t266 = pkin(4) * t212;
t264 = pkin(7) * t219;
t261 = t214 * t219;
t193 = (pkin(3) * t214 + pkin(7)) * t215;
t258 = MDP(18) * t191;
t257 = MDP(25) * t164;
t256 = MDP(32) * t145;
t254 = MDP(38) * t216;
t206 = t217 * pkin(4);
t201 = t206 + pkin(5);
t197 = t216 * t201;
t245 = (-t211 * t266 + t197) * MDP(37);
t244 = t183 * MDP(30);
t243 = (t211 * t201 + t216 * t266) * MDP(38);
t242 = t185 * MDP(31);
t239 = MDP(29) + MDP(36);
t203 = -pkin(3) * t219 - pkin(2);
t237 = MDP(12) * t261;
t236 = MDP(22) + t239;
t166 = pkin(4) * t179 + t193;
t234 = MDP(15) + t236;
t176 = pkin(4) * t190 + t203;
t233 = MDP(13) * t219 - MDP(14) * t214;
t227 = MDP(36) - t243 + t245;
t225 = t239 - t242 + t244;
t209 = t219 ^ 2;
t208 = t215 ^ 2;
t207 = t214 ^ 2;
t205 = t216 * pkin(5);
t174 = t214 * t194 + t238;
t173 = -t214 * t263 + t189;
t148 = pkin(5) * t163 + t176;
t141 = pkin(5) * t152 + t166;
t1 = [-0.2e1 * pkin(1) * t215 * MDP(10) + MDP(1) - (-MDP(18) * t180 + t179 * t277) * t180 + (MDP(25) * t153 + t152 * t274) * t153 + (MDP(32) * t134 + t133 * t271) * t134 + (MDP(11) * t209 + MDP(4) - 0.2e1 * t237) * t208 + t234 * t220 ^ 2 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t233) * t215 - t285 - t286 - t259) * t220 + 0.2e1 * (t174 * t220 + t208 * t264) * MDP(17) + 0.2e1 * (-t173 * t220 + t208 * t265) * MDP(16) + (t147 * t220 - t193 * t180) * t275 + (-t146 * t220 + t193 * t179) * t276 + (-t126 * t220 + t166 * t152) * t273 + (t127 * t220 + t166 * t153) * t272 + (-t117 * t220 + t141 * t133) * t270 + (t118 * t220 + t141 * t134) * t269; -t180 * t258 + t153 * t257 + t134 * t256 + (-t133 * t145 - t134 * t144) * MDP(33) + (-t152 * t164 - t153 * t163) * MDP(26) + (-t179 * t191 + t180 * t190) * MDP(19) + (t148 * t133 + t141 * t144) * MDP(37) + (t148 * t134 + t141 * t145) * MDP(38) + (t176 * t152 + t166 * t163) * MDP(30) + (t176 * t153 + t166 * t164) * MDP(31) + (t203 * t179 + t193 * t190) * MDP(23) + (-t203 * t180 + t193 * t191) * MDP(24) + (-pkin(7) * MDP(9) + MDP(6) + (-t207 + t209) * MDP(12) + (-pkin(2) * t214 - t264) * MDP(16) + (-pkin(2) * t219 + t265) * MDP(17) + MDP(11) * t261) * t215 + (-pkin(7) * MDP(10) + MDP(7) - t289) * t220; 0.2e1 * t237 + t203 * t190 * t276 + t176 * t163 * t273 + t148 * t144 * t270 + t207 * MDP(11) + MDP(8) + 0.2e1 * (MDP(16) * t219 - MDP(17) * t214) * pkin(2) + (t190 * t277 + t203 * t275 + t258) * t191 + (t163 * t274 + t176 * t272 + t257) * t164 + (t144 * t271 + t148 * t269 + t256) * t145; (-MDP(15) - MDP(22) - t225 - t284 - t287) * t220 + t233 * t215 + t173 * MDP(16) - t174 * MDP(17) + t283; t289; t234 - 0.2e1 * t242 + 0.2e1 * t244 - 0.2e1 * t246 + 0.2e1 * t247 + 0.2e1 * t287; (-MDP(22) - MDP(29) - t227 - t288) * t220 + t283; t223; (t198 + t206) * MDP(30) + (t175 + t197) * MDP(37) - t185 * t254 + ((-pkin(4) - t202) * MDP(31) - pkin(4) * t254) * t212 + ((-t185 - t266) * MDP(37) + (-t181 - t201) * MDP(38)) * t211 + (t240 + (-MDP(30) * t212 - MDP(31) * t217 - MDP(24)) * t213) * pkin(3) + t236; t236 - 0.2e1 * t243 + 0.2e1 * t245 + 0.2e1 * t288; (-t239 - t226) * t220 + t282; t224; (t155 + t205) * MDP(37) + (-t260 + (-pkin(5) - t181) * t211) * MDP(38) + t225; (t197 + t205) * MDP(37) + (-pkin(5) - t201) * t255 + (t241 + (-MDP(37) * t211 - MDP(31) - t254) * t212) * pkin(4) + t239; 0.2e1 * t226 + t239; -t220 * MDP(36) + t281; t235; MDP(36) + t284; t227; MDP(36) + t226; MDP(36);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
