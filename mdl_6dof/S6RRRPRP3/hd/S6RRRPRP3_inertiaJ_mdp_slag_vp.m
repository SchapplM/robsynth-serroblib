% Calculate joint inertia matrix for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP3_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:41:34
% EndTime: 2019-03-09 16:41:36
% DurationCPUTime: 1.01s
% Computational Cost: add. (1755->226), mult. (3182->296), div. (0->0), fcn. (3519->8), ass. (0->106)
t211 = sin(qJ(3));
t197 = pkin(2) * t211 + qJ(4);
t208 = sin(pkin(10));
t209 = cos(pkin(10));
t246 = t208 ^ 2 + t209 ^ 2;
t248 = t246 * t197;
t210 = sin(qJ(5));
t266 = cos(qJ(5));
t185 = t208 * t266 + t209 * t210;
t274 = -t208 * t210 + t209 * t266;
t278 = MDP(24) * t185 + MDP(25) * t274;
t277 = pkin(8) + pkin(7);
t276 = -MDP(27) * t274 + MDP(28) * t185;
t275 = -MDP(29) * t274 - MDP(31) * t185;
t273 = MDP(18) * t209 - MDP(19) * t208;
t268 = cos(qJ(2));
t201 = -pkin(2) * t268 - pkin(1);
t272 = 0.2e1 * t201;
t271 = 2 * MDP(20);
t270 = -2 * MDP(23);
t269 = 2 * MDP(30);
t267 = cos(qJ(3));
t212 = sin(qJ(2));
t186 = t211 * t212 - t267 * t268;
t265 = t186 * pkin(5);
t264 = t209 * pkin(4);
t205 = t209 * pkin(9);
t262 = pkin(3) * MDP(21);
t180 = t197 * t209 + t205;
t230 = (-pkin(9) - t197) * t208;
t153 = t180 * t210 - t230 * t266;
t261 = t153 * t186;
t154 = t180 * t266 + t210 * t230;
t260 = t154 * t186;
t190 = qJ(4) * t209 + t205;
t231 = (-pkin(9) - qJ(4)) * t208;
t163 = t190 * t210 - t231 * t266;
t259 = t163 * t186;
t164 = t190 * t266 + t210 * t231;
t258 = t164 * t186;
t193 = t277 * t212;
t194 = t277 * t268;
t167 = t193 * t267 + t194 * t211;
t257 = t167 * t209;
t256 = t186 * qJ(6);
t187 = t211 * t268 + t212 * t267;
t255 = t187 * t208;
t156 = pkin(3) * t186 - qJ(4) * t187 + t201;
t168 = -t193 * t211 + t194 * t267;
t138 = t156 * t209 - t168 * t208;
t135 = pkin(4) * t186 - t187 * t205 + t138;
t139 = t156 * t208 + t168 * t209;
t137 = -pkin(9) * t255 + t139;
t127 = t135 * t210 + t137 * t266;
t124 = t127 + t256;
t228 = -t135 * t266 + t137 * t210;
t125 = t228 - t265;
t253 = t124 * t274 + t125 * t185;
t252 = t153 * t185 + t154 * t274;
t251 = t163 * t185 + t164 * t274;
t198 = -pkin(3) - t264;
t155 = -pkin(5) * t274 - qJ(6) * t185 + t198;
t235 = t267 * pkin(2);
t148 = -t235 + t155;
t250 = t148 + t155;
t200 = -t235 - pkin(3);
t189 = t200 - t264;
t249 = t189 + t198;
t247 = t246 * qJ(4);
t150 = t274 * t187;
t245 = MDP(22) * t150;
t244 = t148 * MDP(32);
t149 = t185 * t187;
t243 = t149 * MDP(25);
t242 = t150 * MDP(24);
t241 = t155 * MDP(32);
t240 = t186 * MDP(26);
t239 = t200 * MDP(21);
t237 = 0.2e1 * t268;
t236 = MDP(28) - MDP(31);
t234 = (-pkin(5) * t185 + qJ(6) * t274) * MDP(30) + t278;
t232 = MDP(15) + (MDP(22) * t185 - t270 * t274) * t185;
t229 = -MDP(32) * pkin(5) - MDP(29);
t227 = t246 * MDP(21);
t226 = -MDP(27) + t229;
t225 = -pkin(3) * t187 - qJ(4) * t186;
t224 = MDP(32) * qJ(6) - t236;
t223 = -t138 * t208 + t139 * t209;
t222 = -t186 * t197 + t187 * t200;
t221 = 0.2e1 * t273;
t220 = MDP(18) * t208 + MDP(19) * t209;
t219 = -MDP(27) * t228 - t127 * MDP(28);
t218 = 0.2e1 * t275;
t145 = pkin(4) * t255 + t167;
t217 = 0.2e1 * t276;
t216 = -t273 + t275 + t276;
t215 = (MDP(16) * t267 - MDP(17) * t211) * pkin(2);
t214 = t223 * MDP(20) + (-t149 * t185 + t150 * t274) * MDP(23) + t185 * t245 - t167 * MDP(16) - t168 * MDP(17) + t187 * MDP(13) + (-MDP(14) + t278) * t186;
t176 = t185 * MDP(30);
t160 = t167 * t208;
t141 = t145 * t185;
t140 = t145 * t274;
t131 = pkin(5) * t149 - qJ(6) * t150 + t145;
t130 = t131 * t185;
t129 = t131 * t274;
t1 = [(t138 ^ 2 + t139 ^ 2 + t167 ^ 2) * MDP(21) + (t124 ^ 2 + t125 ^ 2 + t131 ^ 2) * MDP(32) + pkin(1) * MDP(9) * t237 + MDP(1) + (MDP(11) * t187 + MDP(17) * t272) * t187 + (t149 * t270 + t245) * t150 + (-0.2e1 * MDP(10) * pkin(1) + MDP(4) * t212 + MDP(5) * t237) * t212 + (-0.2e1 * MDP(12) * t187 + MDP(16) * t272 + t240 + 0.2e1 * t242 - 0.2e1 * t243) * t186 + (-t124 * t149 + t125 * t150) * t269 + 0.2e1 * (t149 * MDP(27) + t150 * MDP(28)) * t145 + 0.2e1 * (t149 * MDP(29) - t150 * MDP(31)) * t131 + 0.2e1 * (t138 * MDP(18) - t139 * MDP(19) - t125 * MDP(29) + t124 * MDP(31) + t219) * t186 + 0.2e1 * ((-t138 * t209 - t139 * t208) * MDP(20) + t220 * t167) * t187; (-t149 * t154 + t150 * t153 + t253) * MDP(30) + (-t148 * t150 - t130 + t260) * MDP(31) + (-MDP(10) * t268 - MDP(9) * t212) * pkin(7) + t214 + (t124 * t154 + t125 * t153 + t131 * t148) * MDP(32) + (t148 * t149 - t129 - t261) * MDP(29) + (t149 * t189 - t140 - t261) * MDP(27) + (t150 * t189 + t141 - t260) * MDP(28) + t212 * MDP(6) + (t208 * t222 - t257) * MDP(18) + (t167 * t200 + t197 * t223) * MDP(21) + t268 * MDP(7) + (t209 * t222 + t160) * MDP(19); MDP(8) + (t153 ^ 2 + t154 ^ 2) * MDP(32) + t197 ^ 2 * t227 + (-t221 + t239) * t200 + t189 * t217 + (t218 + t244) * t148 + 0.2e1 * t215 + t248 * t271 + t252 * t269 + t232; (-t149 * t164 + t150 * t163 + t253) * MDP(30) + (-t150 * t155 - t130 + t258) * MDP(31) + t214 + (t124 * t164 + t125 * t163 + t131 * t155) * MDP(32) + (t149 * t155 - t129 - t259) * MDP(29) + (t149 * t198 - t140 - t259) * MDP(27) + (t150 * t198 + t141 - t258) * MDP(28) + (t208 * t225 - t257) * MDP(18) + (-t167 * pkin(3) + qJ(4) * t223) * MDP(21) + (t209 * t225 + t160) * MDP(19); (t247 + t248) * MDP(20) + (-t200 * pkin(3) + qJ(4) * t248) * MDP(21) + (t251 + t252) * MDP(30) + (t148 * t155 + t153 * t163 + t154 * t164) * MDP(32) + t215 + (MDP(28) * t249 - MDP(31) * t250) * t185 - (MDP(27) * t249 + MDP(29) * t250) * t274 + t232 + t273 * (pkin(3) - t200); (t163 ^ 2 + t164 ^ 2) * MDP(32) + qJ(4) ^ 2 * t227 + t198 * t217 + (t218 + t241) * t155 + (t221 + t262) * pkin(3) + t247 * t271 + t251 * t269 + t232; t167 * MDP(21) + t131 * MDP(32) + t220 * t187 + t236 * t150 + (MDP(27) + MDP(29)) * t149; t216 + t239 + t244; t216 + t241 - t262; MDP(21) + MDP(32); t242 - t243 + t240 + (-t228 + 0.2e1 * t265) * MDP(29) + (-pkin(5) * t150 - qJ(6) * t149) * MDP(30) + (t127 + 0.2e1 * t256) * MDP(31) + (-pkin(5) * t125 + qJ(6) * t124) * MDP(32) + t219; t153 * t226 + t154 * t224 + t234; t163 * t226 + t164 * t224 + t234; 0; MDP(26) + 0.2e1 * pkin(5) * MDP(29) + 0.2e1 * qJ(6) * MDP(31) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(32); -MDP(29) * t186 + MDP(30) * t150 + MDP(32) * t125; MDP(32) * t153 + t176; MDP(32) * t163 + t176; 0; t229; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
