% Calculate joint inertia matrix for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR9_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:14:02
% EndTime: 2019-03-09 14:14:07
% DurationCPUTime: 1.48s
% Computational Cost: add. (2508->248), mult. (5605->359), div. (0->0), fcn. (6609->12), ass. (0->123)
t214 = sin(qJ(6));
t218 = cos(qJ(6));
t260 = t214 * MDP(31) + t218 * MDP(32);
t210 = sin(pkin(12));
t212 = cos(pkin(12));
t213 = cos(pkin(6));
t211 = sin(pkin(6));
t217 = sin(qJ(2));
t263 = t211 * t217;
t185 = t210 * t263 - t213 * t212;
t186 = t210 * t213 + t212 * t263;
t216 = sin(qJ(4));
t219 = cos(qJ(4));
t166 = t185 * t219 + t186 * t216;
t167 = -t185 * t216 + t186 * t219;
t215 = sin(qJ(5));
t274 = cos(qJ(5));
t153 = t274 * t166 + t167 * t215;
t154 = -t215 * t166 + t274 * t167;
t282 = t154 * MDP(24) - t153 * MDP(25);
t231 = MDP(34) * t218 - MDP(35) * t214;
t270 = pkin(9) + qJ(3);
t193 = t270 * t210;
t194 = t270 * t212;
t177 = -t193 * t216 + t194 * t219;
t191 = t210 * t216 - t212 * t219;
t164 = -pkin(10) * t191 + t177;
t176 = -t219 * t193 - t194 * t216;
t192 = t210 * t219 + t212 * t216;
t229 = -pkin(10) * t192 + t176;
t148 = t164 * t215 - t274 * t229;
t149 = t274 * t164 + t215 * t229;
t174 = t274 * t191 + t192 * t215;
t175 = -t215 * t191 + t274 * t192;
t281 = -t175 * MDP(24) + t174 * MDP(25) + t148 * MDP(27) + t149 * MDP(28);
t280 = 2 * MDP(13);
t279 = -2 * MDP(16);
t278 = 2 * MDP(21);
t277 = 0.2e1 * MDP(27);
t276 = 0.2e1 * MDP(34);
t275 = 0.2e1 * MDP(35);
t273 = pkin(1) * t217;
t220 = cos(qJ(2));
t272 = pkin(1) * t220;
t269 = pkin(2) * MDP(14);
t262 = t211 * t220;
t241 = pkin(8) * t262;
t182 = t241 + (qJ(3) + t273) * t213;
t183 = (-pkin(2) * t220 - qJ(3) * t217 - pkin(1)) * t211;
t160 = -t182 * t210 + t212 * t183;
t158 = -pkin(3) * t262 - pkin(9) * t186 + t160;
t161 = t212 * t182 + t210 * t183;
t159 = -pkin(9) * t185 + t161;
t138 = t219 * t158 - t159 * t216;
t242 = pkin(4) * t262;
t134 = -pkin(10) * t167 + t138 - t242;
t139 = t158 * t216 + t159 * t219;
t137 = -pkin(10) * t166 + t139;
t129 = t274 * t134 - t215 * t137;
t127 = pkin(5) * t262 - t129;
t268 = t127 * t218;
t267 = t148 * t218;
t266 = t153 * t214;
t265 = t153 * t218;
t264 = t175 * t214;
t261 = t213 * t217;
t258 = MDP(23) * t154;
t257 = MDP(28) * t154;
t256 = MDP(28) * t175;
t142 = t154 * t218 - t214 * t262;
t255 = MDP(31) * t142;
t141 = t154 * t214 + t218 * t262;
t254 = MDP(32) * t141;
t253 = MDP(33) * t153;
t250 = t142 * MDP(29);
t249 = t174 * MDP(33);
t196 = pkin(8) * t263;
t184 = t196 + (-pkin(2) - t272) * t213;
t248 = t184 * MDP(14);
t247 = t191 * MDP(20);
t246 = t192 * MDP(15);
t245 = t210 * MDP(12);
t244 = t212 * MDP(11);
t243 = t218 * MDP(29);
t240 = t274 * pkin(4);
t198 = -pkin(3) * t212 - pkin(2);
t239 = t274 * t137;
t238 = t214 * t218 * MDP(30);
t207 = t214 ^ 2;
t237 = t207 * MDP(29) + MDP(26) + 0.2e1 * t238;
t236 = -pkin(5) * t175 - pkin(11) * t174;
t235 = -t160 * t210 + t161 * t212;
t199 = pkin(4) * t215 + pkin(11);
t200 = -t240 - pkin(5);
t234 = -t174 * t199 + t175 * t200;
t181 = pkin(4) * t191 + t198;
t233 = t167 * MDP(17) - t166 * MDP(18);
t232 = t218 * MDP(31) - t214 * MDP(32);
t230 = t214 * MDP(34) + t218 * MDP(35);
t130 = t215 * t134 + t239;
t228 = -MDP(23) + t232;
t227 = MDP(27) + t231;
t226 = (t274 * MDP(27) - t215 * MDP(28)) * pkin(4);
t208 = t218 ^ 2;
t225 = (-t207 + t208) * t175 * MDP(30) + t243 * t264 - t281 + t260 * t174;
t170 = pkin(3) * t185 + t184;
t224 = -MDP(26) * t262 + (-t141 * t214 + t142 * t218) * MDP(30) + t214 * t250 + t282 + t260 * t153;
t155 = pkin(4) * t166 + t170;
t223 = -t192 * MDP(17) + t191 * MDP(18) - t176 * MDP(20) + t177 * MDP(21);
t128 = -pkin(11) * t262 + t130;
t131 = pkin(5) * t153 - pkin(11) * t154 + t155;
t124 = -t128 * t214 + t131 * t218;
t125 = t128 * t218 + t131 * t214;
t222 = MDP(34) * t124 - MDP(35) * t125 + t253 - t254 + t255;
t205 = t211 ^ 2;
t189 = pkin(1) * t261 + t241;
t188 = t213 * t272 - t196;
t150 = pkin(5) * t174 - pkin(11) * t175 + t181;
t143 = t148 * t214;
t136 = t149 * t218 + t150 * t214;
t135 = -t149 * t214 + t150 * t218;
t126 = t127 * t214;
t1 = [MDP(1) + t154 ^ 2 * MDP(22) + (t160 ^ 2 + t161 ^ 2 + t184 ^ 2) * MDP(14) + t213 ^ 2 * MDP(8) + (t167 * MDP(15) + t166 * t279) * t167 + (-0.2e1 * t141 * MDP(30) + t250) * t142 + ((MDP(4) * t217 + 0.2e1 * MDP(5) * t220) * t217 + (MDP(19) + MDP(26)) * t220 ^ 2) * t205 + (t253 - 0.2e1 * t254 + 0.2e1 * t255 - 0.2e1 * t258) * t153 + 0.2e1 * (MDP(6) * t261 + (MDP(7) * t213 - t233 - t282) * t220) * t211 + (-t129 * t262 + t153 * t155) * t277 + 0.2e1 * (t188 * t213 + t205 * t272) * MDP(9) + 0.2e1 * (-t160 * t262 + t184 * t185) * MDP(11) + 0.2e1 * (t161 * t262 + t184 * t186) * MDP(12) + 0.2e1 * (-t138 * t262 + t166 * t170) * MDP(20) + (t139 * t262 + t167 * t170) * t278 + 0.2e1 * (t130 * t262 + t154 * t155) * MDP(28) + 0.2e1 * (-t189 * t213 - t205 * t273) * MDP(10) + (-t160 * t186 - t161 * t185) * t280 + (-t125 * t153 + t127 * t142) * t275 + (t124 * t153 + t127 * t141) * t276; t167 * t246 + t188 * MDP(9) - t189 * MDP(10) + (-t192 * t166 - t167 * t191) * MDP(16) - pkin(2) * t248 + t235 * MDP(13) + t213 * MDP(8) + (t135 * t153 + t148 * t141) * MDP(34) + (-t136 * t153 + t148 * t142) * MDP(35) + (t198 * t166 + t170 * t191) * MDP(20) + (t198 * t167 + t170 * t192) * MDP(21) + (-pkin(2) * t185 - t184 * t212) * MDP(11) + (-pkin(2) * t186 + t184 * t210) * MDP(12) + (MDP(27) * t153 + t257) * t181 + (t235 * MDP(14) + (-t185 * t212 + t186 * t210) * MDP(13)) * qJ(3) + (MDP(27) * t155 + t222 - t258) * t174 + (t142 * t243 + t154 * MDP(22) + (-t141 * t218 - t142 * t214) * MDP(30) + t155 * MDP(28) + t230 * t127 + t228 * t153) * t175 + (MDP(6) * t217 + (MDP(7) + (MDP(11) * t210 + MDP(12) * t212) * qJ(3) + t223 + t281) * t220) * t211; 0.2e1 * t198 * t247 + 0.2e1 * t181 * t256 + MDP(8) + (0.2e1 * t244 - 0.2e1 * t245 + t269) * pkin(2) + (t191 * t279 + t198 * t278 + t246) * t192 + (MDP(29) * t208 + MDP(22) - 0.2e1 * t238) * t175 ^ 2 + (0.2e1 * t228 * t175 + t181 * t277 + t249) * t174 + (t135 * t174 + t148 * t264) * t276 + (-t136 * t174 + t175 * t267) * t275 + (MDP(14) * qJ(3) + t280) * (t210 ^ 2 + t212 ^ 2) * qJ(3); t185 * MDP(11) + t186 * MDP(12) + t166 * MDP(20) + t167 * MDP(21) + t227 * t153 + t248 + t257; t192 * MDP(21) + t227 * t174 - t244 + t245 + t247 + t256 - t269; MDP(14); -MDP(19) * t262 + t138 * MDP(20) - t139 * MDP(21) + (-t240 * t262 + t129) * MDP(27) + (-t239 + (-t134 + t242) * t215) * MDP(28) + (t141 * t200 - t199 * t266 - t268) * MDP(34) + (t142 * t200 - t199 * t265 + t126) * MDP(35) + t224 + t233; (t234 * t214 - t267) * MDP(34) + (t234 * t218 + t143) * MDP(35) - t223 + t225; 0; -0.2e1 * t200 * t231 + MDP(19) + 0.2e1 * t226 + t237; t129 * MDP(27) - t130 * MDP(28) + (-pkin(5) * t141 - pkin(11) * t266 - t268) * MDP(34) + (-pkin(5) * t142 - pkin(11) * t265 + t126) * MDP(35) + t224; (t236 * t214 - t267) * MDP(34) + (t236 * t218 + t143) * MDP(35) + t225; 0; t226 + t237 + t231 * (pkin(5) - t200); 0.2e1 * pkin(5) * t231 + t237; t222; t135 * MDP(34) - t136 * MDP(35) + t232 * t175 + t249; t231; -t230 * t199 + t260; -t230 * pkin(11) + t260; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
