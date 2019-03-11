% Calculate joint inertia matrix for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRRPP1_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:46:55
% EndTime: 2019-03-09 20:46:58
% DurationCPUTime: 1.02s
% Computational Cost: add. (1898->224), mult. (3392->301), div. (0->0), fcn. (3728->8), ass. (0->94)
t275 = 2 * MDP(28);
t277 = 2 * MDP(25);
t284 = t277 + t275;
t221 = sin(qJ(4));
t224 = cos(qJ(4));
t283 = t221 * MDP(20) + t224 * MDP(21);
t216 = t224 * qJ(5);
t222 = sin(qJ(3));
t242 = pkin(2) * t222 + pkin(9);
t234 = t224 * t242;
t192 = t234 + t216;
t219 = sin(pkin(10));
t220 = cos(pkin(10));
t229 = (-qJ(5) - t242) * t221;
t165 = t192 * t219 - t220 * t229;
t167 = t220 * t192 + t219 * t229;
t194 = t219 * t224 + t220 * t221;
t279 = -t219 * t221 + t220 * t224;
t281 = t165 * t194 + t167 * t279;
t201 = t224 * pkin(9) + t216;
t239 = (-qJ(5) - pkin(9)) * t221;
t175 = t201 * t219 - t220 * t239;
t177 = t220 * t201 + t219 * t239;
t280 = t175 * t194 + t177 * t279;
t254 = -MDP(27) * t279 - t194 * MDP(29);
t231 = t224 * MDP(23) - t221 * MDP(24);
t274 = cos(qJ(2));
t212 = -t274 * pkin(2) - pkin(1);
t278 = 0.2e1 * t212;
t276 = 0.2e1 * MDP(27);
t273 = cos(qJ(3));
t272 = t224 * pkin(4);
t270 = MDP(26) * pkin(4);
t223 = sin(qJ(2));
t202 = (-pkin(8) - pkin(7)) * t223;
t244 = t274 * pkin(7);
t203 = t274 * pkin(8) + t244;
t178 = -t273 * t202 + t222 * t203;
t267 = t178 * t224;
t179 = t222 * t202 + t273 * t203;
t266 = t179 * t224;
t198 = t222 * t274 + t273 * t223;
t265 = t198 * t221;
t262 = t221 * t224;
t197 = t222 * t223 - t273 * t274;
t170 = t197 * pkin(3) - t198 * pkin(9) + t212;
t147 = t224 * t170 - t179 * t221;
t143 = pkin(4) * t197 - t198 * t216 + t147;
t145 = t266 + (-qJ(5) * t198 + t170) * t221;
t137 = t219 * t143 + t220 * t145;
t132 = t197 * qJ(6) + t137;
t136 = t220 * t143 - t145 * t219;
t133 = -pkin(5) * t197 - t136;
t261 = t132 * t279 + t133 * t194;
t260 = -t136 * t194 + t137 * t279;
t211 = -pkin(3) - t272;
t168 = -pkin(5) * t279 - t194 * qJ(6) + t211;
t253 = MDP(30) * t168;
t243 = t273 * pkin(2);
t159 = -t243 + t168;
t252 = t159 * MDP(30);
t205 = pkin(4) * t219 + qJ(6);
t251 = t205 * MDP(29);
t250 = t221 * MDP(23);
t247 = 0.2e1 * t274;
t246 = t165 ^ 2 + t167 ^ 2;
t245 = t175 ^ 2 + t177 ^ 2;
t241 = MDP(19) * t262;
t217 = t221 ^ 2;
t240 = t217 * MDP(18) + MDP(15) + 0.2e1 * t241;
t160 = t194 * t198;
t161 = t279 * t198;
t238 = -t167 * t160 + t161 * t165;
t237 = -t177 * t160 + t161 * t175;
t236 = t165 * t175 + t167 * t177;
t207 = pkin(4) * t220 + pkin(5);
t235 = (-t194 * t207 + t205 * t279) * MDP(28) + (-t194 * t220 + t219 * t279) * pkin(4) * MDP(25) + t283;
t210 = -t243 - pkin(3);
t233 = -pkin(3) * t198 - pkin(9) * t197;
t232 = t224 * MDP(20) - t221 * MDP(21);
t230 = 0.2e1 * t254;
t155 = pkin(4) * t265 + t178;
t228 = -t242 * t197 + t210 * t198;
t227 = (t273 * MDP(16) - t222 * MDP(17)) * pkin(2);
t218 = t224 ^ 2;
t226 = -t178 * MDP(16) - t179 * MDP(17) + ((-t217 + t218) * MDP(19) + MDP(18) * t262 + MDP(13)) * t198 + (-MDP(14) + t283) * t197;
t200 = t210 - t272;
t187 = t194 * MDP(28);
t171 = t178 * t221;
t148 = t170 * t221 + t266;
t140 = t160 * pkin(5) - t161 * qJ(6) + t155;
t139 = t140 * t194;
t138 = t140 * t279;
t1 = [MDP(1) + pkin(1) * MDP(9) * t247 + t198 * MDP(17) * t278 + (t136 ^ 2 + t137 ^ 2 + t155 ^ 2) * MDP(26) + (t132 ^ 2 + t133 ^ 2 + t140 ^ 2) * MDP(30) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t223 + MDP(5) * t247) * t223 + (t218 * MDP(18) + MDP(11) - 0.2e1 * t241) * t198 ^ 2 + (MDP(16) * t278 + t197 * MDP(22) + 0.2e1 * (-MDP(12) + t232) * t198) * t197 + 0.2e1 * (t147 * t197 + t178 * t265) * MDP(23) + 0.2e1 * (-t148 * t197 + t198 * t267) * MDP(24) + (-t136 * t161 - t137 * t160) * t277 + (-t133 * t197 + t140 * t160) * t276 + (-t132 * t160 + t133 * t161) * t275 + 0.2e1 * (t132 * t197 - t140 * t161) * MDP(29); (t159 * t160 - t165 * t197 - t138) * MDP(27) + (-t136 * t165 + t137 * t167 + t155 * t200) * MDP(26) + (t228 * t221 - t267) * MDP(23) + (t228 * t224 + t171) * MDP(24) + (t132 * t167 + t133 * t165 + t140 * t159) * MDP(30) + (-t159 * t161 + t167 * t197 - t139) * MDP(29) + (t238 + t260) * MDP(25) + (t238 + t261) * MDP(28) + t226 + t274 * MDP(7) - MDP(10) * t244 + (-pkin(7) * MDP(9) + MDP(6)) * t223; MDP(8) + (t200 ^ 2 + t246) * MDP(26) + t246 * MDP(30) + (t230 + t252) * t159 + t240 - 0.2e1 * t231 * t210 + 0.2e1 * t227 + t284 * t281; (t132 * t177 + t133 * t175 + t140 * t168) * MDP(30) + (t160 * t168 - t175 * t197 - t138) * MDP(27) + (-t161 * t168 + t177 * t197 - t139) * MDP(29) + (t237 + t260) * MDP(25) + (t237 + t261) * MDP(28) + (-t136 * t175 + t137 * t177 + t155 * t211) * MDP(26) + t226 + (t233 * t221 - t267) * MDP(23) + (t233 * t224 + t171) * MDP(24); (t200 * t211 + t236) * MDP(26) + (t159 * t168 + t236) * MDP(30) + t227 + t240 + t231 * (pkin(3) - t210) + t254 * (t159 + t168) + (MDP(25) + MDP(28)) * (t280 + t281); (t211 ^ 2 + t245) * MDP(26) + t245 * MDP(30) + (t230 + t253) * t168 + 0.2e1 * t231 * pkin(3) + t240 + t284 * t280; t147 * MDP(23) - t148 * MDP(24) + t136 * MDP(27) + (-t160 * t205 - t161 * t207) * MDP(28) + t132 * MDP(29) + (t132 * t205 - t133 * t207) * MDP(30) + t232 * t198 + (MDP(22) + (pkin(5) + t207) * MDP(27) + t251) * t197 + ((-t160 * t219 - t161 * t220) * MDP(25) + (t136 * t220 + t137 * t219) * MDP(26)) * pkin(4); -t242 * t250 - MDP(24) * t234 - t165 * MDP(27) + t167 * MDP(29) + (-t165 * t207 + t167 * t205) * MDP(30) + (-t165 * t220 + t167 * t219) * t270 + t235; -t175 * MDP(27) + t177 * MDP(29) + (-t175 * t207 + t177 * t205) * MDP(30) + (-t224 * MDP(24) - t250) * pkin(9) + (-t175 * t220 + t177 * t219) * t270 + t235; MDP(22) + (t205 ^ 2 + t207 ^ 2) * MDP(30) + (t219 ^ 2 + t220 ^ 2) * MDP(26) * pkin(4) ^ 2 + t207 * t276 + 0.2e1 * t251; MDP(26) * t155 + t160 * MDP(27) - t161 * MDP(29) + t140 * MDP(30); MDP(26) * t200 + t252 + t254; MDP(26) * t211 + t253 + t254; 0; MDP(26) + MDP(30); -t197 * MDP(27) + t161 * MDP(28) + t133 * MDP(30); t165 * MDP(30) + t187; t175 * MDP(30) + t187; -t207 * MDP(30) - MDP(27); 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
