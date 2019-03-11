% Calculate joint inertia matrix for
% S6RRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP3_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:09:32
% EndTime: 2019-03-10 01:09:34
% DurationCPUTime: 0.85s
% Computational Cost: add. (1443->203), mult. (2655->279), div. (0->0), fcn. (2948->8), ass. (0->100)
t213 = sin(qJ(4));
t217 = cos(qJ(4));
t212 = sin(qJ(5));
t216 = cos(qJ(5));
t189 = t212 * t217 + t213 * t216;
t264 = -t212 * t213 + t216 * t217;
t241 = t189 * MDP(27) + MDP(28) * t264;
t270 = t213 * MDP(20) + t217 * MDP(21) + t241;
t268 = pkin(8) + pkin(7);
t214 = sin(qJ(3));
t201 = pkin(2) * t214 + pkin(9);
t183 = (-pkin(10) - t201) * t213;
t209 = t217 * pkin(10);
t184 = t201 * t217 + t209;
t151 = t216 * t183 - t184 * t212;
t152 = t183 * t212 + t184 * t216;
t267 = t151 * MDP(30) - t152 * MDP(31);
t194 = (-pkin(9) - pkin(10)) * t213;
t196 = pkin(9) * t217 + t209;
t163 = t216 * t194 - t196 * t212;
t165 = t194 * t212 + t196 * t216;
t266 = t163 * MDP(30) - t165 * MDP(31);
t265 = -MDP(30) * t264 + MDP(31) * t189;
t223 = MDP(23) * t217 - MDP(24) * t213;
t260 = cos(qJ(2));
t205 = -pkin(2) * t260 - pkin(1);
t263 = 0.2e1 * t205;
t262 = -2 * MDP(26);
t261 = 2 * MDP(32);
t218 = cos(qJ(3));
t259 = pkin(2) * t218;
t215 = sin(qJ(2));
t188 = t214 * t215 - t218 * t260;
t258 = pkin(4) * t188;
t257 = pkin(4) * t212;
t256 = pkin(4) * t216;
t254 = MDP(33) * pkin(5);
t253 = qJ(6) * t189;
t190 = t214 * t260 + t218 * t215;
t154 = t188 * pkin(3) - t190 * pkin(9) + t205;
t195 = t268 * t215;
t197 = t268 * t260;
t166 = -t214 * t195 + t197 * t218;
t250 = t166 * t217;
t127 = t250 + (-pkin(10) * t190 + t154) * t213;
t252 = t127 * t216;
t164 = t195 * t218 + t214 * t197;
t251 = t164 * t217;
t249 = t190 * t213;
t248 = t190 * t217;
t246 = t213 * t217;
t129 = t217 * t154 - t166 * t213;
t126 = -pkin(10) * t248 + t129 + t258;
t122 = t216 * t126 - t127 * t212;
t147 = t264 * t190;
t119 = pkin(5) * t188 - qJ(6) * t147 + t122;
t123 = t126 * t212 + t252;
t146 = t189 * t190;
t121 = -qJ(6) * t146 + t123;
t244 = -t119 * t189 + t121 * t264;
t138 = t151 - t253;
t182 = t264 * qJ(6);
t139 = t152 + t182;
t243 = -t138 * t189 + t139 * t264;
t141 = t163 - t253;
t142 = t165 + t182;
t242 = -t141 * t189 + t142 * t264;
t237 = MDP(25) * t147;
t143 = t146 * MDP(28);
t144 = t147 * MDP(27);
t234 = t216 * MDP(30);
t233 = 0.2e1 * t260;
t232 = MDP(22) + MDP(29);
t231 = t188 * MDP(29) - t143 + t144;
t230 = -pkin(5) * t189 * MDP(32) + t241;
t204 = -t217 * pkin(4) - pkin(3);
t229 = MDP(19) * t246;
t228 = -pkin(3) * t190 - pkin(9) * t188;
t140 = pkin(4) * t249 + t164;
t202 = pkin(5) + t256;
t227 = (-t189 * t202 + t257 * t264) * MDP(32) + t270;
t210 = t213 ^ 2;
t226 = t210 * MDP(18) + MDP(15) + 0.2e1 * t229 + (MDP(25) * t189 - t262 * t264) * t189;
t203 = -pkin(3) - t259;
t225 = -t188 * t201 + t190 * t203;
t168 = -pkin(5) * t264 + t204;
t224 = MDP(20) * t217 - MDP(21) * t213;
t222 = -MDP(23) * t213 - MDP(24) * t217;
t221 = (MDP(16) * t218 - MDP(17) * t214) * pkin(2);
t220 = 0.2e1 * t265;
t211 = t217 ^ 2;
t219 = (-t146 * t189 + t147 * t264) * MDP(26) + t189 * t237 - t164 * MDP(16) - t166 * MDP(17) + ((-t210 + t211) * MDP(19) + MDP(18) * t246 + MDP(13)) * t190 + (-MDP(14) + t270) * t188;
t193 = t204 - t259;
t167 = t168 - t259;
t155 = t164 * t213;
t134 = t140 * t189;
t133 = t140 * t264;
t130 = t154 * t213 + t250;
t128 = t146 * pkin(5) + t140;
t1 = [(t119 ^ 2 + t121 ^ 2 + t128 ^ 2) * MDP(33) + pkin(1) * MDP(9) * t233 + MDP(1) + t232 * t188 ^ 2 + (t146 * t262 + t237) * t147 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t215 + MDP(5) * t233) * t215 + (MDP(16) * t263 - 0.2e1 * t143 + 0.2e1 * t144) * t188 + (-t119 * t147 - t121 * t146) * t261 + 0.2e1 * (t122 * t188 + t140 * t146) * MDP(30) + 0.2e1 * (-t123 * t188 + t140 * t147) * MDP(31) + 0.2e1 * (t129 * t188 + t164 * t249) * MDP(23) + 0.2e1 * (-t130 * t188 + t164 * t248) * MDP(24) + (MDP(17) * t263 + 0.2e1 * (-MDP(12) + t224) * t188 + (t211 * MDP(18) + MDP(11) - 0.2e1 * t229) * t190) * t190; (-t138 * t147 - t139 * t146 + t244) * MDP(32) + t215 * MDP(6) + (t146 * t193 + t151 * t188 - t133) * MDP(30) + (t147 * t193 - t152 * t188 + t134) * MDP(31) + (t119 * t138 + t121 * t139 + t128 * t167) * MDP(33) + (-MDP(10) * t260 - t215 * MDP(9)) * pkin(7) + (t213 * t225 - t251) * MDP(23) + (t217 * t225 + t155) * MDP(24) + t219 + t260 * MDP(7); MDP(8) + t243 * t261 + (t138 ^ 2 + t139 ^ 2 + t167 ^ 2) * MDP(33) + t193 * t220 + t226 - 0.2e1 * t203 * t223 + 0.2e1 * t221; (-t141 * t147 - t142 * t146 + t244) * MDP(32) + (t146 * t204 + t163 * t188 - t133) * MDP(30) + (t147 * t204 - t165 * t188 + t134) * MDP(31) + (t119 * t141 + t121 * t142 + t128 * t168) * MDP(33) + (t213 * t228 - t251) * MDP(23) + (t217 * t228 + t155) * MDP(24) + t219; (t242 + t243) * MDP(32) + (t138 * t141 + t139 * t142 + t167 * t168) * MDP(33) + t221 + t226 + t223 * (pkin(3) - t203) + t265 * (t193 + t204); t242 * t261 + (t141 ^ 2 + t142 ^ 2 + t168 ^ 2) * MDP(33) + t204 * t220 + 0.2e1 * t223 * pkin(3) + t226; t188 * MDP(22) + t129 * MDP(23) - t130 * MDP(24) + (t188 * t256 + t122) * MDP(30) + (-t252 + (-t126 - t258) * t212) * MDP(31) + (-t146 * t257 - t147 * t202) * MDP(32) + (t119 * t202 + t121 * t257) * MDP(33) + t224 * t190 + t231; (t138 * t202 + t139 * t257) * MDP(33) + t222 * t201 + t227 + t267; (t141 * t202 + t142 * t257) * MDP(33) + t222 * pkin(9) + t227 + t266; t202 ^ 2 * MDP(33) + (0.2e1 * t234 + (MDP(33) * t257 - 0.2e1 * MDP(31)) * t212) * pkin(4) + t232; t122 * MDP(30) - t123 * MDP(31) + (-MDP(32) * t147 + MDP(33) * t119) * pkin(5) + t231; t138 * t254 + t230 + t267; t141 * t254 + t230 + t266; t202 * t254 + MDP(29) + (-MDP(31) * t212 + t234) * pkin(4); MDP(33) * pkin(5) ^ 2 + MDP(29); t128 * MDP(33); t167 * MDP(33); t168 * MDP(33); 0; 0; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
