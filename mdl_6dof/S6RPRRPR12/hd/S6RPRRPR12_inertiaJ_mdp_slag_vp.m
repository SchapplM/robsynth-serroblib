% Calculate joint inertia matrix for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRPR12_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:53:01
% EndTime: 2019-03-09 05:53:04
% DurationCPUTime: 1.15s
% Computational Cost: add. (2171->240), mult. (5745->359), div. (0->0), fcn. (6501->12), ass. (0->113)
t191 = sin(pkin(6));
t190 = sin(pkin(7));
t248 = cos(pkin(6));
t218 = t248 * t190;
t192 = cos(pkin(7));
t247 = cos(pkin(12));
t219 = t192 * t247;
t263 = t191 * t219 + t218;
t226 = MDP(15) + MDP(30);
t225 = MDP(20) - MDP(23);
t224 = -MDP(21) + MDP(24);
t220 = t191 * t247;
t189 = sin(pkin(12));
t252 = pkin(1) * t189;
t169 = qJ(2) * t220 + t248 * t252;
t153 = t263 * pkin(9) + t169;
t195 = sin(qJ(3));
t198 = cos(qJ(3));
t223 = pkin(1) * t247;
t181 = t248 * t223;
t245 = t189 * t191;
t157 = t248 * pkin(2) + t181 + (-pkin(9) * t192 - qJ(2)) * t245;
t163 = (-pkin(9) * t189 * t190 - t247 * pkin(2) - pkin(1)) * t191;
t213 = t157 * t192 + t163 * t190;
t143 = -t195 * t153 + t213 * t198;
t260 = 2 * MDP(20);
t259 = 2 * MDP(22);
t258 = 2 * MDP(23);
t257 = 2 * MDP(24);
t256 = 2 * MDP(31);
t255 = 2 * MDP(32);
t254 = pkin(4) + pkin(11);
t253 = pkin(5) + pkin(10);
t155 = t195 * t245 - t263 * t198;
t251 = pkin(4) * t155;
t250 = MDP(25) * pkin(4);
t249 = MDP(25) * pkin(10);
t156 = t195 * t218 + (t189 * t198 + t195 * t219) * t191;
t166 = t190 * t220 - t248 * t192;
t194 = sin(qJ(4));
t197 = cos(qJ(4));
t150 = t156 * t194 + t166 * t197;
t193 = sin(qJ(6));
t196 = cos(qJ(6));
t146 = t150 * t193 + t155 * t196;
t246 = t146 * t193;
t154 = t155 * qJ(5);
t244 = t190 * t195;
t243 = t190 * t198;
t242 = t193 * t197;
t241 = t196 * t197;
t147 = -t157 * t190 + t192 * t163;
t138 = pkin(3) * t155 - pkin(10) * t156 + t147;
t144 = t153 * t198 + t213 * t195;
t142 = -pkin(10) * t166 + t144;
t135 = t194 * t138 + t197 * t142;
t240 = MDP(25) * qJ(5);
t221 = -qJ(5) * t194 - pkin(3);
t175 = -pkin(4) * t197 + t221;
t239 = MDP(25) * t175;
t238 = MDP(25) * pkin(10) ^ 2;
t237 = MDP(26) * t196;
t236 = MDP(31) * t193;
t235 = MDP(32) * t196;
t134 = t138 * t197 - t194 * t142;
t133 = -t134 - t251;
t234 = t133 * MDP(25);
t145 = -t150 * t196 + t155 * t193;
t233 = t145 * MDP(29);
t232 = t146 * MDP(28);
t231 = t150 * MDP(16);
t230 = t155 * MDP(11);
t229 = t155 * MDP(17);
t228 = t155 * MDP(18);
t227 = t156 * MDP(10);
t132 = -t154 - t135;
t222 = t196 * t193 * MDP(27);
t217 = MDP(23) - t250;
t215 = -MDP(20) + t217;
t214 = t224 + t240;
t211 = -MDP(28) * t193 - MDP(29) * t196;
t210 = MDP(31) * t196 - t193 * MDP(32);
t209 = t235 + t236;
t206 = MDP(16) + t211;
t205 = MDP(22) + t210;
t141 = pkin(3) * t166 - t143;
t151 = t156 * t197 - t166 * t194;
t129 = pkin(5) * t151 - t254 * t155 - t134;
t202 = -qJ(5) * t151 + t141;
t131 = t254 * t150 + t202;
t127 = t129 * t196 - t131 * t193;
t128 = t129 * t193 + t131 * t196;
t204 = t127 * MDP(31) - t128 * MDP(32) + t232 - t233;
t203 = (-MDP(31) * t254 + MDP(28)) * t196 + (MDP(32) * t254 - MDP(29)) * t193;
t201 = -pkin(4) * MDP(22) + MDP(17) + t203;
t188 = t197 ^ 2;
t187 = t196 ^ 2;
t186 = t194 ^ 2;
t185 = t193 ^ 2;
t184 = t191 ^ 2;
t178 = t253 * t197;
t177 = t253 * t194;
t173 = -t254 * t197 + t221;
t171 = t192 * t194 + t197 * t244;
t170 = -t197 * t192 + t194 * t244;
t168 = -qJ(2) * t245 + t181;
t161 = -t193 * t170 + t196 * t243;
t160 = t170 * t196 + t193 * t243;
t159 = t173 * t196 + t177 * t193;
t158 = -t173 * t193 + t177 * t196;
t136 = pkin(4) * t150 + t202;
t130 = -pkin(5) * t150 - t132;
t1 = [-0.2e1 * t166 * t227 - 0.2e1 * t150 * t228 + 0.2e1 * t166 * t230 + t166 ^ 2 * MDP(12) - 0.2e1 * t156 * t155 * MDP(9) + 0.2e1 * (-t168 * t189 + t247 * t169) * MDP(6) * t191 + 0.2e1 * (-t143 * t166 + t147 * t155) * MDP(13) + 0.2e1 * (-t135 * t155 + t141 * t151) * MDP(21) + 0.2e1 * (t144 * t166 + t147 * t156) * MDP(14) + 0.2e1 * (t168 * t248 + t184 * t223) * MDP(4) + 0.2e1 * (-t169 * t248 - t184 * t252) * MDP(5) + t226 * t151 ^ 2 + 0.2e1 * (t232 + t229) * t151 - 0.2e1 * (t233 + t231) * t151 + (MDP(26) * t146 - 0.2e1 * t145 * MDP(27) + t130 * t255) * t146 + (pkin(1) ^ 2 * t184 + t168 ^ 2 + t169 ^ 2) * MDP(7) + t156 ^ 2 * MDP(8) + t155 ^ 2 * MDP(19) - t128 * t151 * t255 + (t132 ^ 2 + t133 ^ 2 + t136 ^ 2) * MDP(25) + MDP(1) + (t127 * t151 + t130 * t145) * t256 + (-t132 * t155 - t136 * t151) * t257 + (t133 * t155 - t136 * t150) * t258 + (t132 * t150 + t133 * t151) * t259 + (t134 * t155 + t141 * t150) * t260; (t155 * t192 - t166 * t243) * MDP(13) + (t156 * t192 + t166 * t244) * MDP(14) + (-t150 * t171 + t151 * t170) * MDP(22) + (-t132 * t171 + t133 * t170 - t136 * t243) * MDP(25) + (t145 * t171 + t151 * t160) * MDP(31) + (t146 * t171 + t151 * t161) * MDP(32) + (-t247 * MDP(4) + t189 * MDP(5) - pkin(1) * MDP(7)) * t191 - t225 * (t150 * t243 + t155 * t170) + t224 * (t151 * t243 + t155 * t171); MDP(7) + (t190 ^ 2 * t198 ^ 2 + t170 ^ 2 + t171 ^ 2) * MDP(25); t227 - t230 - t166 * MDP(12) + t143 * MDP(13) - t144 * MDP(14) + (t145 * t178 + t151 * t158) * MDP(31) + (t146 * t178 - t151 * t159) * MDP(32) + (-t150 * MDP(20) - t151 * MDP(21)) * pkin(3) + (-t150 * MDP(23) - t151 * MDP(24) + MDP(25) * t136) * t175 + (-t231 + t229 + t141 * MDP(21) + t133 * MDP(22) - t136 * MDP(24) + t226 * t151 + (t151 * MDP(22) - t225 * t155 + t234) * pkin(10) + t204) * t194 + (t228 - t141 * MDP(20) - t132 * MDP(22) + t136 * MDP(23) - MDP(26) * t246 + (t145 * t193 - t146 * t196) * MDP(27) + t210 * t130 + t206 * t151 + (-t150 * MDP(22) - t132 * MDP(25) + t224 * t155) * pkin(10)) * t197; (t160 * t194 + t171 * t241) * MDP(31) + (t161 * t194 - t171 * t242) * MDP(32) + (-t195 * MDP(14) + (t224 * t194 + t225 * t197 + MDP(13) - t239) * t198) * t190 + (MDP(22) + t249) * (t170 * t194 + t171 * t197); pkin(3) * t197 * t260 + MDP(12) + (t197 * t258 + t239) * t175 + (MDP(26) * t185 + 0.2e1 * t222 + t238) * t188 + (t226 + t238) * t186 + 0.2e1 * (-pkin(3) * MDP(21) - t175 * MDP(24) + t206 * t197) * t194 + (t158 * t194 + t178 * t241) * t256 + (-t159 * t194 - t178 * t242) * t255 + (t186 + t188) * pkin(10) * t259; t155 * MDP(19) + t134 * MDP(20) - t135 * MDP(21) + (-t134 - 0.2e1 * t251) * MDP(23) + (-t132 + t154) * MDP(24) + (-pkin(4) * t133 - qJ(5) * t132) * MDP(25) + t146 * t237 + (-t145 * t196 - t246) * MDP(27) + (qJ(5) * t145 + t130 * t193) * MDP(31) + (qJ(5) * t146 + t130 * t196) * MDP(32) + (-qJ(5) * MDP(22) - MDP(18)) * t150 + t201 * t151; t215 * t170 + (t209 + t214) * t171; t209 * t178 + t201 * t194 + (MDP(18) - t193 * t237 + (t185 - t187) * MDP(27) + t205 * qJ(5)) * t197 + (t215 * t194 + t214 * t197) * pkin(10); -0.2e1 * t222 + t187 * MDP(26) + MDP(19) + (-(2 * MDP(23)) + t250) * pkin(4) + (t257 + 0.2e1 * t235 + 0.2e1 * t236 + t240) * qJ(5); t155 * MDP(23) + t205 * t151 + t234; t170 * MDP(25); (t205 + t249) * t194; t217; MDP(25); t151 * MDP(30) + t204; MDP(31) * t160 + MDP(32) * t161; MDP(30) * t194 + MDP(31) * t158 - MDP(32) * t159 + t211 * t197; t203; t210; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
