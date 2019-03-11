% Calculate joint inertia matrix for
% S6RRRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP7_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:11:11
% EndTime: 2019-03-09 17:11:15
% DurationCPUTime: 1.34s
% Computational Cost: add. (2583->279), mult. (5604->405), div. (0->0), fcn. (6281->10), ass. (0->111)
t186 = sin(pkin(11));
t188 = cos(pkin(11));
t191 = sin(qJ(3));
t194 = cos(qJ(3));
t171 = t186 * t194 + t188 * t191;
t254 = 0.2e1 * t171;
t187 = sin(pkin(6));
t253 = 0.2e1 * t187;
t252 = -(MDP(16) * t191 + MDP(17) * t194) * pkin(9) + t191 * MDP(13) + t194 * MDP(14);
t251 = 2 * MDP(12);
t250 = 0.2e1 * MDP(16);
t249 = 2 * MDP(27);
t248 = 2 * MDP(28);
t247 = 2 * MDP(29);
t192 = sin(qJ(2));
t246 = pkin(1) * t192;
t195 = cos(qJ(2));
t245 = pkin(1) * t195;
t189 = cos(pkin(6));
t238 = t187 * t192;
t165 = t191 * t189 + t194 * t238;
t174 = t191 * t238;
t214 = t189 * t194 - t174;
t150 = t165 * t186 - t188 * t214;
t244 = pkin(5) * t150;
t170 = t186 * t191 - t188 * t194;
t243 = pkin(5) * t170;
t242 = -qJ(4) - pkin(9);
t241 = qJ(6) * t150;
t240 = qJ(6) * t170;
t180 = pkin(3) * t186 + pkin(10);
t239 = t170 * t180;
t237 = t187 * t195;
t236 = t189 * MDP(8);
t235 = t192 * MDP(6);
t219 = pkin(8) * t237;
t162 = t219 + (pkin(9) + t246) * t189;
t163 = (-pkin(2) * t195 - pkin(9) * t192 - pkin(1)) * t187;
t148 = -t162 * t191 + t194 * t163;
t138 = -pkin(3) * t237 - qJ(4) * t165 + t148;
t149 = t194 * t162 + t191 * t163;
t143 = qJ(4) * t214 + t149;
t131 = t186 * t138 + t188 * t143;
t129 = -pkin(10) * t237 + t131;
t151 = t188 * t165 + t186 * t214;
t176 = pkin(8) * t238;
t182 = -pkin(3) * t194 - pkin(2);
t154 = t174 * pkin(3) + t176 + (t182 - t245) * t189;
t134 = t150 * pkin(4) - t151 * pkin(10) + t154;
t190 = sin(qJ(5));
t193 = cos(qJ(5));
t125 = t193 * t129 + t190 * t134;
t155 = pkin(4) * t170 - pkin(10) * t171 + t182;
t173 = t242 * t194;
t218 = t242 * t191;
t159 = -t188 * t173 + t186 * t218;
t140 = t190 * t155 + t193 * t159;
t184 = t190 ^ 2;
t185 = t193 ^ 2;
t234 = t184 + t185;
t233 = MDP(15) * t195;
t232 = MDP(21) * t193;
t231 = MDP(28) * t171;
t230 = MDP(30) * t180;
t216 = t129 * t190 - t193 * t134;
t123 = t216 - t244;
t229 = t123 * MDP(30);
t146 = t151 * t190 + t193 * t237;
t228 = t146 * MDP(21);
t227 = t146 * MDP(23);
t147 = t151 * t193 - t190 * t237;
t226 = t147 * MDP(20);
t225 = t150 * MDP(22);
t224 = t150 * MDP(24);
t223 = t165 * MDP(11);
t222 = t170 * MDP(24);
t221 = MDP(25) + MDP(27);
t220 = -MDP(26) + MDP(29);
t181 = -pkin(3) * t188 - pkin(4);
t217 = -MDP(30) * pkin(5) - MDP(27);
t130 = t138 * t188 - t186 * t143;
t215 = -t193 * t155 + t159 * t190;
t157 = -t173 * t186 - t188 * t218;
t213 = MDP(28) + t230;
t128 = pkin(4) * t237 - t130;
t212 = MDP(25) - t217;
t211 = -pkin(5) * t193 - qJ(6) * t190;
t210 = pkin(5) * t190 - qJ(6) * t193;
t209 = MDP(30) * qJ(6) + t220;
t135 = t140 + t240;
t136 = t215 - t243;
t208 = t135 * t190 - t136 * t193;
t168 = t181 + t211;
t207 = t168 * t171 - t239;
t206 = t171 * t181 - t239;
t203 = t193 * MDP(22) - t190 * MDP(23);
t202 = -MDP(25) * t215 - t140 * MDP(26);
t126 = pkin(5) * t146 - qJ(6) * t147 + t128;
t201 = t128 * MDP(26) - t126 * MDP(29) + t226;
t200 = t147 * MDP(21) - t128 * MDP(25) - t126 * MDP(27);
t199 = -t165 * MDP(13) - MDP(14) * t214;
t198 = t190 * t220 + t193 * t221;
t197 = t147 * MDP(22) - MDP(25) * t216 - t125 * MDP(26) + t224 - t227;
t183 = t187 ^ 2;
t167 = t189 * t246 + t219;
t166 = t189 * t245 - t176;
t161 = t176 + (-pkin(2) - t245) * t189;
t145 = t190 * t146;
t144 = t171 * t210 + t157;
t122 = t125 + t241;
t1 = [(t122 ^ 2 + t123 ^ 2 + t126 ^ 2) * MDP(30) + (t130 ^ 2 + t131 ^ 2 + t154 ^ 2) * MDP(19) + t183 * t192 ^ 2 * MDP(4) + MDP(1) + (t235 * t253 + t236) * t189 + (t214 * t251 + t223) * t165 + (t224 - 0.2e1 * t227) * t150 + (0.2e1 * t225 + t226 - 0.2e1 * t228) * t147 + ((0.2e1 * t192 * MDP(5) + t233) * t183 + (t189 * MDP(7) + t199) * t253) * t195 + (-t148 * t237 - t161 * t214) * t250 + 0.2e1 * (t149 * t237 + t161 * t165) * MDP(17) + 0.2e1 * (t166 * t189 + t183 * t245) * MDP(9) + 0.2e1 * (-t167 * t189 - t183 * t246) * MDP(10) + 0.2e1 * (-t125 * t150 + t128 * t147) * MDP(26) + (-t123 * t150 + t126 * t146) * t249 + 0.2e1 * (-t130 * t151 - t131 * t150) * MDP(18) + (-t122 * t146 + t123 * t147) * t248 + (t122 * t150 - t126 * t147) * t247 + 0.2e1 * (t128 * t146 - t150 * t216) * MDP(25); t236 + t166 * MDP(9) - t167 * MDP(10) + (t165 * t194 + t191 * t214) * MDP(12) + (pkin(2) * t214 - t161 * t194) * MDP(16) + (-pkin(2) * t165 + t161 * t191) * MDP(17) + (-t150 * t159 + t151 * t157) * MDP(18) + (-t130 * t157 + t131 * t159 + t154 * t182) * MDP(19) + (t146 * t157 - t150 * t215) * MDP(25) + (-t140 * t150 + t147 * t157) * MDP(26) + (-t136 * t150 + t144 * t146) * MDP(27) + (-t135 * t146 + t136 * t147) * MDP(28) + (t135 * t150 - t144 * t147) * MDP(29) + (t122 * t135 + t123 * t136 + t126 * t144) * MDP(30) + t191 * t223 + (-t131 * MDP(18) - t123 * MDP(27) + t122 * MDP(29) + t197) * t170 + (t235 + (MDP(7) - t252) * t195) * t187 + (-t130 * MDP(18) + (-t150 * MDP(23) - t122 * MDP(28) - t200) * t190 + (t123 * MDP(28) + t201 + t225 - t228) * t193) * t171; MDP(8) + pkin(2) * t194 * t250 + (t157 ^ 2 + t159 ^ 2 + t182 ^ 2) * MDP(19) + (t135 ^ 2 + t136 ^ 2 + t144 ^ 2) * MDP(30) + (t185 * MDP(20) - 0.2e1 * t190 * t232) * t171 ^ 2 + (MDP(11) * t191 - 0.2e1 * pkin(2) * MDP(17) + t194 * t251) * t191 + (t203 * t254 + t222) * t170 + 0.2e1 * (-t159 * MDP(18) - t136 * MDP(27) + t135 * MDP(29) + t202) * t170 + (-t208 * MDP(28) + (t190 * MDP(27) - t193 * MDP(29)) * t144 + (t190 * MDP(25) + t193 * MDP(26) + MDP(18)) * t157) * t254; -t187 * t233 + t148 * MDP(16) - t149 * MDP(17) - t145 * MDP(21) + (t146 * MDP(25) + t147 * MDP(26)) * t181 + (t146 * MDP(27) - t147 * MDP(29) + t126 * MDP(30)) * t168 + ((-t146 * t180 + t122) * MDP(28) + t122 * t230 + (t180 * t220 + MDP(23)) * t150 + t200) * t193 + ((t147 * t180 + t123) * MDP(28) + t180 * t229 + (-t180 * t221 + MDP(22)) * t150 + t201) * t190 + ((-t150 * t186 - t151 * t188) * MDP(18) + (t130 * t188 + t131 * t186) * MDP(19)) * pkin(3) - t199; t144 * t168 * MDP(30) + (-t184 + t185) * MDP(21) * t171 + (t170 * MDP(23) - t157 * MDP(25) + MDP(26) * t206 - t144 * MDP(27) - MDP(29) * t207 + t135 * t213) * t193 + (t193 * t171 * MDP(20) + t170 * MDP(22) + MDP(25) * t206 + t157 * MDP(26) + MDP(27) * t207 - t144 * MDP(29) + t136 * t213) * t190 + ((-t170 * t186 - t171 * t188) * MDP(18) + (-t157 * t188 + t159 * t186) * MDP(19)) * pkin(3) + t252; MDP(15) + t184 * MDP(20) + (t180 ^ 2 * t234 + t168 ^ 2) * MDP(30) + (t186 ^ 2 + t188 ^ 2) * MDP(19) * pkin(3) ^ 2 + t234 * t180 * t248 + 0.2e1 * (-MDP(25) * t181 - MDP(27) * t168) * t193 + 0.2e1 * (MDP(26) * t181 - MDP(29) * t168 + t232) * t190; t154 * MDP(19) + (-t147 * t193 - t145) * MDP(28) + (t122 * t190 - t123 * t193) * MDP(30) + t198 * t150; t182 * MDP(19) + MDP(30) * t208 + t170 * t198 - t231 * t234; 0; MDP(30) * t234 + MDP(19); (-t216 + 0.2e1 * t244) * MDP(27) + (-pkin(5) * t147 - qJ(6) * t146) * MDP(28) + (t125 + 0.2e1 * t241) * MDP(29) + (-pkin(5) * t123 + qJ(6) * t122) * MDP(30) + t197; t222 + (-t215 + 0.2e1 * t243) * MDP(27) + (t140 + 0.2e1 * t240) * MDP(29) + (-pkin(5) * t136 + qJ(6) * t135) * MDP(30) + (MDP(28) * t211 + t203) * t171 + t202; t190 * MDP(22) + t193 * MDP(23) - t210 * MDP(28) + (-t190 * t212 + t193 * t209) * t180; t190 * t209 + t193 * t212; MDP(24) + pkin(5) * t249 + qJ(6) * t247 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(30); -t150 * MDP(27) + t147 * MDP(28) + t229; -t170 * MDP(27) + t136 * MDP(30) + t193 * t231; t213 * t190; -t193 * MDP(30); t217; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
