% Calculate joint inertia matrix for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP6_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:14:14
% EndTime: 2019-03-09 12:14:18
% DurationCPUTime: 1.39s
% Computational Cost: add. (2124->265), mult. (4961->389), div. (0->0), fcn. (5416->10), ass. (0->103)
t212 = MDP(26) - MDP(29);
t213 = MDP(25) + MDP(27);
t242 = 2 * MDP(27);
t241 = 2 * MDP(28);
t190 = cos(qJ(2));
t240 = pkin(1) * t190;
t181 = sin(pkin(11));
t182 = sin(pkin(6));
t183 = cos(pkin(11));
t187 = sin(qJ(2));
t160 = (t181 * t190 + t183 * t187) * t182;
t184 = cos(pkin(6));
t186 = sin(qJ(4));
t189 = cos(qJ(4));
t151 = t160 * t186 - t184 * t189;
t239 = pkin(5) * t151;
t238 = pkin(8) + qJ(3);
t237 = qJ(6) * t151;
t173 = pkin(2) * t181 + pkin(9);
t185 = sin(qJ(5));
t235 = t173 * t185;
t188 = cos(qJ(5));
t234 = t173 * t188;
t176 = t182 ^ 2;
t233 = t176 * t187;
t232 = t182 * t187;
t231 = t182 * t190;
t230 = t184 * MDP(8);
t229 = t185 * t186;
t228 = t185 * t189;
t227 = t186 * t188;
t170 = t184 * t240;
t153 = pkin(2) * t184 - t232 * t238 + t170;
t211 = pkin(1) * t184 * t187;
t157 = t231 * t238 + t211;
t141 = t153 * t181 + t157 * t183;
t138 = pkin(9) * t184 + t141;
t159 = t181 * t232 - t183 * t231;
t166 = (-pkin(2) * t190 - pkin(1)) * t182;
t144 = pkin(3) * t159 - pkin(9) * t160 + t166;
t135 = t138 * t189 + t144 * t186;
t130 = pkin(10) * t159 + t135;
t140 = t153 * t183 - t157 * t181;
t137 = -pkin(3) * t184 - t140;
t152 = t160 * t189 + t184 * t186;
t133 = pkin(4) * t151 - pkin(10) * t152 + t137;
t126 = t130 * t188 + t133 * t185;
t174 = -pkin(2) * t183 - pkin(3);
t165 = -pkin(4) * t189 - pkin(10) * t186 + t174;
t150 = t165 * t185 + t189 * t234;
t177 = t185 ^ 2;
t179 = t188 ^ 2;
t226 = t177 + t179;
t225 = MDP(21) * t188;
t142 = t152 * t185 - t159 * t188;
t224 = t142 * MDP(23);
t143 = t152 * t188 + t159 * t185;
t223 = t143 * MDP(20);
t222 = t143 * MDP(22);
t221 = t151 * MDP(24);
t220 = t152 * MDP(13);
t219 = t152 * MDP(14);
t218 = t159 * MDP(15);
t205 = -pkin(5) * t188 - qJ(6) * t185;
t168 = -pkin(4) + t205;
t217 = t168 * MDP(30);
t216 = t174 * MDP(19);
t215 = t185 * MDP(22);
t214 = t189 * MDP(18);
t145 = t151 * t227;
t210 = t226 * pkin(10);
t209 = -MDP(30) * pkin(5) - MDP(27);
t208 = t130 * t185 - t133 * t188;
t134 = -t138 * t186 + t144 * t189;
t207 = -t173 * MDP(19) + MDP(16);
t206 = -0.2e1 * qJ(6) * MDP(29) - MDP(24);
t204 = -pkin(5) * t185 + qJ(6) * t188;
t123 = t126 + t237;
t124 = t208 - t239;
t203 = t123 * t188 + t124 * t185;
t146 = -qJ(6) * t189 + t150;
t162 = t188 * t165;
t147 = -t162 + (pkin(5) + t235) * t189;
t202 = t146 * t188 + t147 * t185;
t201 = t188 * MDP(22) - t185 * MDP(23);
t200 = t188 * MDP(23) + t215;
t199 = (-t173 * t228 + t162) * MDP(25) - t150 * MDP(26);
t198 = t185 * MDP(27) - t188 * MDP(29);
t129 = -pkin(4) * t159 - t134;
t196 = (MDP(6) * t187 + MDP(7) * t190) * t182;
t195 = MDP(27) * t147 - MDP(29) * t146 - t199;
t194 = -MDP(25) * t208 - t126 * MDP(26) + t221 + t222 - t224;
t193 = (MDP(30) * qJ(6) - t212) * t188 + (-MDP(25) + t209) * t185;
t180 = t189 ^ 2;
t178 = t186 ^ 2;
t175 = t179 * t186;
t171 = pkin(10) * t228;
t164 = pkin(8) * t231 + t211;
t163 = -pkin(8) * t232 + t170;
t156 = (t173 - t204) * t186;
t139 = t142 * t227;
t127 = pkin(5) * t142 - qJ(6) * t143 + t129;
t1 = [(t140 ^ 2 + t141 ^ 2 + t166 ^ 2) * MDP(12) + t159 ^ 2 * MDP(17) + (t123 ^ 2 + t124 ^ 2 + t127 ^ 2) * MDP(30) + MDP(1) + (MDP(4) * t187 + 0.2e1 * MDP(5) * t190) * t233 + (0.2e1 * t218 + t220) * t152 + (-0.2e1 * t142 * MDP(21) + t223) * t143 + (0.2e1 * t196 + t230) * t184 + (-0.2e1 * MDP(16) * t159 - 0.2e1 * t219 + t221 + 0.2e1 * t222 - 0.2e1 * t224) * t151 + (-t123 * t142 + t124 * t143) * t241 + 0.2e1 * (t129 * t142 - t151 * t208) * MDP(25) + 0.2e1 * (t123 * t151 - t127 * t143) * MDP(29) + (-t124 * t151 + t127 * t142) * t242 + 0.2e1 * (-t126 * t151 + t129 * t143) * MDP(26) + 0.2e1 * (t134 * t159 + t137 * t151) * MDP(18) + 0.2e1 * (-t135 * t159 + t137 * t152) * MDP(19) + 0.2e1 * (-t140 * t160 - t141 * t159) * MDP(11) + 0.2e1 * (-pkin(1) * t233 - t164 * t184) * MDP(10) + 0.2e1 * (t163 * t184 + t176 * t240) * MDP(9); t230 + t163 * MDP(9) - t164 * MDP(10) + t152 * t216 - t139 * MDP(21) + t145 * MDP(22) + (-t142 * t146 + t143 * t147) * MDP(28) + (t123 * t146 + t124 * t147) * MDP(30) + t196 + (MDP(27) * t142 - MDP(29) * t143 + MDP(30) * t127) * t156 + (MDP(18) * t174 - t195) * t151 + (-t137 * MDP(18) + t124 * MDP(27) - t123 * MDP(29) + t159 * t207 - t194 + t219) * t189 + ((-t159 * t181 - t160 * t183) * MDP(11) + (t140 * t183 + t141 * t181) * MDP(12)) * pkin(2) + (t220 - t151 * MDP(14) + t218 + t137 * MDP(19) + (-MDP(18) * t159 + MDP(25) * t142 + MDP(26) * t143) * t173 + (t129 * MDP(26) + MDP(28) * t124 - t127 * MDP(29) + t223) * t188 + (-t143 * MDP(21) - t151 * MDP(23) + t129 * MDP(25) + t127 * MDP(27) - MDP(28) * t123) * t185) * t186; MDP(8) - 0.2e1 * t174 * t214 + t180 * MDP(24) + (t146 ^ 2 + t147 ^ 2 + t156 ^ 2) * MDP(30) + (t181 ^ 2 + t183 ^ 2) * MDP(12) * pkin(2) ^ 2 + 0.2e1 * t195 * t189 + (MDP(20) * t179 - 0.2e1 * t185 * t225 + MDP(13) + 0.2e1 * (t185 * MDP(25) + t188 * MDP(26)) * t173) * t178 + 0.2e1 * (t216 + (MDP(14) - t201) * t189 + (-t146 * t185 + t147 * t188) * MDP(28) + t198 * t156) * t186; t166 * MDP(12) + (t143 * t229 - t139) * MDP(28) + (-t127 * t189 + t186 * t203) * MDP(30) + (-t186 * MDP(19) + t214) * t159 + t213 * (-t142 * t189 - t151 * t229) - t212 * (t143 * t189 + t145); (-t156 * t189 + t186 * t202) * MDP(30); MDP(12) + (t178 * t226 + t180) * MDP(30); t152 * MDP(15) + t159 * MDP(17) + t134 * MDP(18) - t135 * MDP(19) + t185 * t223 + (-t142 * t185 + t143 * t188) * MDP(21) + (-pkin(4) * t142 - t129 * t188) * MDP(25) + (-pkin(4) * t143 + t129 * t185) * MDP(26) + (-t127 * t188 + t142 * t168) * MDP(27) + t203 * MDP(28) + (-t127 * t185 - t143 * t168) * MDP(29) + t127 * t217 + ((-t142 * t188 + t143 * t185) * MDP(28) + t203 * MDP(30)) * pkin(10) + (-MDP(16) + (-t185 * t213 - t188 * t212) * pkin(10) + t200) * t151; t175 * MDP(21) + t171 * MDP(25) + (-t156 * t188 + t171) * MDP(27) + t202 * MDP(28) - t156 * t185 * MDP(29) + (pkin(10) * t202 + t156 * t168) * MDP(30) + (-t215 + (pkin(10) * t212 - MDP(23)) * t188 + t207) * t189 + (MDP(15) - t173 * MDP(18) + t188 * t185 * MDP(20) - t177 * MDP(21) + (-pkin(4) * t185 - t234) * MDP(25) + (-pkin(4) * t188 + t235) * MDP(26) + t198 * t168) * t186; t175 * MDP(28) + (t177 * MDP(28) + MDP(30) * t210 - MDP(19)) * t186 + (-t185 * t212 + t188 * t213 + MDP(18) - t217) * t189; MDP(17) + t177 * MDP(20) + (pkin(10) ^ 2 * t226 + t168 ^ 2) * MDP(30) + t210 * t241 + 0.2e1 * (MDP(25) * pkin(4) - MDP(27) * t168) * t188 + 0.2e1 * (-MDP(26) * pkin(4) - MDP(29) * t168 + t225) * t185; (-t208 + 0.2e1 * t239) * MDP(27) + (-pkin(5) * t143 - qJ(6) * t142) * MDP(28) + (t126 + 0.2e1 * t237) * MDP(29) + (-pkin(5) * t124 + qJ(6) * t123) * MDP(30) + t194; t162 * MDP(27) + t150 * MDP(29) + (-pkin(5) * t147 + qJ(6) * t146) * MDP(30) + ((-0.2e1 * pkin(5) - t235) * MDP(27) + t206) * t189 + (MDP(28) * t205 + t201) * t186 + t199; t193 * t186; MDP(28) * t204 + pkin(10) * t193 + t200; pkin(5) * t242 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(30) - t206; -MDP(27) * t151 + MDP(28) * t143 + MDP(30) * t124; MDP(27) * t189 + MDP(28) * t227 + MDP(30) * t147; MDP(30) * t229; (MDP(30) * pkin(10) + MDP(28)) * t185; t209; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
