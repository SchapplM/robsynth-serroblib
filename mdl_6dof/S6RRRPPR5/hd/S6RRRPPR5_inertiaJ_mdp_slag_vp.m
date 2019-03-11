% Calculate joint inertia matrix for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR5_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:44:17
% EndTime: 2019-03-09 15:44:22
% DurationCPUTime: 1.45s
% Computational Cost: add. (2725->268), mult. (5997->409), div. (0->0), fcn. (6840->12), ass. (0->115)
t211 = sin(qJ(3));
t213 = cos(qJ(3));
t266 = -(MDP(16) * t211 + MDP(17) * t213) * pkin(9) + t211 * MDP(13) + t213 * MDP(14);
t254 = -qJ(4) - pkin(9);
t190 = t254 * t213;
t205 = sin(pkin(11));
t208 = cos(pkin(11));
t233 = t254 * t211;
t171 = -t190 * t205 - t208 * t233;
t265 = t171 ^ 2;
t264 = 2 * MDP(12);
t263 = 0.2e1 * MDP(16);
t262 = 2 * MDP(22);
t261 = -2 * MDP(25);
t260 = 2 * MDP(26);
t259 = 2 * MDP(30);
t258 = cos(qJ(6));
t212 = sin(qJ(2));
t257 = pkin(1) * t212;
t214 = cos(qJ(2));
t256 = pkin(1) * t214;
t196 = pkin(3) * t205 + qJ(5);
t255 = pkin(10) + t196;
t186 = t205 * t213 + t208 * t211;
t204 = sin(pkin(12));
t253 = t186 * t204;
t207 = cos(pkin(12));
t252 = t186 * t207;
t206 = sin(pkin(6));
t251 = t206 * t212;
t250 = t206 * t214;
t209 = cos(pkin(6));
t249 = t209 * MDP(8);
t234 = pkin(8) * t250;
t176 = t234 + (pkin(9) + t257) * t209;
t177 = (-pkin(2) * t214 - pkin(9) * t212 - pkin(1)) * t206;
t157 = -t176 * t211 + t213 * t177;
t179 = t211 * t209 + t213 * t251;
t148 = -pkin(3) * t250 - qJ(4) * t179 + t157;
t158 = t213 * t176 + t211 * t177;
t192 = t211 * t251;
t232 = t209 * t213 - t192;
t153 = t232 * qJ(4) + t158;
t139 = t205 * t148 + t208 * t153;
t136 = -qJ(5) * t250 + t139;
t162 = t179 * t205 - t208 * t232;
t163 = t208 * t179 + t205 * t232;
t193 = pkin(8) * t251;
t200 = -pkin(3) * t213 - pkin(2);
t168 = t192 * pkin(3) + t193 + (t200 - t256) * t209;
t142 = t162 * pkin(4) - t163 * qJ(5) + t168;
t131 = t207 * t136 + t204 * t142;
t184 = t205 * t211 - t208 * t213;
t169 = pkin(4) * t184 - qJ(5) * t186 + t200;
t173 = -t208 * t190 + t205 * t233;
t150 = t204 * t169 + t207 * t173;
t248 = t204 ^ 2 + t207 ^ 2;
t247 = MDP(15) * t214;
t246 = MDP(20) * t207;
t245 = MDP(21) * t204;
t199 = -pkin(3) * t208 - pkin(4);
t244 = MDP(23) * t199;
t210 = sin(qJ(6));
t220 = -t210 * t204 + t258 * t207;
t243 = MDP(29) * t220;
t138 = t148 * t208 - t205 * t153;
t137 = pkin(4) * t250 - t138;
t242 = t137 * MDP(23);
t155 = t163 * t204 + t207 * t250;
t156 = t163 * t207 - t204 * t250;
t143 = t258 * t155 + t156 * t210;
t241 = t143 * MDP(27);
t144 = -t210 * t155 + t258 * t156;
t240 = t144 * MDP(24);
t187 = t258 * t204 + t210 * t207;
t160 = t187 * t186;
t239 = t160 * MDP(27);
t161 = t220 * t186;
t238 = t161 * MDP(24);
t237 = t162 * MDP(28);
t236 = t184 * MDP(28);
t235 = t211 * MDP(11);
t130 = -t136 * t204 + t207 * t142;
t149 = t207 * t169 - t173 * t204;
t231 = t248 * MDP(23);
t230 = t130 * t207 + t131 * t204;
t229 = -t130 * t204 + t131 * t207;
t228 = t149 * t207 + t150 * t204;
t227 = -t149 * t204 + t150 * t207;
t224 = MDP(20) * t204 + MDP(21) * t207;
t145 = pkin(5) * t184 - pkin(10) * t252 + t149;
t146 = -pkin(10) * t253 + t150;
t133 = t258 * t145 - t210 * t146;
t134 = t210 * t145 + t258 * t146;
t223 = t133 * MDP(29) - t134 * MDP(30);
t222 = t160 * MDP(29) + t161 * MDP(30);
t221 = -t187 * MDP(30) + t243;
t219 = -t179 * MDP(13) - t232 * MDP(14);
t218 = t221 - t245 + t246;
t180 = t255 * t204;
t181 = t255 * t207;
t217 = t187 * MDP(26) + t220 * MDP(27) + (-t258 * t180 - t210 * t181) * MDP(29) - (-t210 * t180 + t258 * t181) * MDP(30);
t216 = -t224 * t196 + t217;
t202 = t206 ^ 2;
t189 = -pkin(5) * t207 + t199;
t183 = t209 * t257 + t234;
t182 = t209 * t256 - t193;
t175 = t193 + (-pkin(2) - t256) * t209;
t159 = pkin(5) * t253 + t171;
t132 = pkin(5) * t155 + t137;
t129 = -pkin(10) * t155 + t131;
t128 = pkin(5) * t162 - pkin(10) * t156 + t130;
t127 = t210 * t128 + t258 * t129;
t126 = t258 * t128 - t210 * t129;
t1 = [t202 * t212 ^ 2 * MDP(4) + (t130 ^ 2 + t131 ^ 2 + t137 ^ 2) * MDP(23) + MDP(1) + (t138 ^ 2 + t139 ^ 2 + t168 ^ 2) * MDP(19) + (0.2e1 * MDP(6) * t251 + t249) * t209 + (t179 * MDP(11) + t232 * t264) * t179 + (t237 - 0.2e1 * t241) * t162 + (t143 * t261 + t162 * t260 + t240) * t144 + ((0.2e1 * t212 * MDP(5) + t247) * t202 + 0.2e1 * (t209 * MDP(7) + t219) * t206) * t214 + 0.2e1 * (-t183 * t209 - t202 * t257) * MDP(10) + 0.2e1 * (t182 * t209 + t202 * t256) * MDP(9) + (-t157 * t250 - t175 * t232) * t263 + 0.2e1 * (t158 * t250 + t175 * t179) * MDP(17) + 0.2e1 * (t130 * t162 + t137 * t155) * MDP(20) + 0.2e1 * (t126 * t162 + t132 * t143) * MDP(29) + 0.2e1 * (-t131 * t162 + t137 * t156) * MDP(21) + (-t127 * t162 + t132 * t144) * t259 + 0.2e1 * (-t138 * t163 - t139 * t162) * MDP(18) + (-t130 * t156 - t131 * t155) * t262; t179 * t235 + t144 * t238 + t162 * t236 + t249 + t182 * MDP(9) - t183 * MDP(10) + (t179 * t213 + t211 * t232) * MDP(12) + (pkin(2) * t232 - t175 * t213) * MDP(16) + (-pkin(2) * t179 + t175 * t211) * MDP(17) + (-t138 * t186 - t139 * t184 - t173 * t162 + t171 * t163) * MDP(18) + (-t138 * t171 + t139 * t173 + t168 * t200) * MDP(19) + (t130 * t184 + t137 * t253 + t149 * t162 + t171 * t155) * MDP(20) + (-t131 * t184 + t137 * t252 - t150 * t162 + t171 * t156) * MDP(21) + (-t149 * t156 - t150 * t155 - t230 * t186) * MDP(22) + (t130 * t149 + t131 * t150 + t137 * t171) * MDP(23) + (-t161 * t143 - t144 * t160) * MDP(25) + (t144 * t184 + t161 * t162) * MDP(26) + (-t143 * t184 - t160 * t162) * MDP(27) + (t126 * t184 + t132 * t160 + t133 * t162 + t159 * t143) * MDP(29) + (-t127 * t184 + t132 * t161 - t134 * t162 + t159 * t144) * MDP(30) + (MDP(6) * t212 + (MDP(7) - t266) * t214) * t206; MDP(8) + pkin(2) * t213 * t263 + (t173 ^ 2 + t200 ^ 2 + t265) * MDP(19) + (t149 ^ 2 + t150 ^ 2 + t265) * MDP(23) + (t236 - 0.2e1 * t239) * t184 + (-0.2e1 * pkin(2) * MDP(17) + t213 * t264 + t235) * t211 + (t160 * t261 + t184 * t260 + t238) * t161 + 0.2e1 * t222 * t159 + 0.2e1 * (-t173 * MDP(18) + t149 * MDP(20) - t150 * MDP(21) + t223) * t184 + 0.2e1 * (-t228 * MDP(22) + (MDP(18) + t224) * t171) * t186; -t206 * t247 + t157 * MDP(16) - t158 * MDP(17) + (-t137 * t207 + t155 * t199) * MDP(20) + (t137 * t204 + t156 * t199) * MDP(21) + t229 * MDP(22) + t199 * t242 + t187 * t240 + (-t143 * t187 + t144 * t220) * MDP(25) + (-t132 * t220 + t143 * t189) * MDP(29) + (t132 * t187 + t144 * t189) * MDP(30) + ((-t155 * t207 + t156 * t204) * MDP(22) + t229 * MDP(23)) * t196 + t216 * t162 + ((-t162 * t205 - t163 * t208) * MDP(18) + (t138 * t208 + t139 * t205) * MDP(19)) * pkin(3) - t219; (-t171 * t207 + t199 * t253) * MDP(20) + (t171 * t204 + t199 * t252) * MDP(21) + t227 * MDP(22) + (t171 * t199 + t227 * t196) * MDP(23) + t187 * t238 + (-t160 * t187 + t161 * t220) * MDP(25) + (-t159 * t220 + t160 * t189) * MDP(29) + (t159 * t187 + t161 * t189) * MDP(30) + t216 * t184 + ((-t184 * t205 - t186 * t208) * MDP(18) + (-t171 * t208 + t173 * t205) * MDP(19)) * pkin(3) + t266; -0.2e1 * t189 * t243 + MDP(15) + (t205 ^ 2 + t208 ^ 2) * MDP(19) * pkin(3) ^ 2 + (t244 + 0.2e1 * t245 - 0.2e1 * t246) * t199 + (MDP(24) * t187 + t189 * t259 - t220 * t261) * t187 + (t231 * t196 + t248 * t262) * t196; t168 * MDP(19) + (-t155 * t204 - t156 * t207) * MDP(22) + t230 * MDP(23) + t218 * t162; -t248 * MDP(22) * t186 + t200 * MDP(19) + t228 * MDP(23) + t218 * t184; 0; MDP(19) + t231; t155 * MDP(20) + t156 * MDP(21) + t143 * MDP(29) + t144 * MDP(30) + t242; t171 * MDP(23) + t224 * t186 + t222; -t218 + t244; 0; MDP(23); t144 * MDP(26) + t126 * MDP(29) - t127 * MDP(30) + t237 - t241; t161 * MDP(26) + t223 + t236 - t239; t217; t221; 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
