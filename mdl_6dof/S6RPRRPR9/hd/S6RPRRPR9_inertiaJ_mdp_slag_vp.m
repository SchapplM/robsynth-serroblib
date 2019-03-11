% Calculate joint inertia matrix for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRPR9_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:31:38
% EndTime: 2019-03-09 05:31:41
% DurationCPUTime: 1.28s
% Computational Cost: add. (3163->236), mult. (8316->385), div. (0->0), fcn. (9662->14), ass. (0->104)
t194 = sin(pkin(6));
t193 = sin(pkin(7));
t238 = cos(pkin(6));
t213 = t238 * t193;
t196 = cos(pkin(12));
t237 = cos(pkin(7));
t214 = t196 * t237;
t244 = t194 * t214 + t213;
t192 = sin(pkin(12));
t217 = pkin(1) * t238;
t231 = t194 * t196;
t171 = qJ(2) * t231 + t192 * t217;
t158 = pkin(9) * t244 + t171;
t199 = sin(qJ(3));
t202 = cos(qJ(3));
t182 = t196 * t217;
t234 = t192 * t194;
t161 = t238 * pkin(2) + t182 + (-pkin(9) * t237 - qJ(2)) * t234;
t166 = (-pkin(9) * t192 * t193 - pkin(2) * t196 - pkin(1)) * t194;
t209 = t161 * t237 + t166 * t193;
t142 = -t158 * t199 + t202 * t209;
t243 = 2 * MDP(22);
t242 = 2 * MDP(29);
t241 = 2 * MDP(30);
t188 = t194 ^ 2;
t240 = pkin(1) * t188;
t239 = -qJ(5) - pkin(10);
t191 = sin(pkin(13));
t195 = cos(pkin(13));
t198 = sin(qJ(4));
t201 = cos(qJ(4));
t175 = t191 * t201 + t195 * t198;
t197 = sin(qJ(6));
t236 = t175 * t197;
t200 = cos(qJ(6));
t235 = t175 * t200;
t233 = t193 * t199;
t232 = t193 * t202;
t146 = -t161 * t193 + t166 * t237;
t159 = t199 * t234 - t202 * t244;
t160 = t199 * t213 + (t192 * t202 + t199 * t214) * t194;
t137 = pkin(3) * t159 - pkin(10) * t160 + t146;
t143 = t202 * t158 + t199 * t209;
t168 = t193 * t231 - t237 * t238;
t139 = -t168 * pkin(10) + t143;
t131 = t137 * t201 - t139 * t198;
t149 = t160 * t201 - t168 * t198;
t128 = pkin(4) * t159 - qJ(5) * t149 + t131;
t132 = t137 * t198 + t139 * t201;
t148 = t160 * t198 + t168 * t201;
t130 = -qJ(5) * t148 + t132;
t125 = t128 * t191 + t130 * t195;
t230 = MDP(15) * t198;
t229 = MDP(20) * t201;
t187 = -pkin(4) * t201 - pkin(3);
t228 = MDP(23) * t187;
t141 = -t148 * t191 + t149 * t195;
t135 = t141 * t200 + t159 * t197;
t227 = MDP(24) * t135;
t226 = MDP(24) * t197;
t134 = t141 * t197 - t159 * t200;
t225 = MDP(25) * t134;
t140 = t148 * t195 + t149 * t191;
t224 = MDP(26) * t140;
t223 = MDP(27) * t134;
t222 = MDP(28) * t140;
t221 = t148 * MDP(18);
t220 = t149 * MDP(17);
t219 = t159 * MDP(19);
t174 = t191 * t198 - t195 * t201;
t218 = t174 * MDP(28);
t216 = t197 * t200 * MDP(25);
t215 = t239 * t198;
t124 = t128 * t195 - t130 * t191;
t211 = MDP(29) * t200 - MDP(30) * t197;
t210 = MDP(29) * t197 + MDP(30) * t200;
t208 = (MDP(26) * t200 - MDP(27) * t197) * t175;
t207 = -t198 * t233 + t201 * t237;
t138 = t168 * pkin(3) - t142;
t206 = t198 * MDP(17) + t201 * MDP(18) + (-MDP(20) * t198 - MDP(21) * t201) * pkin(10);
t205 = t197 * MDP(26) + t200 * MDP(27) - (pkin(4) * t191 + pkin(11)) * t210;
t123 = pkin(11) * t159 + t125;
t133 = pkin(4) * t148 + t138;
t126 = pkin(5) * t140 - pkin(11) * t141 + t133;
t120 = -t123 * t197 + t126 * t200;
t121 = t123 * t200 + t126 * t197;
t204 = MDP(26) * t135 + MDP(29) * t120 - MDP(30) * t121 + t222 - t223;
t190 = t200 ^ 2;
t189 = t197 ^ 2;
t186 = -pkin(4) * t195 - pkin(5);
t178 = t239 * t201;
t172 = t198 * t237 + t201 * t233;
t170 = -qJ(2) * t234 + t182;
t164 = -t195 * t178 + t191 * t215;
t162 = -t178 * t191 - t195 * t215;
t156 = pkin(5) * t174 - pkin(11) * t175 + t187;
t154 = t172 * t195 + t191 * t207;
t152 = t172 * t191 - t195 * t207;
t151 = t154 * t200 - t197 * t232;
t150 = -t154 * t197 - t200 * t232;
t145 = t156 * t197 + t164 * t200;
t144 = t156 * t200 - t164 * t197;
t122 = -pkin(5) * t159 - t124;
t1 = [t168 ^ 2 * MDP(12) + MDP(1) + (t124 ^ 2 + t125 ^ 2 + t133 ^ 2) * MDP(23) + (pkin(1) ^ 2 * t188 + t170 ^ 2 + t171 ^ 2) * MDP(7) + (-0.2e1 * MDP(10) * t168 + MDP(8) * t160) * t160 + (MDP(15) * t149 - 0.2e1 * MDP(16) * t148) * t149 + (t222 - 0.2e1 * t223) * t140 + (0.2e1 * t224 - 0.2e1 * t225 + t227) * t135 + (0.2e1 * MDP(11) * t168 - 0.2e1 * MDP(9) * t160 + t219 + 0.2e1 * t220 - 0.2e1 * t221) * t159 + (t120 * t140 + t122 * t134) * t242 + (-t121 * t140 + t122 * t135) * t241 + (-t124 * t141 - t125 * t140) * t243 + 0.2e1 * (t131 * t159 + t138 * t148) * MDP(20) + 0.2e1 * (-t132 * t159 + t138 * t149) * MDP(21) + 0.2e1 * (t143 * t168 + t146 * t160) * MDP(14) + 0.2e1 * (-t142 * t168 + t146 * t159) * MDP(13) + 0.2e1 * (-t171 * t238 - t192 * t240) * MDP(5) + 0.2e1 * (t170 * t238 + t196 * t240) * MDP(4) + 0.2e1 * (-t170 * t192 + t171 * t196) * MDP(6) * t194; (t159 * t237 - t168 * t232) * MDP(13) + (t160 * t237 + t168 * t233) * MDP(14) + (-t148 * t232 + t159 * t207) * MDP(20) + (-t149 * t232 - t159 * t172) * MDP(21) + (-t140 * t154 + t141 * t152) * MDP(22) + (-t124 * t152 + t125 * t154 - t133 * t232) * MDP(23) + (t134 * t152 + t140 * t150) * MDP(29) + (t135 * t152 - t140 * t151) * MDP(30) + (-MDP(4) * t196 + MDP(5) * t192 - MDP(7) * pkin(1)) * t194; MDP(7) + (t193 ^ 2 * t202 ^ 2 + t152 ^ 2 + t154 ^ 2) * MDP(23); t160 * MDP(10) - t168 * MDP(12) + t142 * MDP(13) - t143 * MDP(14) + t149 * t230 + (-t148 * t198 + t149 * t201) * MDP(16) + (-pkin(3) * t148 - t138 * t201) * MDP(20) + (-pkin(3) * t149 + t138 * t198) * MDP(21) + (-t140 * t164 + t141 * t162) * MDP(22) + (-t124 * t162 + t125 * t164 + t133 * t187) * MDP(23) + (t134 * t162 + t140 * t144) * MDP(29) + (t135 * t162 - t140 * t145) * MDP(30) + (-MDP(22) * t125 + t204) * t174 + (-t124 * MDP(22) + (-MDP(25) * t135 - MDP(27) * t140 + MDP(29) * t122) * t197 + (MDP(30) * t122 + t224 - t225 + t227) * t200) * t175 + (-MDP(11) + t206) * t159; (t152 * t175 - t154 * t174) * MDP(22) + (t152 * t162 + t154 * t164) * MDP(23) + (t150 * t174 + t152 * t236) * MDP(29) + (-t151 * t174 + t152 * t235) * MDP(30) + (-t199 * MDP(14) + (-MDP(21) * t198 + MDP(13) - t228 + t229) * t202) * t193; MDP(12) + 0.2e1 * pkin(3) * t229 + (t162 ^ 2 + t164 ^ 2 + t187 ^ 2) * MDP(23) + (MDP(24) * t190 - 0.2e1 * t216) * t175 ^ 2 + (0.2e1 * MDP(16) * t201 - 0.2e1 * MDP(21) * pkin(3) + t230) * t198 + (0.2e1 * t208 + t218) * t174 + (t162 * t175 - t164 * t174) * t243 + (t144 * t174 + t162 * t236) * t242 + (-t145 * t174 + t162 * t235) * t241; t220 - t221 + t219 + t131 * MDP(20) - t132 * MDP(21) + t135 * t226 + (-t134 * t197 + t135 * t200) * MDP(25) + (-t122 * t200 + t134 * t186) * MDP(29) + (t122 * t197 + t135 * t186) * MDP(30) + t205 * t140 + ((-t140 * t191 - t141 * t195) * MDP(22) + (t124 * t195 + t125 * t191) * MDP(23)) * pkin(4); t207 * MDP(20) - t172 * MDP(21) - t211 * t152 + (-t152 * t195 + t154 * t191) * MDP(23) * pkin(4); -t211 * t162 + t205 * t174 + (t200 * t226 + (-t189 + t190) * MDP(25) + t210 * t186) * t175 + ((-t174 * t191 - t175 * t195) * MDP(22) + (-t162 * t195 + t164 * t191) * MDP(23)) * pkin(4) + t206; 0.2e1 * t216 + t189 * MDP(24) + MDP(19) + (t191 ^ 2 + t195 ^ 2) * MDP(23) * pkin(4) ^ 2 - 0.2e1 * t211 * t186; t133 * MDP(23) + t140 * t211; -MDP(23) * t232; t174 * t211 + t228; 0; MDP(23); t204; MDP(29) * t150 - MDP(30) * t151; t144 * MDP(29) - t145 * MDP(30) + t208 + t218; t205; t211; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
