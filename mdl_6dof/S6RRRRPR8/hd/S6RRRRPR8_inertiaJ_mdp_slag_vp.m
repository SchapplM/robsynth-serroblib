% Calculate joint inertia matrix for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR8_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:44:52
% EndTime: 2019-03-09 22:44:56
% DurationCPUTime: 1.31s
% Computational Cost: add. (1257->229), mult. (2353->315), div. (0->0), fcn. (2399->8), ass. (0->102)
t187 = sin(qJ(3));
t191 = cos(qJ(3));
t234 = -pkin(9) - pkin(8);
t166 = t234 * t187;
t167 = t234 * t191;
t186 = sin(qJ(4));
t190 = cos(qJ(4));
t143 = -t190 * t166 - t167 * t186;
t144 = t166 * t186 - t167 * t190;
t227 = t190 * t191;
t158 = t186 * t187 - t227;
t159 = t186 * t191 + t187 * t190;
t127 = pkin(10) * t158 + t144;
t185 = sin(qJ(6));
t189 = cos(qJ(6));
t133 = -t189 * t158 + t159 * t185;
t134 = t158 * t185 + t159 * t189;
t199 = -pkin(10) * t159 + t143;
t195 = -(t189 * t127 + t185 * t199) * MDP(35) - t133 * MDP(32) + t134 * MDP(31) - (t127 * t185 - t189 * t199) * MDP(34);
t212 = MDP(23) + MDP(25);
t251 = MDP(24) - MDP(27);
t196 = t159 * MDP(20) - t158 * MDP(21) - t212 * t143 - t251 * t144 - t195;
t256 = t196 - (MDP(16) * t187 + MDP(17) * t191) * pkin(8) + t187 * MDP(13) + t191 * MDP(14);
t188 = sin(qJ(2));
t152 = t159 * t188;
t229 = t187 * t188;
t153 = -t186 * t229 + t188 * t227;
t254 = -t153 * MDP(20) + t152 * MDP(21);
t253 = MDP(34) * t189 - t185 * MDP(35);
t192 = cos(qJ(2));
t164 = -pkin(2) * t192 - pkin(8) * t188 - pkin(1);
t157 = t191 * t164;
t230 = pkin(9) * t188;
t233 = pkin(7) * t187;
t135 = -t191 * t230 + t157 + (-pkin(3) - t233) * t192;
t231 = pkin(7) * t192;
t211 = t191 * t231;
t138 = t211 + (t164 - t230) * t187;
t119 = t186 * t135 + t190 * t138;
t225 = -t190 * t135 + t186 * t138;
t180 = t192 * pkin(4);
t116 = t180 + t225;
t107 = pkin(5) * t192 - pkin(10) * t153 + t116;
t115 = -qJ(5) * t192 + t119;
t110 = pkin(10) * t152 + t115;
t104 = t107 * t189 - t185 * t110;
t105 = t185 * t107 + t189 * t110;
t124 = -t189 * t152 + t153 * t185;
t125 = t152 * t185 + t153 * t189;
t248 = t125 * MDP(31) - t124 * MDP(32);
t243 = t104 * MDP(34) - t105 * MDP(35) + t248;
t250 = -t225 * MDP(23) - t251 * t119 - t243 - t254;
t214 = t186 * MDP(24);
t249 = (MDP(23) * t190 - t214) * pkin(3);
t246 = MDP(25) + t253;
t242 = -2 * MDP(19);
t241 = 0.2e1 * MDP(24);
t240 = 0.2e1 * MDP(25);
t239 = 2 * MDP(26);
t238 = 2 * MDP(27);
t237 = -2 * MDP(30);
t236 = 0.2e1 * MDP(34);
t235 = 0.2e1 * MDP(35);
t193 = -pkin(4) - pkin(5);
t232 = pkin(7) * t191;
t228 = t187 * t191;
t178 = t186 * pkin(3);
t171 = t178 + qJ(5);
t226 = qJ(5) + t171;
t161 = pkin(3) * t229 + t188 * pkin(7);
t224 = MDP(18) * t159;
t223 = MDP(29) * t134;
t176 = t189 * t193;
t222 = MDP(34) * (qJ(5) * t185 - t176);
t220 = MDP(35) * (t189 * qJ(5) + t185 * t193);
t219 = qJ(5) * MDP(27);
t173 = pkin(3) * t190 + pkin(4);
t170 = -pkin(5) - t173;
t165 = t189 * t170;
t217 = (t171 * t185 - t165) * MDP(34);
t216 = (t185 * t170 + t189 * t171) * MDP(35);
t215 = t173 * MDP(25);
t213 = MDP(22) + MDP(33);
t175 = -pkin(3) * t191 - pkin(2);
t210 = MDP(12) * t228;
t209 = MDP(15) + t213;
t207 = qJ(5) * t153 - t161;
t206 = pkin(4) * t240 + t213;
t204 = MDP(13) * t191 - MDP(14) * t187;
t200 = qJ(5) * t159 - t175;
t198 = -MDP(33) - t216 - t217;
t197 = -MDP(33) - t220 - t222;
t183 = t191 ^ 2;
t182 = t188 ^ 2;
t181 = t187 ^ 2;
t147 = t164 * t187 + t211;
t146 = -t187 * t231 + t157;
t132 = pkin(4) * t158 - t200;
t121 = t193 * t158 + t200;
t120 = pkin(4) * t152 - t207;
t117 = t193 * t152 + t207;
t1 = [(t115 ^ 2 + t116 ^ 2 + t120 ^ 2) * MDP(28) + MDP(1) - 0.2e1 * pkin(1) * t188 * MDP(10) + (MDP(18) * t153 + t152 * t242) * t153 + (MDP(29) * t125 + t124 * t237) * t125 + t209 * t192 ^ 2 + (t183 * MDP(11) + MDP(4) - 0.2e1 * t210) * t182 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t204) * t188 + t248 + t254) * t192 + (-t115 * t152 + t116 * t153) * t239 + (t104 * t192 + t117 * t124) * t236 + (-t105 * t192 + t117 * t125) * t235 + (t116 * t192 + t120 * t152) * t240 + (-t115 * t192 - t120 * t153) * t238 + 0.2e1 * (t152 * t161 + t192 * t225) * MDP(23) + (t119 * t192 + t153 * t161) * t241 + 0.2e1 * (-t146 * t192 + t182 * t233) * MDP(16) + 0.2e1 * (t147 * t192 + t182 * t232) * MDP(17); (-t115 * t158 + t116 * t159 + t143 * t153 - t144 * t152) * MDP(26) + (-t152 * t159 - t153 * t158) * MDP(19) + (t115 * t144 + t116 * t143 + t120 * t132) * MDP(28) + (-t124 * t134 - t125 * t133) * MDP(30) + (-t120 * t159 - t132 * t153) * MDP(27) + (t153 * t175 + t159 * t161) * MDP(24) + (t152 * t175 + t158 * t161) * MDP(23) + (t120 * t158 + t132 * t152) * MDP(25) + (t117 * t133 + t121 * t124) * MDP(34) + (t117 * t134 + t121 * t125) * MDP(35) + t125 * t223 + t153 * t224 + ((-t181 + t183) * MDP(12) + MDP(6) + MDP(11) * t228 + (-pkin(2) * t187 - t232) * MDP(16) + (-pkin(2) * t191 + t233) * MDP(17) - pkin(7) * MDP(9)) * t188 + (-pkin(7) * MDP(10) + MDP(7) - t256) * t192; MDP(8) + t181 * MDP(11) + 0.2e1 * t210 + (t132 ^ 2 + t143 ^ 2 + t144 ^ 2) * MDP(28) + t121 * t133 * t236 + (-0.2e1 * MDP(27) * t132 + t143 * t239 + t158 * t242 + t175 * t241 + t224) * t159 + (t121 * t235 + t133 * t237 + t223) * t134 + 0.2e1 * (MDP(16) * t191 - MDP(17) * t187) * pkin(2) + 0.2e1 * (MDP(23) * t175 + MDP(25) * t132 - MDP(26) * t144) * t158; t146 * MDP(16) - t147 * MDP(17) - t116 * MDP(25) + (-t152 * t171 - t153 * t173) * MDP(26) + (t115 * t171 - t116 * t173) * MDP(28) + t204 * t188 + (-t226 * MDP(27) - MDP(15) - MDP(22) + t198 - t215 - t249) * t192 + t250; (-t158 * t171 - t159 * t173) * MDP(26) + (-t143 * t173 + t144 * t171) * MDP(28) + t256; (t171 ^ 2 + t173 ^ 2) * MDP(28) + 0.2e1 * t249 + 0.2e1 * t215 + t171 * t238 + 0.2e1 * t217 + 0.2e1 * t216 + t209; (-0.2e1 * t180 - t225) * MDP(25) + (-pkin(4) * t153 - qJ(5) * t152) * MDP(26) + (-pkin(4) * t116 + qJ(5) * t115) * MDP(28) + (-MDP(22) + t197 - 0.2e1 * t219) * t192 + t250; (-pkin(4) * t159 - qJ(5) * t158) * MDP(26) + (-pkin(4) * t143 + qJ(5) * t144) * MDP(28) + t196; (0.2e1 * qJ(5) + t178) * MDP(27) + (pkin(4) * t173 + qJ(5) * t171) * MDP(28) + (-t165 - t176) * MDP(34) + t226 * MDP(35) * t189 + (t226 * MDP(34) + (t170 + t193) * MDP(35)) * t185 + (t212 * t190 - t214) * pkin(3) + t206; 0.2e1 * t219 + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(28) + 0.2e1 * t222 + 0.2e1 * t220 + t206; t153 * MDP(26) + MDP(28) * t116 + t192 * t246; MDP(26) * t159 + MDP(28) * t143; -MDP(28) * t173 - t246; -pkin(4) * MDP(28) - t246; MDP(28); t192 * MDP(33) + t243; t195; t198; t197; t253; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
