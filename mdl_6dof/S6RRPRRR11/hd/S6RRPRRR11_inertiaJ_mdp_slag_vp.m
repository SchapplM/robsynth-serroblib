% Calculate joint inertia matrix for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR11_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:32:52
% EndTime: 2019-03-09 14:32:54
% DurationCPUTime: 0.70s
% Computational Cost: add. (953->182), mult. (1705->241), div. (0->0), fcn. (1773->8), ass. (0->99)
t239 = pkin(3) + pkin(7);
t180 = sin(qJ(4));
t184 = cos(qJ(4));
t186 = -pkin(2) - pkin(8);
t223 = -pkin(9) + t186;
t159 = t223 * t180;
t160 = t223 * t184;
t179 = sin(qJ(5));
t183 = cos(qJ(5));
t137 = -t159 * t179 + t183 * t160;
t138 = t159 * t183 + t160 * t179;
t156 = t179 * t184 + t180 * t183;
t157 = -t179 * t180 + t183 * t184;
t120 = -pkin(10) * t157 + t137;
t121 = -pkin(10) * t156 + t138;
t178 = sin(qJ(6));
t182 = cos(qJ(6));
t195 = t182 * t156 + t157 * t178;
t198 = -t156 * t178 + t182 * t157;
t197 = t198 * MDP(31) - t195 * MDP(32) + (t120 * t182 - t121 * t178) * MDP(34) - (t120 * t178 + t121 * t182) * MDP(35);
t190 = t157 * MDP(24) - t156 * MDP(25) + t137 * MDP(27) - t138 * MDP(28) + t197;
t238 = (MDP(20) * t186 + MDP(17)) * t184 + (-MDP(21) * t186 - MDP(18)) * t180 + t190;
t185 = cos(qJ(2));
t216 = t184 * t185;
t217 = t180 * t185;
t143 = t179 * t217 - t183 * t216;
t144 = t156 * t185;
t118 = -t182 * t143 - t144 * t178;
t119 = t143 * t178 - t144 * t182;
t237 = t119 * MDP(31) - t118 * MDP(32);
t236 = -t144 * MDP(24) + t143 * MDP(25);
t215 = MDP(34) * t198 - MDP(35) * t195;
t196 = t157 * MDP(27) - t156 * MDP(28) + t215;
t234 = 0.2e1 * t185;
t233 = 0.2e1 * MDP(27);
t232 = 0.2e1 * MDP(28);
t231 = -2 * MDP(30);
t230 = 0.2e1 * MDP(34);
t229 = 0.2e1 * MDP(35);
t228 = pkin(4) * t179;
t181 = sin(qJ(2));
t227 = pkin(4) * t181;
t226 = pkin(4) * t183;
t225 = pkin(5) * t181;
t224 = pkin(5) * t182;
t222 = MDP(14) * pkin(2);
t162 = t239 * t181;
t158 = t184 * t162;
t202 = -qJ(3) * t181 - pkin(1);
t154 = t185 * t186 + t202;
t201 = pkin(9) * t185 - t154;
t122 = t180 * t201 + t158 + t227;
t219 = t162 * t180;
t123 = -t184 * t201 + t219;
t220 = t123 * t183;
t114 = t122 * t179 + t220;
t108 = pkin(10) * t143 + t114;
t221 = t108 * t182;
t218 = t180 * t184;
t167 = t180 * pkin(4) + qJ(3);
t163 = t239 * t185;
t214 = MDP(20) * t180;
t213 = MDP(21) * t184;
t212 = MDP(22) * t157;
t211 = MDP(29) * t198;
t168 = pkin(5) + t226;
t164 = t182 * t168;
t145 = -t178 * t228 + t164;
t210 = MDP(34) * t145;
t146 = t168 * t178 + t182 * t228;
t209 = MDP(35) * t146;
t208 = MDP(35) * t178;
t207 = t183 * MDP(27);
t206 = pkin(7) ^ 2 * MDP(14);
t205 = MDP(26) + MDP(33);
t204 = t181 * MDP(33) + t237;
t147 = pkin(4) * t216 + t163;
t203 = MDP(16) * t218;
t200 = MDP(19) + t205;
t199 = MDP(12) - t222;
t113 = t183 * t122 - t123 * t179;
t107 = pkin(10) * t144 + t113 + t225;
t104 = t182 * t107 - t108 * t178;
t105 = t107 * t178 + t221;
t194 = -MDP(17) * t180 - MDP(18) * t184;
t193 = MDP(20) * t184 - MDP(21) * t180;
t192 = t181 * MDP(26) + t204 + t236;
t191 = (MDP(34) * t182 - t208) * pkin(5);
t189 = MDP(14) * pkin(7) + MDP(11) + t193;
t177 = t185 ^ 2;
t176 = t184 ^ 2;
t175 = t181 ^ 2;
t174 = t180 ^ 2;
t161 = -pkin(2) * t185 + t202;
t140 = pkin(5) * t156 + t167;
t134 = t154 * t184 + t219;
t133 = -t154 * t180 + t158;
t124 = -pkin(5) * t143 + t147;
t1 = [pkin(1) * MDP(9) * t234 + MDP(1) + (MDP(12) * t234 + MDP(14) * t161) * t161 + (t174 * MDP(15) + 0.2e1 * t203 + t206) * t177 - (-MDP(22) * t144 + 0.2e1 * t143 * MDP(23)) * t144 + (MDP(29) * t119 + t118 * t231) * t119 + (MDP(4) + t200 + t206) * t175 + 0.2e1 * (-pkin(1) * MDP(10) - t161 * MDP(13) + (MDP(5) + t194) * t185 + t236 + t237) * t181 + (t104 * t181 + t118 * t124) * t230 + (-t105 * t181 + t119 * t124) * t229 + (t113 * t181 - t143 * t147) * t233 + (-t114 * t181 - t144 * t147) * t232 + 0.2e1 * (-t134 * t181 - t163 * t217) * MDP(21) + 0.2e1 * (t133 * t181 + t163 * t216) * MDP(20) + 0.2e1 * (t175 + t177) * MDP(11) * pkin(7); (-t118 * t198 - t119 * t195) * MDP(30) + t119 * t211 - t144 * t212 + (t143 * t157 + t144 * t156) * MDP(23) + (-t143 * t167 + t147 * t156) * MDP(27) + (-t144 * t167 + t147 * t157) * MDP(28) + (t118 * t140 + t124 * t195) * MDP(34) + (t119 * t140 + t124 * t198) * MDP(35) + (t213 + t214) * t163 + (-MDP(15) * t218 + MDP(7) + (t174 - t176) * MDP(16) + (-MDP(10) + MDP(13)) * pkin(7) + t189 * qJ(3)) * t185 + (-pkin(2) * MDP(11) + MDP(6) + (-MDP(9) + t199) * pkin(7) + t238) * t181; -0.2e1 * t203 + t167 * t156 * t233 + t140 * t195 * t230 + t176 * MDP(15) + MDP(8) + (-0.2e1 * MDP(12) + t222) * pkin(2) + (-0.2e1 * MDP(23) * t156 + t167 * t232 + t212) * t157 + (t140 * t229 + t195 * t231 + t211) * t198 + (MDP(14) * qJ(3) + 0.2e1 * MDP(13) + 0.2e1 * t213 + 0.2e1 * t214) * qJ(3); (t189 + t196) * t181; t199; MDP(14); t181 * MDP(19) + t133 * MDP(20) - t134 * MDP(21) + (t181 * t226 + t113) * MDP(27) + (-t220 + (-t122 - t227) * t179) * MDP(28) + (t145 * t181 + t104) * MDP(34) + (-t146 * t181 - t105) * MDP(35) + t194 * t185 + t192; t238; t193 + t196; 0.2e1 * (-MDP(28) * t179 + t207) * pkin(4) + 0.2e1 * t210 - 0.2e1 * t209 + t200; t113 * MDP(27) - t114 * MDP(28) + (t181 * t224 + t104) * MDP(34) + (-t221 + (-t107 - t225) * t178) * MDP(35) + t192; t190; t196; (t164 + t224) * MDP(34) + (-pkin(5) - t168) * t208 + (t207 + (-MDP(34) * t178 - MDP(35) * t182 - MDP(28)) * t179) * pkin(4) + t205; 0.2e1 * t191 + t205; t104 * MDP(34) - t105 * MDP(35) + t204; t197; t215; MDP(33) - t209 + t210; MDP(33) + t191; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
