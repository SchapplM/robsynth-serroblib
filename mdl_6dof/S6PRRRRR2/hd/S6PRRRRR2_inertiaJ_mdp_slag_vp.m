% Calculate joint inertia matrix for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR2_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:45:25
% EndTime: 2019-03-09 00:45:27
% DurationCPUTime: 0.70s
% Computational Cost: add. (804->168), mult. (1641->237), div. (0->0), fcn. (1887->12), ass. (0->96)
t198 = sin(qJ(5));
t203 = cos(qJ(5));
t197 = sin(qJ(6));
t202 = cos(qJ(6));
t173 = t197 * t198 - t202 * t203;
t175 = t197 * t203 + t202 * t198;
t229 = t175 * MDP(28) - t173 * MDP(29);
t217 = t198 * MDP(21) + t203 * MDP(22) + t229;
t212 = t203 * MDP(24) - t198 * MDP(25);
t246 = t173 * MDP(31) + t175 * MDP(32);
t199 = sin(qJ(4));
t185 = t199 * pkin(3) + pkin(10);
t169 = (-pkin(11) - t185) * t198;
t192 = t203 * pkin(11);
t170 = t203 * t185 + t192;
t142 = t202 * t169 - t197 * t170;
t143 = t197 * t169 + t202 * t170;
t248 = t142 * MDP(31) - t143 * MDP(32);
t178 = (-pkin(10) - pkin(11)) * t198;
t180 = t203 * pkin(10) + t192;
t154 = t202 * t178 - t197 * t180;
t156 = t197 * t178 + t202 * t180;
t247 = t154 * MDP(31) - t156 * MDP(32);
t204 = cos(qJ(3));
t188 = -t204 * pkin(3) - pkin(2);
t245 = 0.2e1 * t188;
t244 = -2 * MDP(27);
t243 = pkin(8) + pkin(9);
t242 = cos(qJ(4));
t200 = sin(qJ(3));
t174 = t199 * t200 - t242 * t204;
t241 = t174 * pkin(5);
t240 = t203 * pkin(5);
t179 = t243 * t200;
t181 = t243 * t204;
t155 = t242 * t179 + t199 * t181;
t238 = t155 * t203;
t176 = t199 * t204 + t242 * t200;
t237 = t176 * t198;
t236 = t176 * t203;
t195 = sin(pkin(6));
t201 = sin(qJ(2));
t235 = t195 * t201;
t205 = cos(qJ(2));
t234 = t195 * t205;
t233 = t198 * t203;
t145 = t174 * pkin(4) - t176 * pkin(10) + t188;
t157 = -t199 * t179 + t242 * t181;
t231 = t203 * t157;
t119 = t231 + (-pkin(11) * t176 + t145) * t198;
t232 = t202 * t119;
t196 = cos(pkin(6));
t161 = t196 * t204 - t200 * t235;
t162 = t196 * t200 + t204 * t235;
t136 = t199 * t161 + t242 * t162;
t127 = -t198 * t136 - t203 * t234;
t128 = t203 * t136 - t198 * t234;
t114 = t202 * t127 - t197 * t128;
t115 = t197 * t127 + t202 * t128;
t230 = t114 * MDP(31) - t115 * MDP(32);
t227 = MDP(10) * t204;
t226 = MDP(18) * t176;
t138 = t173 * t176;
t225 = MDP(26) * t138;
t137 = t175 * t176;
t133 = t137 * MDP(29);
t134 = t138 * MDP(28);
t220 = MDP(23) + MDP(30);
t219 = t174 * MDP(30) - t133 - t134;
t218 = MDP(20) * t233;
t120 = t203 * t145 - t198 * t157;
t118 = -pkin(11) * t236 + t120 + t241;
t110 = t202 * t118 - t197 * t119;
t186 = -t242 * pkin(3) - pkin(4);
t216 = -pkin(4) * t176 - pkin(10) * t174;
t193 = t198 ^ 2;
t215 = t193 * MDP(19) + MDP(16) + 0.2e1 * t218 + (MDP(26) * t175 + t173 * t244) * t175;
t214 = -t174 * t185 + t176 * t186;
t213 = t203 * MDP(21) - t198 * MDP(22);
t211 = -MDP(24) * t198 - MDP(25) * t203;
t210 = (MDP(31) * t202 - MDP(32) * t197) * pkin(5);
t209 = 0.2e1 * t246;
t208 = (t242 * MDP(17) - t199 * MDP(18)) * pkin(3);
t135 = -t242 * t161 + t199 * t162;
t207 = -t136 * MDP(18) + (-MDP(17) + t246 - t212) * t135;
t194 = t203 ^ 2;
t206 = (-t175 * t137 + t138 * t173) * MDP(27) - t175 * t225 - t155 * MDP(17) - t157 * MDP(18) + ((-t193 + t194) * MDP(20) + MDP(19) * t233 + MDP(14)) * t176 + (-MDP(15) + t217) * t174;
t187 = -pkin(4) - t240;
t177 = t186 - t240;
t146 = t155 * t198;
t130 = pkin(5) * t237 + t155;
t123 = t130 * t175;
t122 = t130 * t173;
t121 = t198 * t145 + t231;
t111 = t197 * t118 + t232;
t1 = [MDP(1); (t127 * t174 + t135 * t237) * MDP(24) + (-t128 * t174 + t135 * t236) * MDP(25) + (t114 * t174 + t135 * t137) * MDP(31) + (-t115 * t174 - t135 * t138) * MDP(32) + (-t201 * MDP(4) + (-MDP(11) * t200 - MDP(17) * t174 + MDP(3) - t226 + t227) * t205) * t195; 0.2e1 * pkin(2) * t227 + t226 * t245 + MDP(2) + t220 * t174 ^ 2 - (t137 * t244 - t225) * t138 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t200 + 0.2e1 * t204 * MDP(6)) * t200 + (t194 * MDP(19) + MDP(12) - 0.2e1 * t218) * t176 ^ 2 + (MDP(17) * t245 - 0.2e1 * t134 - 0.2e1 * t133 + 0.2e1 * (-MDP(13) + t213) * t176) * t174 + 0.2e1 * (t120 * t174 + t155 * t237) * MDP(24) + 0.2e1 * (-t121 * t174 + t155 * t236) * MDP(25) + 0.2e1 * (t110 * t174 + t130 * t137) * MDP(31) + 0.2e1 * (-t111 * t174 - t130 * t138) * MDP(32); t161 * MDP(10) - t162 * MDP(11) + t207; (t177 * t137 + t142 * t174 + t122) * MDP(31) + (-t177 * t138 - t143 * t174 + t123) * MDP(32) + t200 * MDP(7) + t204 * MDP(8) + t206 + (-t200 * MDP(10) - t204 * MDP(11)) * pkin(8) + (t214 * t198 - t238) * MDP(24) + (t214 * t203 + t146) * MDP(25); t177 * t209 - 0.2e1 * t186 * t212 + MDP(9) + 0.2e1 * t208 + t215; t207; (t187 * t137 + t154 * t174 + t122) * MDP(31) + (-t187 * t138 - t156 * t174 + t123) * MDP(32) + t206 + (t216 * t198 - t238) * MDP(24) + (t216 * t203 + t146) * MDP(25); t208 + t215 + t212 * (pkin(4) - t186) + t246 * (t177 + t187); 0.2e1 * pkin(4) * t212 + t187 * t209 + t215; t127 * MDP(24) - t128 * MDP(25) + t230; t174 * MDP(23) + t120 * MDP(24) - t121 * MDP(25) + (t202 * t241 + t110) * MDP(31) + (-t232 + (-t118 - t241) * t197) * MDP(32) + t213 * t176 + t219; t211 * t185 + t217 + t248; t211 * pkin(10) + t217 + t247; 0.2e1 * t210 + t220; t230; t110 * MDP(31) - t111 * MDP(32) + t219; t229 + t248; t229 + t247; MDP(30) + t210; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
