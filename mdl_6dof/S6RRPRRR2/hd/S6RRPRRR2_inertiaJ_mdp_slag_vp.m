% Calculate joint inertia matrix for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR2_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:19:12
% EndTime: 2019-03-09 13:19:14
% DurationCPUTime: 0.71s
% Computational Cost: add. (1374->166), mult. (2582->230), div. (0->0), fcn. (3089->10), ass. (0->95)
t194 = sin(qJ(5));
t198 = cos(qJ(5));
t193 = sin(qJ(6));
t197 = cos(qJ(6));
t173 = t193 * t194 - t197 * t198;
t174 = t193 * t198 + t194 * t197;
t222 = t174 * MDP(29) - t173 * MDP(30);
t212 = t194 * MDP(22) + t198 * MDP(23) + t222;
t206 = MDP(25) * t198 - MDP(26) * t194;
t240 = t173 * MDP(32) + t174 * MDP(33);
t243 = t206 - t240;
t192 = cos(pkin(11));
t181 = pkin(2) * t192 + pkin(3);
t195 = sin(qJ(4));
t191 = sin(pkin(11));
t235 = pkin(2) * t191;
t236 = cos(qJ(4));
t162 = -t195 * t181 - t236 * t235;
t160 = pkin(9) - t162;
t153 = (-pkin(10) - t160) * t194;
t188 = t198 * pkin(10);
t154 = t160 * t198 + t188;
t132 = t153 * t197 - t154 * t193;
t133 = t153 * t193 + t154 * t197;
t242 = t132 * MDP(32) - t133 * MDP(33);
t177 = (-pkin(9) - pkin(10)) * t194;
t178 = pkin(9) * t198 + t188;
t155 = t177 * t197 - t178 * t193;
t156 = t177 * t193 + t178 * t197;
t241 = t155 * MDP(32) - t156 * MDP(33);
t199 = cos(qJ(2));
t184 = -t199 * pkin(2) - pkin(1);
t196 = sin(qJ(2));
t208 = t191 * t196 - t192 * t199;
t158 = t208 * pkin(3) + t184;
t239 = 0.2e1 * t158;
t238 = 0.2e1 * t199;
t237 = -2 * MDP(28);
t169 = t191 * t199 + t192 * t196;
t146 = t169 * t195 + t236 * t208;
t234 = pkin(5) * t146;
t233 = t198 * pkin(5);
t231 = -qJ(3) - pkin(7);
t175 = t231 * t196;
t176 = t231 * t199;
t151 = t192 * t175 + t176 * t191;
t136 = -pkin(8) * t169 + t151;
t152 = t191 * t175 - t192 * t176;
t137 = -t208 * pkin(8) + t152;
t123 = -t236 * t136 + t195 * t137;
t230 = t123 * t198;
t147 = t236 * t169 - t195 * t208;
t229 = t147 * t194;
t228 = t147 * t198;
t227 = t194 * t198;
t122 = t146 * pkin(4) - t147 * pkin(9) + t158;
t124 = t195 * t136 + t236 * t137;
t225 = t198 * t124;
t110 = t225 + (-pkin(10) * t147 + t122) * t194;
t226 = t197 * t110;
t128 = t173 * t147;
t219 = MDP(27) * t128;
t127 = t174 * t147;
t125 = t127 * MDP(30);
t126 = t128 * MDP(29);
t218 = t147 * MDP(19);
t161 = t236 * t181 - t195 * t235;
t217 = t161 * MDP(18);
t216 = t162 * MDP(19);
t215 = MDP(24) + MDP(31);
t214 = t146 * MDP(31) - t125 - t126;
t213 = MDP(21) * t227;
t112 = t198 * t122 - t124 * t194;
t109 = -pkin(10) * t228 + t112 + t234;
t106 = t197 * t109 - t110 * t193;
t211 = -pkin(4) * t147 - pkin(9) * t146;
t189 = t194 ^ 2;
t210 = t189 * MDP(20) + MDP(17) + 0.2e1 * t213 + (MDP(27) * t174 + t173 * t237) * t174;
t159 = -pkin(4) - t161;
t209 = -t146 * t160 + t147 * t159;
t207 = MDP(22) * t198 - MDP(23) * t194;
t205 = -MDP(25) * t194 - MDP(26) * t198;
t203 = (MDP(32) * t197 - MDP(33) * t193) * pkin(5);
t202 = 0.2e1 * t240;
t190 = t198 ^ 2;
t201 = (-t127 * t174 + t128 * t173) * MDP(28) - t174 * t219 - t123 * MDP(18) - t124 * MDP(19) + ((-t189 + t190) * MDP(21) + MDP(20) * t227 + MDP(15)) * t147 + (-MDP(16) + t212) * t146;
t183 = -pkin(4) - t233;
t157 = t159 - t233;
t119 = t123 * t194;
t116 = pkin(5) * t229 + t123;
t115 = t116 * t174;
t114 = t116 * t173;
t113 = t122 * t194 + t225;
t107 = t109 * t193 + t226;
t1 = [(t151 ^ 2 + t152 ^ 2 + t184 ^ 2) * MDP(12) + t218 * t239 + pkin(1) * MDP(9) * t238 + MDP(1) + t215 * t146 ^ 2 - (t127 * t237 - t219) * t128 + (-0.2e1 * MDP(10) * pkin(1) + MDP(4) * t196 + MDP(5) * t238) * t196 + (MDP(20) * t190 + MDP(13) - 0.2e1 * t213) * t147 ^ 2 + (MDP(18) * t239 - 0.2e1 * t126 - 0.2e1 * t125 + 0.2e1 * (-MDP(14) + t207) * t147) * t146 + 0.2e1 * (t106 * t146 + t116 * t127) * MDP(32) + 0.2e1 * (-t107 * t146 - t116 * t128) * MDP(33) + 0.2e1 * (t112 * t146 + t123 * t229) * MDP(25) + 0.2e1 * (-t113 * t146 + t123 * t228) * MDP(26) + 0.2e1 * (-t151 * t169 - t152 * t208) * MDP(11); t201 + (t209 * t194 - t230) * MDP(25) + t199 * MDP(7) + t196 * MDP(6) + (t127 * t157 + t132 * t146 + t114) * MDP(32) + (-t128 * t157 - t133 * t146 + t115) * MDP(33) + (t209 * t198 + t119) * MDP(26) + (-t199 * MDP(10) - t196 * MDP(9)) * pkin(7) + ((t151 * t192 + t152 * t191) * MDP(12) + (-t192 * t169 - t191 * t208) * MDP(11)) * pkin(2); MDP(8) + (t191 ^ 2 + t192 ^ 2) * MDP(12) * pkin(2) ^ 2 - 0.2e1 * t206 * t159 + t157 * t202 + 0.2e1 * t217 + 0.2e1 * t216 + t210; t184 * MDP(12) + t218 + (MDP(18) + t243) * t146; 0; MDP(12); t201 + (t211 * t198 + t119) * MDP(26) + (t127 * t183 + t146 * t155 + t114) * MDP(32) + (-t128 * t183 - t146 * t156 + t115) * MDP(33) + (t211 * t194 - t230) * MDP(25); t210 + t216 + t217 + t206 * (pkin(4) - t159) + t240 * (t157 + t183); 0; 0.2e1 * pkin(4) * t206 + t183 * t202 + t210; t146 * MDP(24) + t112 * MDP(25) - t113 * MDP(26) + (t197 * t234 + t106) * MDP(32) + (-t226 + (-t109 - t234) * t193) * MDP(33) + t207 * t147 + t214; t205 * t160 + t212 + t242; t243; t205 * pkin(9) + t212 + t241; 0.2e1 * t203 + t215; MDP(32) * t106 - MDP(33) * t107 + t214; t222 + t242; -t240; t222 + t241; MDP(31) + t203; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
