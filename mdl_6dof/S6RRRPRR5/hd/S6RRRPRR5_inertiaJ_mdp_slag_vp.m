% Calculate joint inertia matrix for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR5_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:23:12
% EndTime: 2019-03-09 18:23:14
% DurationCPUTime: 0.65s
% Computational Cost: add. (885->182), mult. (1575->240), div. (0->0), fcn. (1699->8), ass. (0->99)
t184 = sin(qJ(6));
t185 = sin(qJ(5));
t188 = cos(qJ(6));
t189 = cos(qJ(5));
t162 = t184 * t189 + t185 * t188;
t164 = -t184 * t185 + t188 * t189;
t219 = t164 * MDP(31) - t162 * MDP(32);
t242 = MDP(16) - MDP(19);
t241 = -MDP(17) + MDP(20);
t215 = MDP(27) * t189;
t198 = -MDP(28) * t185 + t215;
t220 = t164 * MDP(34) - t162 * MDP(35);
t240 = t198 + t220;
t239 = MDP(34) * t162 + MDP(35) * t164;
t191 = cos(qJ(2));
t238 = 0.2e1 * t191;
t237 = -2 * MDP(19);
t236 = 2 * MDP(20);
t235 = -2 * MDP(30);
t234 = pkin(3) + pkin(9);
t233 = -pkin(8) - pkin(7);
t186 = sin(qJ(3));
t187 = sin(qJ(2));
t190 = cos(qJ(3));
t165 = t186 * t191 + t187 * t190;
t232 = pkin(5) * t165;
t177 = -pkin(2) * t190 - pkin(3);
t172 = -pkin(9) + t177;
t231 = -pkin(10) + t172;
t230 = -pkin(10) - t234;
t229 = MDP(21) * pkin(3);
t163 = t186 * t187 - t190 * t191;
t228 = qJ(4) * t163;
t178 = -pkin(2) * t191 - pkin(1);
t197 = -qJ(4) * t165 + t178;
t117 = t234 * t163 + t197;
t205 = pkin(10) * t163 + t117;
t169 = t233 * t187;
t170 = t233 * t191;
t144 = -t169 * t190 - t170 * t186;
t123 = pkin(4) * t165 + t144;
t226 = t123 * t185;
t109 = t205 * t189 + t226;
t227 = t109 * t188;
t145 = t169 * t186 - t170 * t190;
t124 = -pkin(4) * t163 + t145;
t119 = t124 * t185;
t120 = t124 * t189;
t173 = pkin(2) * t186 + qJ(4);
t225 = t163 * t173;
t223 = t185 * t189;
t217 = MDP(21) * t177;
t216 = MDP(27) * t185;
t214 = MDP(28) * t189;
t128 = t162 * t163;
t213 = MDP(29) * t128;
t127 = t164 * t163;
t121 = t127 * MDP(32);
t122 = t128 * MDP(31);
t210 = t165 * MDP(25);
t180 = t189 * MDP(24);
t209 = MDP(26) + MDP(33);
t208 = 0.2e1 * t163;
t207 = t165 * MDP(33) + t121 + t122;
t206 = MDP(23) * t223;
t118 = t189 * t123;
t108 = -t205 * t185 + t118 + t232;
t105 = t188 * t108 - t109 * t184;
t150 = t231 * t185;
t151 = t231 * t189;
t129 = -t150 * t184 + t151 * t188;
t130 = t150 * t188 + t151 * t184;
t204 = t129 * MDP(34) - t130 * MDP(35) + t219;
t166 = t230 * t185;
t167 = t230 * t189;
t135 = -t166 * t184 + t167 * t188;
t136 = t166 * t188 + t167 * t184;
t203 = t135 * MDP(34) - t136 * MDP(35) + t219;
t183 = t189 ^ 2;
t202 = t183 * MDP(22) + MDP(15) - 0.2e1 * t206 + (MDP(29) * t164 + t162 * t235) * t164;
t201 = t165 * t234 + t228;
t200 = -t165 * t172 + t225;
t199 = MDP(24) * t185 + MDP(25) * t189;
t196 = t236 + 0.2e1 * t214 + 0.2e1 * t216;
t195 = (MDP(34) * t188 - MDP(35) * t184) * pkin(5);
t194 = 0.2e1 * t239;
t182 = t185 ^ 2;
t193 = t120 * MDP(28) + (t127 * t164 - t128 * t162) * MDP(30) + t164 * t213 + t241 * t145 - t242 * t144 + ((-t182 + t183) * MDP(23) + MDP(22) * t223 - MDP(14)) * t163 + (t180 + MDP(13) + t219) * t165;
t181 = t185 * pkin(5);
t175 = qJ(4) + t181;
t168 = t173 + t181;
t132 = pkin(3) * t163 + t197;
t115 = (-pkin(5) * t189 - pkin(4)) * t163 + t145;
t114 = t115 * t164;
t113 = t115 * t162;
t112 = t117 * t189 + t226;
t111 = -t117 * t185 + t118;
t106 = t108 * t184 + t227;
t1 = [(t132 ^ 2 + t144 ^ 2 + t145 ^ 2) * MDP(21) + pkin(1) * MDP(9) * t238 + MDP(1) + (t178 * MDP(16) - t132 * MDP(19)) * t208 + (t182 * MDP(22) + 0.2e1 * t206) * t163 ^ 2 + (-t127 * t235 + t213) * t128 + (-0.2e1 * MDP(10) * pkin(1) + MDP(4) * t187 + MDP(5) * t238) * t187 + (MDP(11) + t209) * t165 ^ 2 + (0.2e1 * MDP(17) * t178 - 0.2e1 * MDP(20) * t132 + 0.2e1 * t122 + 0.2e1 * t121 + (-MDP(12) + t199) * t208) * t165 + 0.2e1 * (t105 * t165 - t115 * t127) * MDP(34) + 0.2e1 * (t144 * t165 - t145 * t163) * MDP(18) + 0.2e1 * (-t106 * t165 + t115 * t128) * MDP(35) + 0.2e1 * (-t112 * t165 + t163 * t119) * MDP(28) + 0.2e1 * (t111 * t165 - t163 * t120) * MDP(27); (-MDP(10) * t191 - MDP(9) * t187) * pkin(7) + t193 + (t200 * MDP(28) - t210) * t185 + (-t127 * t168 + t129 * t165 + t113) * MDP(34) + (t128 * t168 - t130 * t165 + t114) * MDP(35) + (t165 * t177 - t225) * MDP(18) + (t144 * t177 + t145 * t173) * MDP(21) + t187 * MDP(6) + t191 * MDP(7) + (-t200 * t189 + t119) * MDP(27); MDP(8) + ((2 * MDP(19)) + t217) * t177 + t168 * t194 + 0.2e1 * (t190 * MDP(16) - t186 * MDP(17)) * pkin(2) + (MDP(21) * t173 + t196) * t173 + t202; (t201 * MDP(28) - t210) * t185 + t193 + (-pkin(3) * t144 + qJ(4) * t145) * MDP(21) + (-pkin(3) * t165 - t228) * MDP(18) + (-t127 * t175 + t135 * t165 + t113) * MDP(34) + (t128 * t175 - t136 * t165 + t114) * MDP(35) + (-t201 * t189 + t119) * MDP(27); pkin(3) * t237 + qJ(4) * t236 + (-pkin(3) * t177 + qJ(4) * t173) * MDP(21) + (t241 * t186 + t242 * t190) * pkin(2) + t202 + (t214 + t216) * (qJ(4) + t173) + t239 * (t168 + t175); t175 * t194 + (t237 + t229) * pkin(3) + (MDP(21) * qJ(4) + t196) * qJ(4) + t202; MDP(21) * t144 + (MDP(18) + t240) * t165; MDP(19) + t217; MDP(19) - t229; MDP(21); t165 * MDP(26) + t111 * MDP(27) - t112 * MDP(28) + (t188 * t232 + t105) * MDP(34) + (-t227 + (-t108 - t232) * t184) * MDP(35) + t199 * t163 + t207; -t185 * MDP(25) + t198 * t172 + t180 + t204; -t234 * t215 + t180 + (MDP(28) * t234 - MDP(25)) * t185 + t203; t240; 0.2e1 * t195 + t209; t105 * MDP(34) - t106 * MDP(35) + t207; t204; t203; t220; MDP(33) + t195; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
