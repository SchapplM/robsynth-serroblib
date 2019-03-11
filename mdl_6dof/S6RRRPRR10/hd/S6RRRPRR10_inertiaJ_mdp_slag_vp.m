% Calculate joint inertia matrix for
% S6RRRPRR10
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
%   see S6RRRPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR10_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:21:59
% EndTime: 2019-03-09 19:22:03
% DurationCPUTime: 1.15s
% Computational Cost: add. (998->208), mult. (1849->278), div. (0->0), fcn. (1862->8), ass. (0->100)
t258 = 2 * MDP(19);
t196 = cos(qJ(2));
t257 = 0.2e1 * t196;
t190 = sin(qJ(5));
t192 = sin(qJ(2));
t191 = sin(qJ(3));
t194 = cos(qJ(5));
t234 = t194 * t191;
t195 = cos(qJ(3));
t237 = t192 * t195;
t152 = t190 * t237 - t192 * t234;
t162 = t190 * t191 + t194 * t195;
t153 = t162 * t192;
t189 = sin(qJ(6));
t193 = cos(qJ(6));
t127 = t193 * t152 + t189 * t153;
t128 = -t189 * t152 + t193 * t153;
t254 = t128 * MDP(31) - t127 * MDP(32);
t256 = t153 * MDP(24) - t152 * MDP(25) + t254;
t170 = -t196 * pkin(2) - t192 * pkin(8) - pkin(1);
t241 = pkin(7) * t195;
t150 = t191 * t170 + t196 * t241;
t184 = t196 * pkin(3);
t238 = t191 * t196;
t231 = pkin(7) * t238 - t195 * t170;
t147 = t184 + t231;
t133 = t196 * pkin(4) - pkin(9) * t237 + t147;
t146 = -t196 * qJ(4) + t150;
t139 = t191 * t192 * pkin(9) + t146;
t214 = -t194 * t133 + t190 * t139;
t240 = t196 * pkin(5);
t116 = -t153 * pkin(10) - t214 + t240;
t123 = t190 * t133 + t194 * t139;
t117 = -t152 * pkin(10) + t123;
t215 = -t193 * t116 + t189 * t117;
t236 = t193 * t117;
t206 = -t215 * MDP(34) - (t189 * t116 + t236) * MDP(35);
t208 = -t214 * MDP(27) - t123 * MDP(28);
t255 = -MDP(16) * t231 - t150 * MDP(17) - t206 - t208;
t253 = -t195 * pkin(3) - t191 * qJ(4);
t243 = pkin(8) - pkin(9);
t171 = t243 * t191;
t172 = t243 * t195;
t144 = -t194 * t171 + t190 * t172;
t145 = t190 * t171 + t194 * t172;
t164 = -t190 * t195 + t234;
t130 = -t162 * pkin(10) + t145;
t137 = t193 * t162 + t189 * t164;
t138 = -t189 * t162 + t193 * t164;
t204 = -t164 * pkin(10) - t144;
t213 = t138 * MDP(31) - t137 * MDP(32) - (t189 * t130 - t193 * t204) * MDP(34) - (t193 * t130 + t189 * t204) * MDP(35);
t252 = t164 * MDP(24) - t162 * MDP(25) - t144 * MDP(27) - t145 * MDP(28) + t213;
t232 = -(t189 * t190 - t193 * t194) * MDP(34) - (t189 * t194 + t193 * t190) * MDP(35);
t251 = t194 * MDP(27) - t190 * MDP(28) + t232;
t180 = t196 * MDP(33);
t249 = t180 + t256;
t169 = -pkin(2) + t253;
t158 = t195 * pkin(4) - t169;
t142 = t162 * pkin(5) + t158;
t248 = 0.2e1 * t142;
t247 = 0.2e1 * t158;
t246 = -2 * MDP(23);
t245 = -2 * MDP(30);
t244 = -pkin(3) - pkin(4);
t242 = pkin(3) * t191;
t239 = pkin(3) * MDP(21);
t168 = t194 * qJ(4) + t190 * t244;
t235 = t193 * t168;
t233 = t195 * qJ(4);
t185 = t191 ^ 2;
t187 = t195 ^ 2;
t230 = t185 + t187;
t229 = qJ(4) * MDP(20);
t227 = t128 * MDP(29);
t167 = t190 * qJ(4) - t194 * t244;
t166 = -pkin(5) - t167;
t140 = -t193 * t166 + t189 * t168;
t226 = t140 * MDP(34);
t225 = (t189 * t166 + t235) * MDP(35);
t223 = t153 * MDP(22);
t221 = t167 * MDP(27);
t220 = t168 * MDP(28);
t219 = t195 * MDP(12);
t218 = MDP(17) - MDP(20);
t217 = MDP(26) + MDP(33);
t216 = MDP(15) + t217;
t212 = t146 * t195 + t147 * t191;
t211 = t195 * MDP(13) - t191 * MDP(14);
t209 = t191 * MDP(18) - t195 * MDP(20);
t205 = -t225 - t226;
t202 = -MDP(26) - t220 - t221;
t201 = (t193 * MDP(34) - t189 * MDP(35)) * pkin(5);
t174 = t192 * t233;
t143 = t174 + (t191 * t244 - pkin(7)) * t192;
t199 = MDP(18) + t251;
t198 = -t191 * MDP(13) + t252;
t176 = pkin(8) * t238;
t151 = -t174 + (pkin(7) + t242) * t192;
t124 = t152 * pkin(5) + t143;
t1 = [MDP(1) + (t146 ^ 2 + t147 ^ 2 + t151 ^ 2) * MDP(21) + (t152 * t246 + t223) * t153 + (t127 * t245 + t227) * t128 + t216 * t196 ^ 2 + 0.2e1 * (t152 * MDP(27) + t153 * MDP(28)) * t143 + 0.2e1 * (t127 * MDP(34) + t128 * MDP(35)) * t124 + (t147 * MDP(18) - t146 * MDP(20) + pkin(1) * MDP(9) - t255 + t256) * t257 + (-0.2e1 * pkin(1) * MDP(10) + (MDP(5) - t211) * t257 + (-t146 * t191 + t147 * t195) * t258 + 0.2e1 * t209 * t151 + (t187 * MDP(11) - 0.2e1 * t191 * t219 + MDP(4) + 0.2e1 * (t191 * MDP(16) + t195 * MDP(17)) * pkin(7)) * t192) * t192; t164 * t223 + t138 * t227 + (-t138 * t127 - t128 * t137) * MDP(30) + (-t164 * t152 - t153 * t162) * MDP(23) + t212 * MDP(19) + (pkin(8) * t212 + t151 * t169) * MDP(21) + t176 * MDP(16) + (-t151 * t195 + t176) * MDP(18) + (t124 * t137 + t142 * t127) * MDP(34) + (t124 * t138 + t142 * t128) * MDP(35) + (t143 * t162 + t158 * t152) * MDP(27) + (t143 * t164 + t158 * t153) * MDP(28) - t151 * t191 * MDP(20) + (-pkin(7) * MDP(9) + t195 * t191 * MDP(11) + MDP(6) + (-t185 + t187) * MDP(12) + (-pkin(2) * t191 - t241) * MDP(16) + (-pkin(2) * t195 + pkin(7) * t191) * MDP(17) + t209 * t169) * t192 + (-pkin(7) * MDP(10) + MDP(7) + (pkin(8) * t218 - MDP(14)) * t195 + t198) * t196; MDP(8) + t185 * MDP(11) + (pkin(8) ^ 2 * t230 + t169 ^ 2) * MDP(21) + t162 * MDP(27) * t247 + t137 * MDP(34) * t248 + t230 * pkin(8) * t258 + (MDP(22) * t164 + MDP(28) * t247 + t162 * t246) * t164 + (MDP(29) * t138 + MDP(35) * t248 + t137 * t245) * t138 + 0.2e1 * (pkin(2) * MDP(16) - t169 * MDP(18)) * t195 + 0.2e1 * (-pkin(2) * MDP(17) - t169 * MDP(20) + t219) * t191; (-0.2e1 * t184 - t231) * MDP(18) + t150 * MDP(20) + (-t147 * pkin(3) + t146 * qJ(4)) * MDP(21) + (-MDP(15) + t202 + t205 - 0.2e1 * t229) * t196 + (t253 * MDP(19) + t211) * t192 - t249 + t255; t195 * MDP(14) + (t233 - t242) * MDP(19) + ((MDP(21) * qJ(4) - t218) * t195 + (-MDP(16) - MDP(18) - t239) * t191) * pkin(8) - t198; 0.2e1 * pkin(3) * MDP(18) + 0.2e1 * t229 + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(21) + 0.2e1 * t221 + 0.2e1 * t220 + 0.2e1 * t226 + 0.2e1 * t225 + t216; MDP(19) * t237 + t147 * MDP(21) + t196 * t199; (MDP(21) * pkin(8) + MDP(19)) * t191; -t199 - t239; MDP(21); t196 * MDP(26) + (t193 * t240 - t215) * MDP(34) + (-t236 + (-t116 - t240) * t189) * MDP(35) + t208 + t249; t252; -MDP(33) + (-t193 * pkin(5) - t140) * MDP(34) + (-t235 + (pkin(5) - t166) * t189) * MDP(35) + t202; t251; 0.2e1 * t201 + t217; t180 + t206 + t254; t213; -MDP(33) + t205; t232; MDP(33) + t201; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
