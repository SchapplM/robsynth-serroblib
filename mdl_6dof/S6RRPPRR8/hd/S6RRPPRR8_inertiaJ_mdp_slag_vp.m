% Calculate joint inertia matrix for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR8_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:26:09
% EndTime: 2019-03-09 09:26:11
% DurationCPUTime: 0.68s
% Computational Cost: add. (787->176), mult. (1542->233), div. (0->0), fcn. (1579->8), ass. (0->86)
t166 = sin(qJ(2));
t169 = cos(qJ(2));
t143 = -pkin(2) * t169 - qJ(3) * t166 - pkin(1);
t162 = sin(pkin(10));
t207 = pkin(7) * t169;
t149 = t162 * t207;
t157 = t169 * pkin(3);
t163 = cos(pkin(10));
t206 = pkin(8) * t166;
t113 = pkin(4) * t169 + t149 + t157 + (-t143 - t206) * t163;
t126 = t162 * t143 + t163 * t207;
t123 = -qJ(4) * t169 + t126;
t118 = t162 * t206 + t123;
t165 = sin(qJ(5));
t168 = cos(qJ(5));
t103 = t168 * t113 - t118 * t165;
t104 = t113 * t165 + t118 * t168;
t201 = t163 * t166;
t202 = t162 * t168;
t129 = t165 * t201 - t166 * t202;
t137 = t162 * t165 + t163 * t168;
t130 = t137 * t166;
t220 = t130 * MDP(21) - t129 * MDP(22) + t103 * MDP(24) - t104 * MDP(25);
t218 = MDP(14) * pkin(7);
t164 = sin(qJ(6));
t167 = cos(qJ(6));
t108 = t167 * t129 + t130 * t164;
t109 = -t129 * t164 + t130 * t167;
t217 = t109 * MDP(28) - t108 * MDP(29);
t205 = -pkin(8) + qJ(3);
t144 = t205 * t162;
t145 = t205 * t163;
t120 = t168 * t144 - t145 * t165;
t121 = t144 * t165 + t145 * t168;
t138 = -t163 * t165 + t202;
t110 = -pkin(9) * t138 + t120;
t111 = -pkin(9) * t137 + t121;
t116 = t167 * t137 + t138 * t164;
t117 = -t137 * t164 + t138 * t167;
t183 = t117 * MDP(28) - t116 * MDP(29) + (t110 * t167 - t111 * t164) * MDP(31) - (t110 * t164 + t111 * t167) * MDP(32);
t216 = t138 * MDP(21) - t137 * MDP(22) + t120 * MDP(24) - t121 * MDP(25) + t183;
t200 = (-t164 * t165 + t167 * t168) * MDP(31) - (t164 * t168 + t165 * t167) * MDP(32);
t215 = MDP(24) * t168 - MDP(25) * t165 + t200;
t188 = MDP(12) - MDP(17);
t189 = MDP(11) + MDP(15);
t214 = t189 * t162 + t188 * t163;
t213 = MDP(14) + MDP(18);
t141 = -t163 * pkin(3) - t162 * qJ(4) - pkin(2);
t132 = t163 * pkin(4) - t141;
t119 = pkin(5) * t137 + t132;
t212 = 0.2e1 * t119;
t211 = 0.2e1 * t132;
t210 = -2 * MDP(20);
t209 = -2 * MDP(27);
t208 = pkin(5) * t169;
t204 = pkin(2) * MDP(14);
t98 = -pkin(9) * t129 + t104;
t203 = t167 * t98;
t198 = MDP(11) * t163;
t197 = MDP(12) * t162;
t196 = MDP(15) * t163;
t195 = MDP(17) * t162;
t194 = MDP(18) * t141;
t193 = MDP(19) * t138;
t192 = MDP(24) * t137;
t191 = MDP(26) * t117;
t190 = MDP(31) * t116;
t187 = MDP(23) + MDP(30);
t186 = t169 * MDP(30) + t217;
t97 = -pkin(9) * t130 + t103 + t208;
t94 = -t164 * t98 + t167 * t97;
t125 = t143 * t163 - t149;
t182 = MDP(31) * t94 - MDP(32) * (t164 * t97 + t203);
t124 = -t125 + t157;
t181 = t123 * t163 + t124 * t162;
t180 = -t125 * t162 + t126 * t163;
t179 = MDP(11) * t162 + MDP(12) * t163;
t178 = MDP(15) * t162 - MDP(17) * t163;
t175 = MDP(24) * t129 + MDP(25) * t130;
t173 = MDP(31) * t108 + MDP(32) * t109;
t172 = (MDP(31) * t167 - MDP(32) * t164) * pkin(5);
t148 = qJ(4) * t201;
t122 = t148 + (-pkin(7) + (-pkin(3) - pkin(4)) * t162) * t166;
t128 = -t148 + (pkin(3) * t162 + pkin(7)) * t166;
t105 = pkin(5) * t129 + t122;
t1 = [-0.2e1 * pkin(1) * t166 * MDP(10) + MDP(1) + (t125 ^ 2 + t126 ^ 2) * MDP(14) + (t123 ^ 2 + t124 ^ 2 + t128 ^ 2) * MDP(18) + t187 * t169 ^ 2 + (MDP(19) * t130 + t129 * t210) * t130 + (MDP(26) * t109 + t108 * t209) * t109 + 0.2e1 * t175 * t122 + 0.2e1 * t173 * t105 + 0.2e1 * ((-t123 * t162 + t124 * t163) * MDP(16) + (-t125 * t163 - t126 * t162) * MDP(13) + t178 * t128) * t166 + (MDP(4) + (0.2e1 * t179 + t218) * pkin(7)) * t166 ^ 2 + 0.2e1 * (-MDP(11) * t125 + MDP(12) * t126 + MDP(15) * t124 - MDP(17) * t123 + MDP(5) * t166 + MDP(9) * pkin(1) + t182 + t217 + t220) * t169; t180 * MDP(13) + t181 * MDP(16) + (-t129 * t138 - t130 * t137) * MDP(20) + (t122 * t137 + t129 * t132) * MDP(24) + (t122 * t138 + t130 * t132) * MDP(25) + (-t108 * t117 - t109 * t116) * MDP(27) + (t105 * t116 + t108 * t119) * MDP(31) + (t105 * t117 + t109 * t119) * MDP(32) + t130 * t193 + t109 * t191 + (t194 - t195 - t196) * t128 + (MDP(6) + t178 * t141 - t179 * pkin(2) + (-MDP(9) + t197 - t198 - t204) * pkin(7)) * t166 + (-pkin(7) * MDP(10) + MDP(7) + t216) * t169 + (t180 * MDP(14) + t181 * MDP(18) + t214 * t169) * qJ(3); MDP(8) + t192 * t211 + t190 * t212 + (t194 - 0.2e1 * t195 - 0.2e1 * t196) * t141 + (-0.2e1 * t197 + 0.2e1 * t198 + t204) * pkin(2) + (MDP(25) * t211 + t137 * t210 + t193) * t138 + (MDP(32) * t212 + t116 * t209 + t191) * t117 + (t213 * qJ(3) + 0.2e1 * MDP(13) + 0.2e1 * MDP(16)) * (t162 ^ 2 + t163 ^ 2) * qJ(3); MDP(18) * t128 + (t214 + t218) * t166 - t173 - t175; -MDP(25) * t138 - MDP(32) * t117 + t188 * t162 - t189 * t163 - t190 - t192 + t194 - t204; t213; MDP(16) * t201 + MDP(18) * t124 + (MDP(15) + t215) * t169; (MDP(18) * qJ(3) + MDP(16)) * t162; 0; MDP(18); t169 * MDP(23) + (t167 * t208 + t94) * MDP(31) + (-t203 + (-t97 - t208) * t164) * MDP(32) + t186 + t220; t216; 0; t215; 0.2e1 * t172 + t187; t182 + t186; t183; 0; t200; MDP(30) + t172; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
