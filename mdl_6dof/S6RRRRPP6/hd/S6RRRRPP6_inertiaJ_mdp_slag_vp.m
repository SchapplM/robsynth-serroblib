% Calculate joint inertia matrix for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP6_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:15:07
% EndTime: 2019-03-09 21:15:10
% DurationCPUTime: 1.33s
% Computational Cost: add. (1262->250), mult. (2292->329), div. (0->0), fcn. (2236->6), ass. (0->88)
t163 = sin(qJ(3));
t166 = cos(qJ(3));
t206 = -pkin(9) - pkin(8);
t141 = t206 * t163;
t142 = t206 * t166;
t162 = sin(qJ(4));
t165 = cos(qJ(4));
t123 = -t165 * t141 - t162 * t142;
t137 = t162 * t166 + t165 * t163;
t110 = t137 * pkin(5) + t123;
t124 = t162 * t141 - t165 * t142;
t192 = t165 * t166;
t136 = t162 * t163 - t192;
t111 = -t136 * pkin(5) + t124;
t170 = t137 * MDP(20) - t136 * MDP(21) + t111 * MDP(30) - t110 * MDP(31) - (MDP(23) - MDP(26)) * t123 - (MDP(24) - MDP(27)) * t124;
t223 = t170 - (t163 * MDP(16) + t166 * MDP(17)) * pkin(8) + t163 * MDP(13) + t166 * MDP(14);
t160 = pkin(4) + qJ(6);
t207 = 0.2e1 * MDP(31);
t220 = -0.2e1 * pkin(4) * MDP(26) + t160 * t207 + MDP(22);
t164 = sin(qJ(2));
t167 = cos(qJ(2));
t140 = -t167 * pkin(2) - t164 * pkin(8) - pkin(1);
t135 = t166 * t140;
t202 = pkin(9) * t164;
t205 = pkin(7) * t163;
t115 = -t166 * t202 + t135 + (-pkin(3) - t205) * t167;
t203 = pkin(7) * t167;
t181 = t166 * t203;
t118 = t181 + (t140 - t202) * t163;
t105 = t162 * t115 + t165 * t118;
t155 = t167 * pkin(4);
t189 = -t165 * t115 + t162 * t118;
t102 = t155 + t189;
t194 = t163 * t164;
t132 = -t162 * t194 + t164 * t192;
t178 = -t132 * pkin(5) - t102;
t131 = t137 * t164;
t216 = -t132 * MDP(20) + t131 * MDP(21);
t219 = -t189 * MDP(23) - t105 * MDP(24) + MDP(31) * t178 - t216;
t182 = MDP(27) + MDP(30);
t218 = 0.2e1 * t182;
t186 = t162 * MDP(24);
t217 = (t165 * MDP(23) - t186) * pkin(3);
t200 = t162 * pkin(3);
t147 = qJ(5) + t200;
t212 = t147 ^ 2;
t211 = -2 * MDP(19);
t210 = 0.2e1 * MDP(24);
t209 = 2 * MDP(25);
t208 = 2 * MDP(29);
t204 = pkin(7) * t166;
t201 = t131 * pkin(5);
t199 = qJ(5) * t131;
t198 = qJ(5) * t136;
t197 = t147 * qJ(5);
t196 = t147 * t131;
t195 = t147 * t136;
t193 = t163 * t166;
t191 = t167 * qJ(5);
t139 = pkin(3) * t194 + t164 * pkin(7);
t188 = t132 * MDP(18);
t151 = -t165 * pkin(3) - pkin(4);
t187 = t151 * MDP(26);
t185 = MDP(15) + MDP(22);
t184 = MDP(25) + MDP(29);
t183 = MDP(26) - MDP(31);
t152 = -t166 * pkin(3) - pkin(2);
t180 = MDP(12) * t193;
t177 = t105 - t201;
t101 = -t105 + t191;
t176 = -t132 * qJ(5) + t139;
t175 = t166 * MDP(13) - t163 * MDP(14);
t171 = -t137 * qJ(5) + t152;
t168 = qJ(5) ^ 2;
t158 = t166 ^ 2;
t157 = t164 ^ 2;
t156 = t163 ^ 2;
t153 = -0.2e1 * t191;
t144 = qJ(6) - t151;
t127 = t163 * t140 + t181;
t126 = -t163 * t203 + t135;
t114 = t136 * pkin(4) + t171;
t107 = t136 * t160 + t171;
t106 = t131 * pkin(4) + t176;
t103 = t131 * t160 + t176;
t100 = -t101 - t201;
t99 = t167 * qJ(6) - t178;
t1 = [(t100 ^ 2 + t103 ^ 2 + t99 ^ 2) * MDP(32) + (t101 ^ 2 + t102 ^ 2 + t106 ^ 2) * MDP(28) + MDP(1) - 0.2e1 * pkin(1) * t164 * MDP(10) + t185 * t167 ^ 2 + (t131 * t211 + t188) * t132 + (t158 * MDP(11) + MDP(4) - 0.2e1 * t180) * t157 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t175) * t164 + t216) * t167 + (t101 * t131 + t102 * t132) * t209 + (-t100 * t131 + t99 * t132) * t208 + (t103 * t131 + t99 * t167) * t207 + 0.2e1 * (-t102 * t167 - t106 * t131) * MDP(26) + 0.2e1 * (-t100 * t167 - t103 * t132) * MDP(30) + 0.2e1 * (t101 * t167 - t106 * t132) * MDP(27) + 0.2e1 * (t139 * t131 + t167 * t189) * MDP(23) + (t105 * t167 + t139 * t132) * t210 + 0.2e1 * (-t126 * t167 + t157 * t205) * MDP(16) + 0.2e1 * (t127 * t167 + t157 * t204) * MDP(17); t137 * t188 + (t100 * t111 + t103 * t107 + t99 * t110) * MDP(32) + (-t137 * t131 - t132 * t136) * MDP(19) + (t152 * t131 + t139 * t136) * MDP(23) + (t152 * t132 + t139 * t137) * MDP(24) + (t101 * t136 + t102 * t137 + t123 * t132 - t124 * t131) * MDP(25) + (-t106 * t136 - t114 * t131) * MDP(26) + (-t106 * t137 - t114 * t132) * MDP(27) + (-t101 * t124 + t102 * t123 + t106 * t114) * MDP(28) + (-t100 * t136 + t110 * t132 - t111 * t131 + t99 * t137) * MDP(29) + (-t103 * t137 - t107 * t132) * MDP(30) + (t103 * t136 + t107 * t131) * MDP(31) + (-pkin(7) * MDP(10) + MDP(7) - t223) * t167 + (-pkin(7) * MDP(9) + MDP(6) + (-t156 + t158) * MDP(12) + (-pkin(2) * t163 - t204) * MDP(16) + (-pkin(2) * t166 + t205) * MDP(17) + MDP(11) * t193) * t164; MDP(8) + t156 * MDP(11) + 0.2e1 * t180 + (t114 ^ 2 + t123 ^ 2 + t124 ^ 2) * MDP(28) + (t107 ^ 2 + t110 ^ 2 + t111 ^ 2) * MDP(32) + (MDP(18) * t137 - 0.2e1 * t114 * MDP(27) - 0.2e1 * t107 * MDP(30) + t136 * t211 + t152 * t210) * t137 + (t123 * t137 - t124 * t136) * t209 + (t110 * t137 - t111 * t136) * t208 + 0.2e1 * (t166 * MDP(16) - t163 * MDP(17)) * pkin(2) + 0.2e1 * (t152 * MDP(23) - t114 * MDP(26) + t107 * MDP(31)) * t136; t126 * MDP(16) - t127 * MDP(17) + (t151 * t132 - t196) * MDP(25) + t102 * MDP(26) + t105 * MDP(27) + (-t101 * t147 + t102 * t151) * MDP(28) + (-t144 * t132 - t196) * MDP(29) + t177 * MDP(30) + (t100 * t147 - t99 * t144) * MDP(32) + t175 * t164 + (-t187 + (-qJ(6) - t144) * MDP(31) - t217 - t185 + t182 * (-qJ(5) - t147)) * t167 + t219; (t151 * t137 - t195) * MDP(25) + (t123 * t151 + t124 * t147) * MDP(28) + (-t144 * t137 - t195) * MDP(29) + (-t110 * t144 + t111 * t147) * MDP(32) + t223; (t151 ^ 2 + t212) * MDP(28) + (t144 ^ 2 + t212) * MDP(32) + 0.2e1 * t217 + 0.2e1 * t187 + t144 * t207 + t147 * t218 + t185; (-pkin(4) * t132 - t199) * MDP(25) + (0.2e1 * t155 + t189) * MDP(26) + (t153 + t105) * MDP(27) + (-t102 * pkin(4) - t101 * qJ(5)) * MDP(28) + (-t160 * t132 - t199) * MDP(29) + (t153 + t177) * MDP(30) + (t100 * qJ(5) - t99 * t160) * MDP(32) + (-MDP(22) + (-qJ(6) - t160) * MDP(31)) * t167 + t219; (-pkin(4) * t137 - t198) * MDP(25) + (-t123 * pkin(4) + t124 * qJ(5)) * MDP(28) + (-t160 * t137 - t198) * MDP(29) + (t111 * qJ(5) - t110 * t160) * MDP(32) + t170; (-t151 * pkin(4) + t197) * MDP(28) + (t144 * t160 + t197) * MDP(32) + t182 * (0.2e1 * qJ(5) + t200) + (-t186 + (MDP(23) - t183) * t165) * pkin(3) + t220; (pkin(4) ^ 2 + t168) * MDP(28) + (t160 ^ 2 + t168) * MDP(32) + qJ(5) * t218 + t220; t102 * MDP(28) + t99 * MDP(32) + t132 * t184 - t167 * t183; t123 * MDP(28) + t110 * MDP(32) + t137 * t184; t151 * MDP(28) - t144 * MDP(32) + t183; -pkin(4) * MDP(28) - t160 * MDP(32) + t183; MDP(28) + MDP(32); -t131 * MDP(29) - t167 * MDP(30) + t100 * MDP(32); -t136 * MDP(29) + t111 * MDP(32); t147 * MDP(32) + MDP(30); MDP(32) * qJ(5) + MDP(30); 0; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
