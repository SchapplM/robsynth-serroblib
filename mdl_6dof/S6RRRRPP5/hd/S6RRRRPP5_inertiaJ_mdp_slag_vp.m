% Calculate joint inertia matrix for
% S6RRRRPP5
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
%   see S6RRRRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP5_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:09:10
% EndTime: 2019-03-09 21:09:14
% DurationCPUTime: 1.31s
% Computational Cost: add. (1258->252), mult. (2292->330), div. (0->0), fcn. (2244->6), ass. (0->89)
t163 = sin(qJ(3));
t166 = cos(qJ(3));
t206 = -pkin(9) - pkin(8);
t141 = t206 * t166;
t162 = sin(qJ(4));
t165 = cos(qJ(4));
t183 = t206 * t163;
t123 = -t141 * t162 - t165 * t183;
t137 = t162 * t166 + t163 * t165;
t110 = -qJ(6) * t137 + t123;
t124 = -t165 * t141 + t162 * t183;
t195 = t165 * t166;
t136 = t162 * t163 - t195;
t111 = t136 * qJ(6) + t124;
t173 = t137 * MDP(20) - t136 * MDP(21) - t110 * MDP(29) + t111 * MDP(30) - (MDP(23) + MDP(25)) * t123 - (MDP(24) - MDP(27)) * t124;
t223 = t173 - (MDP(16) * t163 + MDP(17) * t166) * pkin(8) + t163 * MDP(13) + t166 * MDP(14);
t164 = sin(qJ(2));
t167 = cos(qJ(2));
t140 = -pkin(2) * t167 - pkin(8) * t164 - pkin(1);
t135 = t166 * t140;
t202 = pkin(9) * t164;
t205 = pkin(7) * t163;
t115 = -t166 * t202 + t135 + (-pkin(3) - t205) * t167;
t203 = pkin(7) * t167;
t184 = t166 * t203;
t118 = t184 + (t140 - t202) * t163;
t105 = t162 * t115 + t165 * t118;
t157 = t167 * pkin(4);
t192 = -t165 * t115 + t162 * t118;
t102 = t157 + t192;
t197 = t163 * t164;
t132 = -t162 * t197 + t164 * t195;
t180 = t132 * qJ(6) - t102;
t131 = t137 * t164;
t217 = -t132 * MDP(20) + t131 * MDP(21);
t220 = -t192 * MDP(23) - t105 * MDP(24) + t180 * MDP(29) - t217;
t185 = MDP(27) + MDP(30);
t219 = 0.2e1 * t185;
t189 = t162 * MDP(24);
t218 = (t165 * MDP(23) - t189) * pkin(3);
t213 = -2 * MDP(19);
t212 = 0.2e1 * MDP(24);
t211 = 0.2e1 * MDP(25);
t210 = 2 * MDP(26);
t209 = 0.2e1 * MDP(29);
t208 = 0.2e1 * MDP(30);
t207 = 2 * MDP(31);
t168 = pkin(4) + pkin(5);
t204 = pkin(7) * t166;
t201 = qJ(5) * t131;
t200 = qJ(5) * t136;
t155 = t162 * pkin(3);
t148 = t155 + qJ(5);
t199 = t131 * t148;
t198 = t136 * t148;
t196 = t163 * t166;
t194 = t167 * qJ(5);
t139 = pkin(3) * t197 + t164 * pkin(7);
t191 = t132 * MDP(18);
t151 = pkin(3) * t165 + pkin(4);
t190 = t151 * MDP(25);
t188 = MDP(15) + MDP(22);
t187 = MDP(25) + MDP(29);
t186 = MDP(26) - MDP(31);
t182 = -0.2e1 * t194 + t105;
t153 = -pkin(3) * t166 - pkin(2);
t181 = MDP(12) * t196;
t101 = t105 - t194;
t179 = qJ(5) * t132 - t139;
t178 = t166 * MDP(13) - t163 * MDP(14);
t174 = qJ(5) * t137 - t153;
t171 = 0.2e1 * pkin(4);
t170 = qJ(5) ^ 2;
t160 = t166 ^ 2;
t159 = t164 ^ 2;
t158 = t163 ^ 2;
t146 = pkin(5) + t151;
t145 = t148 ^ 2;
t143 = t148 * qJ(5);
t129 = t131 * qJ(6);
t126 = t140 * t163 + t184;
t125 = -t163 * t203 + t135;
t114 = pkin(4) * t136 - t174;
t107 = -t168 * t136 + t174;
t106 = pkin(4) * t131 - t179;
t103 = -t168 * t131 + t179;
t100 = t101 + t129;
t99 = pkin(5) * t167 - t180;
t1 = [(t101 ^ 2 + t102 ^ 2 + t106 ^ 2) * MDP(28) + (t100 ^ 2 + t103 ^ 2 + t99 ^ 2) * MDP(32) - 0.2e1 * pkin(1) * t164 * MDP(10) + MDP(1) + t188 * t167 ^ 2 + (t131 * t213 + t191) * t132 + (t160 * MDP(11) + MDP(4) - 0.2e1 * t181) * t159 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t178) * t164 + t217) * t167 + 0.2e1 * (t126 * t167 + t159 * t204) * MDP(17) + 0.2e1 * (-t125 * t167 + t159 * t205) * MDP(16) + 0.2e1 * (-t101 * t167 - t106 * t132) * MDP(27) + 0.2e1 * (t131 * t139 + t167 * t192) * MDP(23) + (t105 * t167 + t132 * t139) * t212 + (t102 * t167 + t106 * t131) * t211 + (-t103 * t131 + t167 * t99) * t209 + (-t100 * t167 + t103 * t132) * t208 + (-t101 * t131 + t102 * t132) * t210 + (t100 * t131 - t132 * t99) * t207; (-t131 * t137 - t132 * t136) * MDP(19) + (t131 * t153 + t136 * t139) * MDP(23) + (t132 * t153 + t137 * t139) * MDP(24) + (t106 * t136 + t114 * t131) * MDP(25) + (-t101 * t136 + t102 * t137 + t123 * t132 - t124 * t131) * MDP(26) + (-t106 * t137 - t114 * t132) * MDP(27) + (t101 * t124 + t102 * t123 + t106 * t114) * MDP(28) + (-t103 * t136 - t107 * t131) * MDP(29) + (t103 * t137 + t107 * t132) * MDP(30) + (t100 * t136 - t110 * t132 + t111 * t131 - t137 * t99) * MDP(31) + (t100 * t111 + t103 * t107 + t110 * t99) * MDP(32) + t137 * t191 + (-pkin(7) * MDP(10) + MDP(7) - t223) * t167 + (MDP(6) + (-t158 + t160) * MDP(12) + (-pkin(2) * t163 - t204) * MDP(16) + (-pkin(2) * t166 + t205) * MDP(17) + MDP(11) * t196 - pkin(7) * MDP(9)) * t164; MDP(8) + t158 * MDP(11) + 0.2e1 * t181 + (t114 ^ 2 + t123 ^ 2 + t124 ^ 2) * MDP(28) + (t107 ^ 2 + t110 ^ 2 + t111 ^ 2) * MDP(32) + (MDP(18) * t137 - 0.2e1 * MDP(27) * t114 + t107 * t208 + t136 * t213 + t153 * t212) * t137 + (t123 * t137 - t124 * t136) * t210 + (-t110 * t137 + t111 * t136) * t207 + 0.2e1 * (t166 * MDP(16) - t163 * MDP(17)) * pkin(2) + 0.2e1 * (t153 * MDP(23) + t114 * MDP(25) - t107 * MDP(29)) * t136; t125 * MDP(16) - t126 * MDP(17) - t102 * MDP(25) + (-t132 * t151 - t199) * MDP(26) + t105 * MDP(27) + (t101 * t148 - t102 * t151) * MDP(28) + (t129 + t105) * MDP(30) + (t132 * t146 + t199) * MDP(31) + (t100 * t148 - t146 * t99) * MDP(32) + t178 * t164 + (-t190 + (-pkin(5) - t146) * MDP(29) - t218 - t188 + t185 * (-qJ(5) - t148)) * t167 + t220; (-t137 * t151 - t198) * MDP(26) + (-t123 * t151 + t124 * t148) * MDP(28) + (t137 * t146 + t198) * MDP(31) + (-t110 * t146 + t111 * t148) * MDP(32) + t223; (t151 ^ 2 + t145) * MDP(28) + (t146 ^ 2 + t145) * MDP(32) + 0.2e1 * t218 + 0.2e1 * t190 + t146 * t209 + t148 * t219 + t188; (-0.2e1 * t157 - t192) * MDP(25) + (-pkin(4) * t132 - t201) * MDP(26) + t182 * MDP(27) + (-pkin(4) * t102 + qJ(5) * t101) * MDP(28) + (t129 + t182) * MDP(30) + (t132 * t168 + t201) * MDP(31) + (qJ(5) * t100 - t168 * t99) * MDP(32) + (-MDP(22) + (-pkin(5) - t168) * MDP(29)) * t167 + t220; (-pkin(4) * t137 - t200) * MDP(26) + (-pkin(4) * t123 + qJ(5) * t124) * MDP(28) + (t137 * t168 + t200) * MDP(31) + (qJ(5) * t111 - t110 * t168) * MDP(32) + t173; MDP(22) + t171 * MDP(25) + (pkin(4) * t151 + t143) * MDP(28) + (0.2e1 * pkin(5) + t171) * MDP(29) + (t146 * t168 + t143) * MDP(32) + t185 * (0.2e1 * qJ(5) + t155) + (-t189 + (MDP(23) + t187) * t165) * pkin(3); MDP(22) + pkin(4) * t211 + (pkin(4) ^ 2 + t170) * MDP(28) + t168 * t209 + (t168 ^ 2 + t170) * MDP(32) + qJ(5) * t219; t102 * MDP(28) + t99 * MDP(32) + t186 * t132 + t187 * t167; t123 * MDP(28) + t110 * MDP(32) + t186 * t137; -MDP(28) * t151 - MDP(32) * t146 - t187; -pkin(4) * MDP(28) - MDP(32) * t168 - t187; MDP(28) + MDP(32); -t131 * MDP(29) + t132 * MDP(30) + t103 * MDP(32); -t136 * MDP(29) + t137 * MDP(30) + MDP(32) * t107; 0; 0; 0; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
