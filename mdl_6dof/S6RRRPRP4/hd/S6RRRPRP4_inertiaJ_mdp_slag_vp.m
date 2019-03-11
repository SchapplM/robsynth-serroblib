% Calculate joint inertia matrix for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP4_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:45:59
% EndTime: 2019-03-09 16:46:01
% DurationCPUTime: 0.68s
% Computational Cost: add. (970->186), mult. (1610->237), div. (0->0), fcn. (1618->6), ass. (0->88)
t214 = MDP(16) - MDP(19);
t213 = -MDP(17) + MDP(20);
t205 = 2 * MDP(20);
t156 = cos(qJ(5));
t208 = t156 * MDP(28);
t153 = sin(qJ(5));
t209 = t153 * MDP(27);
t212 = 0.2e1 * t208 + 0.2e1 * t209 + t205;
t211 = -t153 * pkin(5) + qJ(6) * t156;
t150 = t153 ^ 2;
t151 = t156 ^ 2;
t138 = t150 + t151;
t173 = MDP(32) * t138;
t174 = t138 * MDP(30);
t157 = cos(qJ(3));
t144 = -pkin(2) * t157 - pkin(3);
t210 = t144 * MDP(21);
t158 = cos(qJ(2));
t207 = 0.2e1 * t158;
t206 = -2 * MDP(19);
t204 = 2 * MDP(29);
t203 = -0.2e1 * MDP(30);
t202 = pkin(3) + pkin(9);
t201 = -pkin(8) - pkin(7);
t154 = sin(qJ(3));
t200 = pkin(2) * t154;
t155 = sin(qJ(2));
t129 = t154 * t158 + t155 * t157;
t199 = pkin(5) * t129;
t198 = pkin(3) * MDP(21);
t128 = t154 * t155 - t157 * t158;
t197 = qJ(4) * t128;
t196 = qJ(6) * t129;
t141 = qJ(4) + t200;
t194 = t128 * t141;
t140 = -pkin(9) + t144;
t193 = t129 * t140;
t192 = t129 * t202;
t191 = t153 * t156;
t190 = qJ(4) + t141;
t145 = -pkin(2) * t158 - pkin(1);
t163 = -qJ(4) * t129 + t145;
t104 = t202 * t128 + t163;
t135 = t201 * t155;
t136 = t201 * t158;
t117 = -t135 * t157 - t136 * t154;
t109 = pkin(4) * t129 + t117;
t100 = t156 * t104 + t153 * t109;
t132 = qJ(4) - t211;
t123 = t132 + t200;
t189 = t123 + t132;
t188 = MDP(32) * t140;
t118 = t135 * t154 - t136 * t157;
t133 = pkin(5) * t156 + qJ(6) * t153;
t102 = (-pkin(4) - t133) * t128 + t118;
t187 = t102 * MDP(31);
t186 = t123 * MDP(32);
t185 = t129 * MDP(25);
t184 = t132 * MDP(32);
t148 = t156 * MDP(24);
t183 = t156 * MDP(32);
t182 = t202 * MDP(32);
t181 = MDP(19) - t174;
t180 = -MDP(28) + MDP(31);
t178 = MDP(23) * t191;
t177 = t151 * MDP(22) + MDP(15) - 0.2e1 * t178;
t176 = -MDP(32) * pkin(5) - MDP(29);
t175 = t104 * t153 - t156 * t109;
t172 = -t153 * MDP(25) - t133 * MDP(30) + t148;
t171 = t192 + t197;
t170 = -t123 * t128 + t193;
t169 = -t128 * t132 - t192;
t168 = -t193 + t194;
t167 = -t175 * MDP(27) - t100 * MDP(28);
t166 = t153 * MDP(24) + t156 * MDP(25);
t164 = -0.2e1 * t156 * MDP(31) + t153 * t204;
t162 = (MDP(27) + MDP(29)) * t156 + t180 * t153;
t161 = (MDP(27) - t176) * t156 + (MDP(32) * qJ(6) + t180) * t153;
t110 = -pkin(4) * t128 + t118;
t97 = t100 + t196;
t98 = t175 - t199;
t96 = t153 * t97 - t156 * t98;
t160 = t102 * t153 * MDP(29) - t96 * MDP(30) + (t148 + MDP(13)) * t129 + t213 * t118 - t214 * t117 + (t208 + t209) * t110 + ((-t150 + t151) * MDP(23) + MDP(22) * t191 - MDP(14)) * t128;
t147 = t156 * MDP(30);
t130 = t138 * t202;
t120 = t138 * t140;
t112 = pkin(3) * t128 + t163;
t1 = [(t112 ^ 2 + t117 ^ 2 + t118 ^ 2) * MDP(21) + (t102 ^ 2 + t97 ^ 2 + t98 ^ 2) * MDP(32) + pkin(1) * MDP(9) * t207 + MDP(1) + (MDP(11) + MDP(26)) * t129 ^ 2 + (t150 * MDP(22) + 0.2e1 * t178) * t128 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t155 + MDP(5) * t207) * t155 + (0.2e1 * t145 * MDP(16) + t112 * t206) * t128 + 0.2e1 * (-t118 * MDP(18) + (t153 * t98 + t156 * t97) * MDP(30) + (-t156 * MDP(27) + t153 * MDP(28)) * t110 + (-t156 * MDP(29) - t153 * MDP(31)) * t102) * t128 + 0.2e1 * (t145 * MDP(17) - t112 * MDP(20) + (-MDP(12) + t166) * t128 + t117 * MDP(18) - t98 * MDP(29) + t97 * MDP(31) + t167) * t129; (-t158 * MDP(10) - t155 * MDP(9)) * pkin(7) + (t129 * t144 - t194) * MDP(18) + (t117 * t144 + t118 * t141) * MDP(21) + (t168 * MDP(28) + t170 * MDP(31) + t97 * t188 - t185) * t153 + (-t168 * MDP(27) + t170 * MDP(29) - t98 * t188 - t187) * t156 + t160 + t102 * t186 + t155 * MDP(6) + t158 * MDP(7); MDP(8) + t140 ^ 2 * t173 + (t164 + t186) * t123 + 0.2e1 * (t157 * MDP(16) - t154 * MDP(17)) * pkin(2) + t120 * t203 + t177 + ((2 * MDP(19)) + t210) * t144 + (t141 * MDP(21) + t212) * t141; (-pkin(3) * t117 + qJ(4) * t118) * MDP(21) + (-pkin(3) * t129 - t197) * MDP(18) + (t171 * MDP(28) + t169 * MDP(31) - t97 * t182 - t185) * t153 + (-t171 * MDP(27) + t169 * MDP(29) + t98 * t182 - t187) * t156 + t160 + t102 * t184; pkin(3) * t206 + qJ(4) * t205 + (-pkin(3) * t144 + qJ(4) * t141) * MDP(21) + t123 * t184 + t202 * t174 + (-t173 * t202 - t174) * t140 + (t213 * t154 + t214 * t157) * pkin(2) + (t190 * MDP(28) - t189 * MDP(31)) * t156 + (t190 * MDP(27) + t189 * MDP(29)) * t153 + t177; -t130 * t203 + t202 ^ 2 * t173 + (t164 + t184) * t132 + (t206 + t198) * pkin(3) + (qJ(4) * MDP(21) + t212) * qJ(4) + t177; t117 * MDP(21) + MDP(32) * t96 + (MDP(18) + t162) * t129; MDP(32) * t120 + t181 + t210; -MDP(32) * t130 + t181 - t198; MDP(21) + t173; t129 * MDP(26) + (-t175 + 0.2e1 * t199) * MDP(29) + (t100 + 0.2e1 * t196) * MDP(31) + (-pkin(5) * t98 + qJ(6) * t97) * MDP(32) + (t211 * MDP(30) + t166) * t128 + t167; t161 * t140 + t172; -t161 * t202 + t172; MDP(32) * t133 + t162; MDP(26) + pkin(5) * t204 + 0.2e1 * qJ(6) * MDP(31) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(32); t128 * t153 * MDP(30) - t129 * MDP(29) + t98 * MDP(32); -t140 * t183 + t147; t156 * t182 + t147; -t183; t176; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
