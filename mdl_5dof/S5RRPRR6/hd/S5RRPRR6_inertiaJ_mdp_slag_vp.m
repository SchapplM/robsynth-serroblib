% Calculate joint inertia matrix for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRR6_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:17:26
% EndTime: 2022-01-20 11:17:28
% DurationCPUTime: 0.39s
% Computational Cost: add. (445->113), mult. (877->151), div. (0->0), fcn. (780->8), ass. (0->82)
t159 = sin(qJ(2));
t147 = t159 * pkin(1) + qJ(3);
t155 = sin(pkin(9));
t153 = t155 ^ 2;
t156 = cos(pkin(9));
t154 = t156 ^ 2;
t178 = t153 + t154;
t181 = t178 * t147;
t158 = sin(qJ(4));
t161 = cos(qJ(4));
t157 = sin(qJ(5));
t160 = cos(qJ(5));
t129 = -t157 * t158 + t160 * t161;
t130 = t157 * t161 + t160 * t158;
t183 = t129 * MDP(22) - t130 * MDP(23);
t206 = t161 * MDP(15) - t158 * MDP(16) + t183;
t174 = MDP(14) + MDP(21);
t204 = 2 * MDP(8);
t201 = 0.2e1 * MDP(15);
t200 = 0.2e1 * MDP(16);
t199 = 0.2e1 * MDP(22);
t198 = 0.2e1 * MDP(23);
t197 = pkin(8) * t155;
t162 = cos(qJ(2));
t196 = t162 * pkin(1);
t120 = t129 * t155;
t191 = t155 * t158;
t145 = pkin(4) * t191;
t122 = t155 * t147 + t145;
t136 = -t156 * pkin(3) - t155 * pkin(7) - pkin(2);
t124 = t136 - t196;
t187 = t156 * t161;
t171 = t147 * t187;
t100 = t171 + (t124 - t197) * t158;
t121 = t161 * t124;
t190 = t155 * t161;
t173 = pkin(8) * t190;
t98 = -t173 + t121 + (-t147 * t158 - pkin(4)) * t156;
t93 = t160 * t100 + t157 * t98;
t195 = t122 * t120 + t93 * t156;
t128 = t155 * qJ(3) + t145;
t127 = t161 * t136;
t105 = -t173 + t127 + (-qJ(3) * t158 - pkin(4)) * t156;
t172 = qJ(3) * t187;
t110 = t172 + (t136 - t197) * t158;
t96 = t157 * t105 + t160 * t110;
t194 = t128 * t120 + t96 * t156;
t193 = t153 * t158;
t192 = t153 * t161;
t189 = t156 * MDP(7);
t188 = t156 * t158;
t92 = -t157 * t100 + t160 * t98;
t186 = t92 * MDP(22);
t95 = t160 * t105 - t157 * t110;
t185 = t95 * MDP(22);
t119 = t130 * t155;
t116 = t119 * MDP(20);
t118 = t120 * MDP(19);
t184 = t118 - t116;
t109 = t158 * t124 + t171;
t182 = t109 * t156 + t147 * t192;
t115 = t158 * t136 + t172;
t180 = qJ(3) * t192 + t115 * t156;
t179 = t178 * qJ(3);
t177 = MDP(22) * t160;
t170 = MDP(13) * t191;
t141 = MDP(12) * t190;
t169 = -t156 * MDP(21) + t184;
t168 = (t162 * MDP(5) - t159 * MDP(6)) * pkin(1);
t167 = (-MDP(23) * t157 + t177) * pkin(4);
t166 = (-MDP(7) - t206) * t156;
t165 = t161 ^ 2 * t153 * MDP(10) - 0.2e1 * t158 * MDP(11) * t192 + MDP(4) - 0.2e1 * (t118 + t141) * t156 + 0.2e1 * (t116 + t170) * t156 + t174 * t154 + (MDP(17) * t120 - 0.2e1 * MDP(18) * t119) * t120;
t164 = t141 + (-pkin(4) * t177 - t174) * t156 - t170 + t184;
t148 = -pkin(2) - t196;
t144 = t157 * pkin(4) * t156;
t142 = qJ(3) * t193;
t131 = t147 * t193;
t114 = -qJ(3) * t188 + t127;
t108 = -t147 * t188 + t121;
t106 = t128 * t119;
t102 = t122 * t119;
t1 = [t165 + t181 * t204 + (-t108 * t156 + t131) * t201 + t182 * t200 + (-t92 * t156 + t102) * t199 + t195 * t198 + (t147 ^ 2 * t178 + t148 ^ 2) * MDP(9) + 0.2e1 * t168 - 0.2e1 * t148 * t189 + MDP(1); t165 + (t180 + t182) * MDP(16) + (t179 + t181) * MDP(8) + t168 + (t131 + t142) * MDP(15) + (t102 + t106) * MDP(22) + ((pkin(2) - t148) * MDP(7) + (-t108 - t114) * MDP(15) + (-t92 - t95) * MDP(22)) * t156 + (-t148 * pkin(2) + qJ(3) * t181) * MDP(9) + (t194 + t195) * MDP(23); t165 + t179 * t204 + (-t114 * t156 + t142) * t201 + t180 * t200 + (-t95 * t156 + t106) * t199 + t194 * t198 + (qJ(3) ^ 2 * t178 + pkin(2) ^ 2) * MDP(9) + 0.2e1 * pkin(2) * t189; t148 * MDP(9) + t166; -pkin(2) * MDP(9) + t166; MDP(9); t108 * MDP(15) - t109 * MDP(16) + t186 + (t144 - t93) * MDP(23) + t164; t114 * MDP(15) - t115 * MDP(16) + t185 + (t144 - t96) * MDP(23) + t164; t206; 0.2e1 * t167 + t174; -t93 * MDP(23) + t169 + t186; -t96 * MDP(23) + t169 + t185; t183; MDP(21) + t167; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
