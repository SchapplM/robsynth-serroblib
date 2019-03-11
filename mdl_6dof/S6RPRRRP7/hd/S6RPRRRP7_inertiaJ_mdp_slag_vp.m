% Calculate joint inertia matrix for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRP7_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:21:02
% EndTime: 2019-03-09 06:21:05
% DurationCPUTime: 0.94s
% Computational Cost: add. (1413->195), mult. (2605->274), div. (0->0), fcn. (2923->8), ass. (0->84)
t156 = sin(qJ(4));
t159 = cos(qJ(4));
t194 = -pkin(9) - pkin(8);
t140 = t194 * t159;
t155 = sin(qJ(5));
t158 = cos(qJ(5));
t173 = t194 * t156;
t116 = -t140 * t155 - t158 * t173;
t117 = -t158 * t140 + t155 * t173;
t136 = t155 * t159 + t156 * t158;
t174 = -MDP(28) + MDP(31);
t175 = MDP(27) + MDP(29);
t200 = -t155 * t156 + t158 * t159;
t163 = t136 * MDP(24) + MDP(25) * t200 - t116 * t175 + t117 * t174;
t164 = -MDP(20) * t156 - MDP(21) * t159;
t205 = t156 * MDP(17) + t159 * MDP(18) + pkin(8) * t164 + t163;
t201 = -MDP(27) * t200 + t136 * MDP(28);
t154 = cos(pkin(10));
t143 = -pkin(2) * t154 - pkin(1);
t199 = 0.2e1 * t143;
t198 = -2 * MDP(23);
t197 = 2 * MDP(29);
t196 = 2 * MDP(30);
t195 = 2 * MDP(31);
t193 = cos(qJ(3));
t192 = pkin(1) * MDP(7);
t191 = pkin(4) * t158;
t153 = sin(pkin(10));
t157 = sin(qJ(3));
t132 = t153 * t157 - t154 * t193;
t122 = t132 * pkin(5);
t148 = t155 * pkin(4);
t190 = pkin(7) + qJ(2);
t133 = t153 * t193 + t157 * t154;
t186 = t133 * t159;
t108 = pkin(3) * t132 - pkin(8) * t133 + t143;
t138 = t190 * t153;
t139 = t190 * t154;
t111 = -t157 * t138 + t139 * t193;
t99 = t159 * t108 - t111 * t156;
t96 = pkin(4) * t132 - pkin(9) * t186 + t99;
t188 = t111 * t159;
t98 = t188 + (-pkin(9) * t133 + t108) * t156;
t92 = t155 * t96 + t158 * t98;
t106 = t200 * t133;
t189 = t106 * t200;
t187 = t133 * t156;
t185 = t153 * MDP(5);
t184 = t154 * MDP(4);
t182 = t156 * t159;
t179 = MDP(22) * t106;
t105 = t136 * t133;
t103 = t105 * MDP(25);
t104 = t106 * MDP(24);
t178 = t133 * MDP(14);
t123 = t200 * MDP(29);
t126 = t136 * MDP(31);
t177 = t155 * MDP(28);
t176 = MDP(19) + MDP(26);
t121 = t132 * qJ(6);
t89 = t121 + t92;
t172 = t132 * MDP(26) - t103 + t104;
t147 = -pkin(4) * t159 - pkin(3);
t171 = MDP(16) * t182;
t170 = t155 * t98 - t158 * t96;
t169 = pkin(5) * t197 + MDP(26);
t168 = t126 + t123 - t201;
t90 = -t122 + t170;
t110 = t138 * t193 + t157 * t139;
t167 = pkin(5) * t200 + qJ(6) * t136;
t166 = MDP(17) * t159 - MDP(18) * t156;
t165 = t159 * MDP(20) - t156 * MDP(21);
t102 = pkin(4) * t187 + t110;
t162 = MDP(13) + t165;
t152 = t159 ^ 2;
t151 = t156 ^ 2;
t145 = pkin(5) + t191;
t142 = t148 + qJ(6);
t134 = t136 ^ 2;
t109 = t147 - t167;
t101 = t136 * t105;
t100 = t108 * t156 + t188;
t93 = t105 * pkin(5) - t106 * qJ(6) + t102;
t1 = [(t89 ^ 2 + t90 ^ 2 + t93 ^ 2) * MDP(32) + t178 * t199 + MDP(1) + t176 * t132 ^ 2 + (t105 * t198 + t179) * t106 + (0.2e1 * t184 - 0.2e1 * t185 + t192) * pkin(1) + (t152 * MDP(15) + MDP(8) - 0.2e1 * t171) * t133 ^ 2 + (MDP(13) * t199 + 0.2e1 * t104 - 0.2e1 * t103 + 0.2e1 * (-MDP(9) + t166) * t133) * t132 + (-t105 * t89 + t106 * t90) * t196 + 0.2e1 * (t102 * t105 - t132 * t170) * MDP(27) + (-t106 * t93 + t132 * t89) * t195 + (t105 * t93 - t132 * t90) * t197 + 0.2e1 * (t102 * t106 - t132 * t92) * MDP(28) + 0.2e1 * (t110 * t187 + t132 * t99) * MDP(20) + 0.2e1 * (-t100 * t132 + t110 * t186) * MDP(21) + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t153 ^ 2 + t154 ^ 2) * qJ(2); -t184 + t185 - t192 + t178 + (-t101 - t189) * MDP(30) + (t136 * t89 - t200 * t90) * MDP(32) + (t136 * t174 + t175 * t200 + t162) * t132; MDP(7) + (t200 ^ 2 + t134) * MDP(32); -t111 * MDP(14) + t136 * t179 + (-t101 + t189) * MDP(23) + (-t102 * t200 + t105 * t147) * MDP(27) + (t102 * t136 + t106 * t147) * MDP(28) + (t105 * t109 - t200 * t93) * MDP(29) + (-t105 * t117 + t106 * t116 + t136 * t90 + t200 * t89) * MDP(30) + (-t106 * t109 - t136 * t93) * MDP(31) + (t109 * t93 + t116 * t90 + t117 * t89) * MDP(32) - t162 * t110 + (MDP(10) + MDP(15) * t182 + (-t151 + t152) * MDP(16) + t164 * pkin(3)) * t133 + (-MDP(11) + t205) * t132; (-t116 * t200 + t117 * t136) * MDP(32); MDP(12) + t151 * MDP(15) + 0.2e1 * t171 + t134 * MDP(22) - t136 * t200 * t198 + (t116 * t136 + t117 * t200) * t196 + (t116 ^ 2 + t117 ^ 2) * MDP(32) + (MDP(32) * t109 - 0.2e1 * t123 - 0.2e1 * t126) * t109 + 0.2e1 * t201 * t147 + 0.2e1 * t165 * pkin(3); t132 * MDP(19) + t99 * MDP(20) - t100 * MDP(21) + (t132 * t191 - t170) * MDP(27) + (-t132 * t148 - t92) * MDP(28) + (t132 * t145 - t90) * MDP(29) + (-t105 * t142 - t106 * t145) * MDP(30) + (t132 * t142 + t89) * MDP(31) + (t142 * t89 - t145 * t90) * MDP(32) + t166 * t133 + t172; (t136 * t142 + t145 * t200) * MDP(32) + t165 + t168; (-t136 * t145 + t142 * t200) * MDP(30) + (-t116 * t145 + t117 * t142) * MDP(32) + t205; (t142 ^ 2 + t145 ^ 2) * MDP(32) + 0.2e1 * (MDP(27) * t158 - t177) * pkin(4) + t145 * t197 + t142 * t195 + t176; -t170 * MDP(27) - t92 * MDP(28) + (-t90 + t122) * MDP(29) + (-pkin(5) * t106 - qJ(6) * t105) * MDP(30) + (0.2e1 * t121 + t92) * MDP(31) + (-pkin(5) * t90 + qJ(6) * t89) * MDP(32) + t172; MDP(32) * t167 + t168; (-pkin(5) * t136 + qJ(6) * t200) * MDP(30) + (-pkin(5) * t116 + qJ(6) * t117) * MDP(32) + t163; (0.2e1 * qJ(6) + t148) * MDP(31) + (pkin(5) * t145 + qJ(6) * t142) * MDP(32) + (t158 * t175 - t177) * pkin(4) + t169; qJ(6) * t195 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(32) + t169; -t132 * MDP(29) + t106 * MDP(30) + MDP(32) * t90; -t200 * MDP(32); MDP(30) * t136 + MDP(32) * t116; -MDP(32) * t145 - MDP(29); -MDP(32) * pkin(5) - MDP(29); MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
