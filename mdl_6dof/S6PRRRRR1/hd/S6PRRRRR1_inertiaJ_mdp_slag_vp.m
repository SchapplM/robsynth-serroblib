% Calculate joint inertia matrix for
% S6PRRRRR1
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
%   see S6PRRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR1_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:39:53
% EndTime: 2019-03-09 00:39:54
% DurationCPUTime: 0.52s
% Computational Cost: add. (820->146), mult. (1643->202), div. (0->0), fcn. (1928->12), ass. (0->89)
t150 = sin(qJ(6));
t155 = cos(qJ(6));
t187 = t150 * MDP(28) + t155 * MDP(29);
t183 = MDP(31) * t155;
t200 = -MDP(32) * t150 + t183;
t163 = 0.2e1 * t200;
t152 = sin(qJ(4));
t153 = sin(qJ(3));
t157 = cos(qJ(4));
t158 = cos(qJ(3));
t130 = t152 * t153 - t157 * t158;
t141 = -t158 * pkin(3) - pkin(2);
t120 = t130 * pkin(4) + t141;
t199 = 0.2e1 * t120;
t198 = 0.2e1 * t141;
t197 = pkin(8) + pkin(9);
t196 = pkin(3) * t152;
t195 = pkin(5) * t150;
t156 = cos(qJ(5));
t194 = t156 * pkin(4);
t131 = t152 * t158 + t157 * t153;
t133 = t197 * t153;
t134 = t197 * t158;
t174 = -t157 * t133 - t152 * t134;
t106 = -t131 * pkin(10) + t174;
t169 = t152 * t133 - t157 * t134;
t107 = -t130 * pkin(10) - t169;
t151 = sin(qJ(5));
t98 = -t156 * t106 + t151 * t107;
t193 = t98 * t155;
t116 = -t151 * t130 + t156 * t131;
t192 = t116 * t150;
t191 = t116 * t155;
t148 = sin(pkin(6));
t154 = sin(qJ(2));
t190 = t148 * t154;
t159 = cos(qJ(2));
t189 = t148 * t159;
t188 = t150 * t155;
t186 = MDP(10) * t158;
t185 = MDP(17) * t130;
t184 = MDP(25) * t116;
t115 = t156 * t130 + t151 * t131;
t181 = t115 * MDP(30);
t140 = t157 * pkin(3) + pkin(4);
t135 = t156 * t140;
t123 = -t151 * t196 + t135;
t180 = t123 * MDP(24);
t124 = -t151 * t140 - t156 * t196;
t179 = t124 * MDP(25);
t178 = t151 * MDP(25);
t177 = t157 * MDP(17);
t176 = MDP(27) * t188;
t146 = t150 ^ 2;
t175 = t146 * MDP(26) + MDP(23) + 0.2e1 * t176;
t173 = MDP(16) + t175;
t172 = -pkin(5) * t116 - pkin(11) * t115;
t121 = -pkin(5) - t123;
t122 = pkin(11) - t124;
t171 = -t115 * t122 + t116 * t121;
t138 = t151 * pkin(4) + pkin(11);
t139 = -pkin(5) - t194;
t170 = -t115 * t138 + t116 * t139;
t168 = MDP(28) * t155 - MDP(29) * t150;
t166 = -MDP(31) * t150 - MDP(32) * t155;
t149 = cos(pkin(6));
t125 = t149 * t158 - t153 * t190;
t126 = t149 * t153 + t158 * t190;
t108 = t157 * t125 - t152 * t126;
t109 = t152 * t125 + t157 * t126;
t100 = -t156 * t108 + t151 * t109;
t101 = t151 * t108 + t156 * t109;
t165 = -t101 * MDP(25) + (-MDP(24) - t200) * t100;
t164 = (t156 * MDP(24) - t178) * pkin(4);
t147 = t155 ^ 2;
t99 = t151 * t106 + t156 * t107;
t162 = -t98 * MDP(24) - t99 * MDP(25) + ((-t146 + t147) * MDP(27) + MDP(26) * t188 + MDP(21)) * t116 + (-MDP(22) + t187) * t115;
t161 = t108 * MDP(17) - t109 * MDP(18) + t165;
t160 = t131 * MDP(14) - t130 * MDP(15) + t174 * MDP(17) + t169 * MDP(18) + t162;
t145 = pkin(5) * t155;
t136 = t139 * t150;
t119 = t121 * t150;
t95 = t115 * pkin(5) - t116 * pkin(11) + t120;
t92 = t98 * t150;
t90 = t155 * t101 - t150 * t189;
t89 = -t150 * t101 - t155 * t189;
t88 = t150 * t95 + t155 * t99;
t87 = -t150 * t99 + t155 * t95;
t1 = [MDP(1); (t100 * t192 + t89 * t115) * MDP(31) + (t100 * t191 - t90 * t115) * MDP(32) + (-t154 * MDP(4) + (-MDP(11) * t153 - MDP(18) * t131 - MDP(24) * t115 + MDP(3) - t184 - t185 + t186) * t159) * t148; 0.2e1 * pkin(2) * t186 + t185 * t198 + t184 * t199 + MDP(2) + (MDP(24) * t199 + t181 + 0.2e1 * (-MDP(20) + t168) * t116) * t115 + 0.2e1 * (t87 * t115 + t192 * t98) * MDP(31) + 0.2e1 * (-t88 * t115 + t191 * t98) * MDP(32) + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t153 + 0.2e1 * t158 * MDP(6)) * t153 + (MDP(12) * t131 - 0.2e1 * t130 * MDP(13) + MDP(18) * t198) * t131 + (t147 * MDP(26) + MDP(19) - 0.2e1 * t176) * t116 ^ 2; t125 * MDP(10) - t126 * MDP(11) + t161; (t150 * t171 - t193) * MDP(31) + (t155 * t171 + t92) * MDP(32) + t160 + (-t153 * MDP(10) - t158 * MDP(11)) * pkin(8) + t153 * MDP(7) + t158 * MDP(8); MDP(9) - t121 * t163 + 0.2e1 * (-t152 * MDP(18) + t177) * pkin(3) + 0.2e1 * t180 + 0.2e1 * t179 + t173; t161; (t150 * t170 - t193) * MDP(31) + (t155 * t170 + t92) * MDP(32) + t160; (t135 + t194) * MDP(24) + (t136 + t119) * MDP(32) + (-t121 - t139) * t183 + (-pkin(4) - t140) * t178 + (t177 + (-MDP(24) * t151 - MDP(25) * t156 - MDP(18)) * t152) * pkin(3) + t173; -t139 * t163 + 0.2e1 * t164 + t173; t165; (t150 * t172 - t193) * MDP(31) + (t155 * t172 + t92) * MDP(32) + t162; t180 + t179 + (-t121 * t155 + t145) * MDP(31) + (t119 - t195) * MDP(32) + t175; (-t139 * t155 + t145) * MDP(31) + (t136 - t195) * MDP(32) + t164 + t175; pkin(5) * t163 + t175; t89 * MDP(31) - t90 * MDP(32); t87 * MDP(31) - t88 * MDP(32) + t116 * t168 + t181; t122 * t166 + t187; t138 * t166 + t187; pkin(11) * t166 + t187; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
