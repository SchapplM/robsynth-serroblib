% Calculate joint inertia matrix for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR1_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:04:08
% EndTime: 2019-03-09 18:04:10
% DurationCPUTime: 0.58s
% Computational Cost: add. (1541->155), mult. (2848->214), div. (0->0), fcn. (3404->10), ass. (0->86)
t157 = sin(qJ(6));
t161 = cos(qJ(6));
t190 = t157 * MDP(29) + t161 * MDP(30);
t183 = t161 * MDP(32);
t171 = -t157 * MDP(33) + t183;
t168 = 0.2e1 * t171;
t159 = sin(qJ(3));
t160 = sin(qJ(2));
t163 = cos(qJ(2));
t197 = cos(qJ(3));
t137 = t159 * t160 - t197 * t163;
t148 = -t163 * pkin(2) - pkin(1);
t127 = t137 * pkin(3) + t148;
t138 = t159 * t163 + t197 * t160;
t155 = sin(pkin(11));
t156 = cos(pkin(11));
t173 = t156 * t137 + t155 * t138;
t110 = t173 * pkin(4) + t127;
t201 = 0.2e1 * t110;
t200 = 0.2e1 * t148;
t199 = 0.2e1 * t163;
t198 = pkin(7) + pkin(8);
t196 = pkin(2) * t159;
t195 = pkin(5) * t157;
t194 = t155 * pkin(3);
t158 = sin(qJ(5));
t162 = cos(qJ(5));
t142 = t198 * t160;
t143 = t198 * t163;
t179 = -t197 * t142 - t159 * t143;
t119 = -t138 * qJ(4) + t179;
t169 = t159 * t142 - t197 * t143;
t120 = -t137 * qJ(4) - t169;
t100 = t156 * t119 - t155 * t120;
t122 = -t155 * t137 + t156 * t138;
t98 = -t122 * pkin(9) + t100;
t101 = t155 * t119 + t156 * t120;
t99 = -t173 * pkin(9) + t101;
t94 = t158 * t99 - t162 * t98;
t91 = t94 * t157;
t193 = t94 * t161;
t192 = t157 * t161;
t145 = t156 * pkin(3) + pkin(4);
t191 = t158 * t145;
t108 = t158 * t122 + t162 * t173;
t189 = t108 * MDP(31);
t109 = t162 * t122 - t158 * t173;
t188 = t109 * MDP(26);
t147 = t197 * pkin(2) + pkin(3);
t131 = t156 * t147 - t155 * t196;
t128 = pkin(4) + t131;
t133 = t155 * t147 + t156 * t196;
t180 = -t162 * t128 + t158 * t133;
t187 = t180 * MDP(25);
t174 = t158 * t128 + t162 * t133;
t186 = t174 * MDP(26);
t141 = t162 * t145;
t132 = -t158 * t194 + t141;
t185 = t132 * MDP(25);
t134 = -t162 * t194 - t191;
t184 = t134 * MDP(26);
t182 = MDP(28) * t192;
t153 = t157 ^ 2;
t181 = t153 * MDP(27) + MDP(24) + 0.2e1 * t182;
t178 = MDP(15) + t181;
t177 = -pkin(5) * t109 - pkin(10) * t108;
t112 = -pkin(5) + t180;
t113 = pkin(10) + t174;
t176 = -t108 * t113 + t109 * t112;
t129 = -pkin(5) - t132;
t130 = pkin(10) - t134;
t175 = -t108 * t130 + t109 * t129;
t172 = MDP(29) * t161 - MDP(30) * t157;
t170 = -MDP(32) * t157 - MDP(33) * t161;
t154 = t161 ^ 2;
t95 = t158 * t98 + t162 * t99;
t167 = -t94 * MDP(25) - t95 * MDP(26) + (MDP(27) * t192 + MDP(22) + (-t153 + t154) * MDP(28)) * t109 + (-MDP(23) + t190) * t108;
t166 = (t197 * MDP(16) - t159 * MDP(17)) * pkin(2);
t165 = t138 * MDP(13) - t137 * MDP(14) + t179 * MDP(16) + t169 * MDP(17) + t167;
t152 = pkin(5) * t161;
t126 = t129 * t157;
t111 = t112 * t157;
t96 = t108 * pkin(5) - t109 * pkin(10) + t110;
t90 = t157 * t96 + t161 * t95;
t89 = -t157 * t95 + t161 * t96;
t1 = [MDP(1) + pkin(1) * MDP(9) * t199 + t137 * MDP(16) * t200 + (t100 ^ 2 + t101 ^ 2 + t127 ^ 2) * MDP(19) + t188 * t201 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t160 + MDP(5) * t199) * t160 + (MDP(11) * t138 - 0.2e1 * t137 * MDP(12) + MDP(17) * t200) * t138 + (t154 * MDP(27) + MDP(20) - 0.2e1 * t182) * t109 ^ 2 + (MDP(25) * t201 + t189 + 0.2e1 * (-MDP(21) + t172) * t109) * t108 + 0.2e1 * (-t100 * t122 - t101 * t173) * MDP(18) + 0.2e1 * (t89 * t108 + t109 * t91) * MDP(32) + 0.2e1 * (-t90 * t108 + t109 * t193) * MDP(33); (t176 * t161 + t91) * MDP(33) + (t176 * t157 - t193) * MDP(32) + (-t163 * MDP(10) - t160 * MDP(9)) * pkin(7) + (t100 * t131 + t101 * t133) * MDP(19) + (-t131 * t122 - t133 * t173) * MDP(18) + t160 * MDP(6) + t163 * MDP(7) + t165; MDP(8) + (t131 ^ 2 + t133 ^ 2) * MDP(19) - t112 * t168 + 0.2e1 * t166 - 0.2e1 * t187 - 0.2e1 * t186 + t178; (t175 * t161 + t91) * MDP(33) + (t175 * t157 - t193) * MDP(32) + t165 + ((t100 * t156 + t101 * t155) * MDP(19) + (-t156 * t122 - t155 * t173) * MDP(18)) * pkin(3); (t141 - t180) * MDP(25) + (-t174 - t191) * MDP(26) + (t126 + t111) * MDP(33) + (-t112 - t129) * t183 + (t131 * t156 * MDP(19) + (t133 * MDP(19) - t158 * MDP(25) - t162 * MDP(26)) * t155) * pkin(3) + t166 + t178; (t155 ^ 2 + t156 ^ 2) * MDP(19) * pkin(3) ^ 2 - t129 * t168 + 0.2e1 * t185 + 0.2e1 * t184 + t178; t127 * MDP(19) + t188 + (MDP(25) + t171) * t108; 0; 0; MDP(19); (t177 * t157 - t193) * MDP(32) + (t177 * t161 + t91) * MDP(33) + t167; -t187 - t186 + (-t112 * t161 + t152) * MDP(32) + (t111 - t195) * MDP(33) + t181; t185 + t184 + (-t129 * t161 + t152) * MDP(32) + (t126 - t195) * MDP(33) + t181; 0; pkin(5) * t168 + t181; t89 * MDP(32) - t90 * MDP(33) + t172 * t109 + t189; t170 * t113 + t190; t170 * t130 + t190; t171; t170 * pkin(10) + t190; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
