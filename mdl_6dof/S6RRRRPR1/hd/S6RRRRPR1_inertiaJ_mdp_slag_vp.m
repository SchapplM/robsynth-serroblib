% Calculate joint inertia matrix for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR1_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:55:06
% EndTime: 2019-03-09 21:55:08
% DurationCPUTime: 0.62s
% Computational Cost: add. (1671->165), mult. (3071->235), div. (0->0), fcn. (3632->10), ass. (0->86)
t157 = sin(qJ(6));
t161 = cos(qJ(6));
t194 = t157 * MDP(29) + t161 * MDP(30);
t159 = sin(qJ(3));
t160 = sin(qJ(2));
t163 = cos(qJ(3));
t164 = cos(qJ(2));
t133 = t159 * t160 - t163 * t164;
t187 = -t164 * pkin(2) - pkin(1);
t123 = t133 * pkin(3) + t187;
t202 = 0.2e1 * t123;
t201 = 0.2e1 * t164;
t200 = pkin(7) + pkin(8);
t199 = pkin(2) * t159;
t158 = sin(qJ(4));
t198 = t158 * pkin(3);
t197 = MDP(26) * pkin(4);
t155 = sin(pkin(11));
t156 = cos(pkin(11));
t134 = t159 * t164 + t163 * t160;
t162 = cos(qJ(4));
t117 = -t158 * t133 + t162 * t134;
t138 = t200 * t160;
t139 = t200 * t164;
t178 = t159 * t138 - t163 * t139;
t172 = t133 * pkin(9) + t178;
t184 = -t163 * t138 - t159 * t139;
t175 = -t134 * pkin(9) + t184;
t168 = t158 * t172 + t162 * t175;
t166 = -t117 * qJ(5) + t168;
t167 = -t158 * t175 + t162 * t172;
t179 = t162 * t133 + t158 * t134;
t99 = -t179 * qJ(5) - t167;
t93 = t155 * t99 - t156 * t166;
t92 = t93 * t157;
t196 = t93 * t161;
t195 = t157 * t161;
t147 = t163 * pkin(2) + pkin(3);
t140 = t162 * t147;
t129 = -t158 * t199 + t140;
t126 = pkin(4) + t129;
t130 = t158 * t147 + t162 * t199;
t113 = t155 * t126 + t156 * t130;
t152 = t162 * pkin(3);
t146 = t152 + pkin(4);
t128 = t155 * t146 + t156 * t198;
t106 = t155 * t117 + t156 * t179;
t193 = t106 * MDP(31);
t192 = t129 * MDP(23);
t191 = t130 * MDP(24);
t190 = t158 * MDP(24);
t189 = t161 * MDP(32);
t188 = t163 * MDP(16);
t186 = MDP(28) * t195;
t153 = t157 ^ 2;
t185 = t153 * MDP(27) + MDP(22) + 0.2e1 * t186;
t183 = MDP(15) + t185;
t107 = t156 * t117 - t155 * t179;
t112 = t156 * t126 - t155 * t130;
t110 = -pkin(5) - t112;
t111 = pkin(10) + t113;
t182 = -t106 * t111 + t107 * t110;
t127 = t156 * t146 - t155 * t198;
t124 = -pkin(5) - t127;
t125 = pkin(10) + t128;
t181 = -t106 * t125 + t107 * t124;
t143 = t155 * pkin(4) + pkin(10);
t144 = -t156 * pkin(4) - pkin(5);
t180 = -t106 * t143 + t107 * t144;
t177 = -t157 * MDP(33) + t189;
t176 = -MDP(32) * t157 - MDP(33) * t161;
t174 = (t162 * MDP(23) - t190) * pkin(3);
t173 = (MDP(29) * t161 - MDP(30) * t157) * t107;
t171 = 0.2e1 * t177;
t154 = t161 ^ 2;
t170 = t117 * MDP(20) - t179 * MDP(21) + t168 * MDP(23) + t167 * MDP(24) + (MDP(27) * t195 + (-t153 + t154) * MDP(28)) * t107 + t194 * t106;
t169 = t134 * MDP(13) - t133 * MDP(14) + t184 * MDP(16) + t178 * MDP(17) + t170;
t108 = t179 * pkin(4) + t123;
t137 = t144 * t157;
t120 = t124 * t157;
t109 = t110 * t157;
t96 = t106 * pkin(5) - t107 * pkin(10) + t108;
t95 = t155 * t166 + t156 * t99;
t91 = t157 * t96 + t161 * t95;
t90 = -t157 * t95 + t161 * t96;
t1 = [MDP(1) + pkin(1) * MDP(9) * t201 + t179 * MDP(23) * t202 + (t108 ^ 2 + t93 ^ 2 + t95 ^ 2) * MDP(26) + (t154 * MDP(27) - 0.2e1 * t186) * t107 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t160 + MDP(5) * t201) * t160 + (MDP(18) * t117 - 0.2e1 * t179 * MDP(19) + MDP(24) * t202) * t117 + t193 * t106 + 0.2e1 * (-t95 * t106 + t93 * t107) * MDP(25) + 0.2e1 * (t90 * t106 + t107 * t92) * MDP(32) + 0.2e1 * (-t91 * t106 + t107 * t196) * MDP(33) + (MDP(11) * t134 - 0.2e1 * t133 * MDP(12)) * t134 + 0.2e1 * (t133 * MDP(16) + t134 * MDP(17)) * t187 + 0.2e1 * t173 * t106; t169 + t164 * MDP(7) + t160 * MDP(6) + (-t113 * t106 - t112 * t107) * MDP(25) + (-t93 * t112 + t95 * t113) * MDP(26) + (t182 * t157 - t196) * MDP(32) + (t182 * t161 + t92) * MDP(33) + (-t164 * MDP(10) - t160 * MDP(9)) * pkin(7); MDP(8) + (t112 ^ 2 + t113 ^ 2) * MDP(26) - t110 * t171 + 0.2e1 * (-t159 * MDP(17) + t188) * pkin(2) + 0.2e1 * t192 - 0.2e1 * t191 + t183; t169 + (-t128 * t106 - t127 * t107) * MDP(25) + (-t93 * t127 + t95 * t128) * MDP(26) + (t181 * t157 - t196) * MDP(32) + (t181 * t161 + t92) * MDP(33); (t140 + t152) * MDP(23) + (t112 * t127 + t113 * t128) * MDP(26) + (t120 + t109) * MDP(33) + (-t110 - t124) * t189 + (-pkin(3) - t147) * t190 + (t188 + (-MDP(23) * t158 - MDP(24) * t162 - MDP(17)) * t159) * pkin(2) + t183; (t127 ^ 2 + t128 ^ 2) * MDP(26) - t124 * t171 + 0.2e1 * t174 + t183; (t180 * t157 - t196) * MDP(32) + (t180 * t161 + t92) * MDP(33) + ((-t106 * t155 - t107 * t156) * MDP(25) + (t155 * t95 - t156 * t93) * MDP(26)) * pkin(4) + t170; t192 - t191 + (t137 + t109) * MDP(33) + (-t110 - t144) * t189 + (t112 * t156 + t113 * t155) * t197 + t185; (t137 + t120) * MDP(33) + (-t124 - t144) * t189 + (t127 * t156 + t128 * t155) * t197 + t174 + t185; (t155 ^ 2 + t156 ^ 2) * MDP(26) * pkin(4) ^ 2 - t144 * t171 + t185; t108 * MDP(26) + t177 * t106; 0; 0; 0; MDP(26); t90 * MDP(32) - t91 * MDP(33) + t173 + t193; t176 * t111 + t194; t176 * t125 + t194; t176 * t143 + t194; t177; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
