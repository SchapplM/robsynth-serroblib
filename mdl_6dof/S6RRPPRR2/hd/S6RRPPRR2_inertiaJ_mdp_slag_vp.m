% Calculate joint inertia matrix for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPPRR2_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:52:38
% EndTime: 2019-03-09 08:52:40
% DurationCPUTime: 0.66s
% Computational Cost: add. (1284->153), mult. (2448->223), div. (0->0), fcn. (2856->10), ass. (0->86)
t158 = sin(pkin(10));
t160 = cos(pkin(10));
t163 = sin(qJ(2));
t197 = cos(qJ(2));
t143 = t158 * t197 + t160 * t163;
t157 = sin(pkin(11));
t159 = cos(pkin(11));
t162 = sin(qJ(5));
t165 = cos(qJ(5));
t144 = t157 * t165 + t159 * t162;
t116 = t144 * t143;
t142 = t157 * t162 - t165 * t159;
t117 = t142 * t143;
t141 = t158 * t163 - t160 * t197;
t154 = -t197 * pkin(2) - pkin(1);
t126 = t141 * pkin(3) - t143 * qJ(4) + t154;
t207 = -qJ(3) - pkin(7);
t147 = t207 * t163;
t148 = t207 * t197;
t132 = t147 * t158 - t148 * t160;
t110 = t159 * t126 - t132 * t157;
t192 = t143 * t159;
t103 = pkin(4) * t141 - pkin(8) * t192 + t110;
t111 = t157 * t126 + t159 * t132;
t193 = t143 * t157;
t106 = -pkin(8) * t193 + t111;
t96 = t165 * t103 - t106 * t162;
t97 = t103 * t162 + t106 * t165;
t209 = -t117 * MDP(19) - t116 * MDP(20) + t96 * MDP(22) - t97 * MDP(23);
t161 = sin(qJ(6));
t164 = cos(qJ(6));
t107 = t164 * t116 - t117 * t161;
t108 = -t116 * t161 - t117 * t164;
t206 = t108 * MDP(26) - t107 * MDP(27);
t150 = pkin(2) * t158 + qJ(4);
t195 = pkin(8) + t150;
t137 = t195 * t157;
t138 = t195 * t159;
t120 = -t165 * t137 - t138 * t162;
t121 = -t137 * t162 + t138 * t165;
t112 = -pkin(9) * t144 + t120;
t113 = -pkin(9) * t142 + t121;
t127 = t164 * t142 + t144 * t161;
t128 = -t142 * t161 + t144 * t164;
t179 = t128 * MDP(26) - t127 * MDP(27) + (t112 * t164 - t113 * t161) * MDP(29) - (t112 * t161 + t113 * t164) * MDP(30);
t205 = t144 * MDP(19) - t142 * MDP(20) + t120 * MDP(22) - t121 * MDP(23) + t179;
t186 = t142 * MDP(22);
t122 = t127 * MDP(29);
t191 = -t128 * MDP(30) - t122;
t204 = -t144 * MDP(23) - t186 + t191;
t190 = t157 ^ 2 + t159 ^ 2;
t203 = t190 * MDP(15);
t130 = -t160 * t147 - t148 * t158;
t202 = t130 ^ 2;
t153 = -pkin(2) * t160 - pkin(3);
t146 = -pkin(4) * t159 + t153;
t133 = pkin(5) * t142 + t146;
t201 = 0.2e1 * t133;
t200 = 0.2e1 * t146;
t199 = -2 * MDP(18);
t198 = -2 * MDP(25);
t196 = pkin(5) * t141;
t95 = -pkin(9) * t116 + t97;
t194 = t164 * t95;
t189 = MDP(16) * t153;
t188 = t108 * MDP(24);
t187 = t117 * MDP(17);
t185 = t157 * MDP(14);
t184 = t159 * MDP(13);
t183 = 0.2e1 * t197;
t182 = MDP(21) + MDP(28);
t181 = t141 * MDP(28) + t206;
t94 = pkin(9) * t117 + t196 + t96;
t91 = -t161 * t95 + t164 * t94;
t180 = t190 * MDP(16);
t114 = pkin(4) * t193 + t130;
t177 = t91 * MDP(29) - (t161 * t94 + t194) * MDP(30);
t176 = t110 * t159 + t111 * t157;
t175 = -t110 * t157 + t111 * t159;
t174 = MDP(13) * t157 + MDP(14) * t159;
t172 = t116 * MDP(22) - t117 * MDP(23);
t170 = t107 * MDP(29) + t108 * MDP(30);
t169 = (MDP(29) * t164 - MDP(30) * t161) * pkin(5);
t167 = t184 - t185 + t204;
t109 = pkin(5) * t116 + t114;
t1 = [pkin(1) * MDP(9) * t183 + (t132 ^ 2 + t154 ^ 2 + t202) * MDP(12) + (t110 ^ 2 + t111 ^ 2 + t202) * MDP(16) + MDP(1) + t182 * t141 ^ 2 - (t116 * t199 - t187) * t117 + (t107 * t198 + t188) * t108 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t163 + MDP(5) * t183) * t163 + 0.2e1 * t172 * t114 + 0.2e1 * t170 * t109 + 0.2e1 * (-t176 * MDP(15) + (MDP(11) + t174) * t130) * t143 + 0.2e1 * (-t132 * MDP(11) + t110 * MDP(13) - t111 * MDP(14) + t177 + t206 + t209) * t141; t163 * MDP(6) + t197 * MDP(7) + (-t130 * t159 + t153 * t193) * MDP(13) + (t130 * t157 + t153 * t192) * MDP(14) + t175 * MDP(15) + (t130 * t153 + t150 * t175) * MDP(16) - t144 * t187 + (-t116 * t144 + t117 * t142) * MDP(18) + (t114 * t142 + t116 * t146) * MDP(22) + (t114 * t144 - t117 * t146) * MDP(23) + t128 * t188 + (-t107 * t128 - t108 * t127) * MDP(25) + (t107 * t133 + t109 * t127) * MDP(29) + (t108 * t133 + t109 * t128) * MDP(30) + (-t197 * MDP(10) - t163 * MDP(9)) * pkin(7) + (-t150 * t174 + t205) * t141 + ((-t141 * t158 - t143 * t160) * MDP(11) + (-t130 * t160 + t132 * t158) * MDP(12)) * pkin(2); t186 * t200 + t122 * t201 + MDP(8) + (t158 ^ 2 + t160 ^ 2) * MDP(12) * pkin(2) ^ 2 + (-0.2e1 * t184 + 0.2e1 * t185 + t189) * t153 + (MDP(17) * t144 + MDP(23) * t200 + t142 * t199) * t144 + (MDP(24) * t128 + MDP(30) * t201 + t127 * t198) * t128 + (t180 * t150 + 0.2e1 * t203) * t150; t154 * MDP(12) + MDP(16) * t176 + t141 * t167 - t143 * t203; 0; MDP(12) + t180; t130 * MDP(16) + t143 * t174 + t170 + t172; -t167 + t189; 0; MDP(16); t141 * MDP(21) + (t164 * t196 + t91) * MDP(29) + (-t194 + (-t94 - t196) * t161) * MDP(30) + t181 + t209; t205; t204; 0; 0.2e1 * t169 + t182; t177 + t181; t179; t191; 0; MDP(28) + t169; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
