% Calculate joint inertia matrix for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRPR5_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:28:31
% EndTime: 2019-03-08 23:28:33
% DurationCPUTime: 0.84s
% Computational Cost: add. (1346->216), mult. (3236->348), div. (0->0), fcn. (3749->14), ass. (0->96)
t163 = sin(pkin(7));
t214 = 0.2e1 * t163;
t169 = sin(qJ(4));
t173 = cos(qJ(4));
t213 = -(t169 * MDP(17) + t173 * MDP(18)) * pkin(10) + t169 * MDP(14) + t173 * MDP(15);
t212 = 2 * MDP(19);
t211 = 2 * MDP(26);
t210 = 2 * MDP(27);
t170 = sin(qJ(3));
t209 = pkin(2) * t170;
t174 = cos(qJ(3));
t208 = pkin(2) * t174;
t207 = -qJ(5) - pkin(10);
t162 = sin(pkin(13));
t165 = cos(pkin(13));
t149 = t162 * t173 + t165 * t169;
t168 = sin(qJ(6));
t206 = t149 * t168;
t172 = cos(qJ(6));
t205 = t149 * t172;
t204 = t163 * t170;
t203 = t163 * t174;
t166 = cos(pkin(7));
t202 = t166 * MDP(9);
t175 = cos(qJ(2));
t201 = t166 * t175;
t200 = t170 * MDP(7);
t188 = pkin(9) * t203;
t139 = t188 + (pkin(10) + t209) * t166;
t140 = (-pkin(3) * t174 - pkin(10) * t170 - pkin(2)) * t163;
t126 = -t169 * t139 + t173 * t140;
t144 = t169 * t166 + t173 * t204;
t117 = -pkin(4) * t203 - t144 * qJ(5) + t126;
t127 = t173 * t139 + t169 * t140;
t143 = -t173 * t166 + t169 * t204;
t121 = -t143 * qJ(5) + t127;
t111 = t162 * t117 + t165 * t121;
t199 = MDP(16) * t174;
t129 = -t162 * t143 + t165 * t144;
t124 = t168 * t129 + t172 * t203;
t198 = t124 * MDP(22);
t197 = t124 * MDP(24);
t125 = t172 * t129 - t168 * t203;
t196 = t125 * MDP(21);
t128 = t165 * t143 + t162 * t144;
t195 = t128 * MDP(23);
t194 = t128 * MDP(25);
t193 = t144 * MDP(12);
t148 = t162 * t169 - t165 * t173;
t192 = t148 * MDP(25);
t158 = -t173 * pkin(4) - pkin(3);
t191 = t158 * MDP(20);
t190 = t168 * MDP(21);
t189 = t173 * MDP(17);
t187 = t168 * t172 * MDP(22);
t186 = t207 * t169;
t110 = t165 * t117 - t162 * t121;
t164 = sin(pkin(6));
t167 = cos(pkin(6));
t171 = sin(qJ(2));
t133 = t167 * t204 + (t170 * t201 + t171 * t174) * t164;
t142 = -t164 * t175 * t163 + t167 * t166;
t185 = t133 * t169 - t142 * t173;
t184 = t144 * MDP(14) - t143 * MDP(15);
t181 = t172 * MDP(26) - t168 * MDP(27);
t180 = MDP(26) * t168 + MDP(27) * t172;
t154 = pkin(9) * t204;
t138 = t154 + (-pkin(3) - t208) * t166;
t179 = (MDP(23) * t172 - MDP(24) * t168) * t149;
t130 = t143 * pkin(4) + t138;
t178 = t168 * MDP(23) + t172 * MDP(24) - t180 * (t162 * pkin(4) + pkin(11));
t109 = -pkin(11) * t203 + t111;
t112 = t128 * pkin(5) - t129 * pkin(11) + t130;
t104 = -t168 * t109 + t172 * t112;
t105 = t172 * t109 + t168 * t112;
t177 = t125 * MDP(23) + t104 * MDP(26) - t105 * MDP(27) + t194 - t197;
t161 = t172 ^ 2;
t160 = t168 ^ 2;
t159 = t163 ^ 2;
t157 = -t165 * pkin(4) - pkin(5);
t152 = t207 * t173;
t146 = t166 * t209 + t188;
t145 = t166 * t208 - t154;
t136 = -t165 * t152 + t162 * t186;
t134 = -t162 * t152 - t165 * t186;
t132 = -t167 * t203 + (t170 * t171 - t174 * t201) * t164;
t131 = t148 * pkin(5) - t149 * pkin(11) + t158;
t123 = t133 * t173 + t142 * t169;
t119 = t168 * t131 + t172 * t136;
t118 = t172 * t131 - t168 * t136;
t115 = t165 * t123 - t162 * t185;
t113 = t162 * t123 + t165 * t185;
t108 = pkin(5) * t203 - t110;
t107 = t172 * t115 + t132 * t168;
t106 = -t168 * t115 + t132 * t172;
t1 = [MDP(1) + (t113 ^ 2 + t115 ^ 2 + t132 ^ 2) * MDP(20); (-t132 * t166 - t142 * t203) * MDP(10) + (-t133 * t166 + t142 * t204) * MDP(11) + (t132 * t143 + t185 * t203) * MDP(17) + (t123 * t203 + t132 * t144) * MDP(18) + (t113 * t129 - t115 * t128) * MDP(19) + (-t113 * t110 + t115 * t111 + t132 * t130) * MDP(20) + (t106 * t128 + t113 * t124) * MDP(26) + (-t107 * t128 + t113 * t125) * MDP(27) + (t175 * MDP(3) - t171 * MDP(4)) * t164; (t110 ^ 2 + t111 ^ 2 + t130 ^ 2) * MDP(20) + t159 * t170 ^ 2 * MDP(5) + MDP(2) + (t200 * t214 + t202) * t166 + (-0.2e1 * t143 * MDP(13) + t193) * t144 + (t194 - 0.2e1 * t197) * t128 + (0.2e1 * t195 + t196 - 0.2e1 * t198) * t125 + ((0.2e1 * MDP(6) * t170 + t199) * t159 + (MDP(8) * t166 - t184) * t214) * t174 + 0.2e1 * (-t126 * t203 + t138 * t143) * MDP(17) + 0.2e1 * (t127 * t203 + t138 * t144) * MDP(18) + (-t110 * t129 - t111 * t128) * t212 + (t104 * t128 + t108 * t124) * t211 + (-t105 * t128 + t108 * t125) * t210 + 0.2e1 * (t145 * t166 + t159 * t208) * MDP(10) + 0.2e1 * (-t146 * t166 - t159 * t209) * MDP(11); -t133 * MDP(11) + (t113 * t149 - t115 * t148) * MDP(19) + (t113 * t134 + t115 * t136) * MDP(20) + (t106 * t148 + t113 * t206) * MDP(26) + (-t107 * t148 + t113 * t205) * MDP(27) + (t169 * MDP(18) - MDP(10) - t189 + t191) * t132; t202 + t145 * MDP(10) - t146 * MDP(11) + t169 * t193 + (-t169 * t143 + t144 * t173) * MDP(13) + (-pkin(3) * t143 - t138 * t173) * MDP(17) + (-pkin(3) * t144 + t138 * t169) * MDP(18) + (-t136 * t128 + t134 * t129) * MDP(19) + (-t110 * t134 + t111 * t136 + t130 * t158) * MDP(20) + (t118 * t128 + t134 * t124) * MDP(26) + (-t119 * t128 + t134 * t125) * MDP(27) + (-t111 * MDP(19) + t177) * t148 + (-t110 * MDP(19) + (-t125 * MDP(22) - t128 * MDP(24) + t108 * MDP(26)) * t168 + (t108 * MDP(27) + t195 + t196 - t198) * t172) * t149 + (t200 + (MDP(8) - t213) * t174) * t163; MDP(9) + 0.2e1 * pkin(3) * t189 + (t134 ^ 2 + t136 ^ 2 + t158 ^ 2) * MDP(20) + (t161 * MDP(21) - 0.2e1 * t187) * t149 ^ 2 + (MDP(12) * t169 + 0.2e1 * t173 * MDP(13) - 0.2e1 * pkin(3) * MDP(18)) * t169 + (0.2e1 * t179 + t192) * t148 + (t134 * t149 - t136 * t148) * t212 + (t118 * t148 + t134 * t206) * t211 + (-t119 * t148 + t134 * t205) * t210; -t185 * MDP(17) - t123 * MDP(18) - t181 * t113 + (-t113 * t165 + t115 * t162) * MDP(20) * pkin(4); -t163 * t199 + t126 * MDP(17) - t127 * MDP(18) + t125 * t190 + (-t168 * t124 + t125 * t172) * MDP(22) + (-t108 * t172 + t157 * t124) * MDP(26) + (t108 * t168 + t157 * t125) * MDP(27) + t178 * t128 + ((-t128 * t162 - t129 * t165) * MDP(19) + (t110 * t165 + t111 * t162) * MDP(20)) * pkin(4) + t184; -t181 * t134 + t178 * t148 + (t172 * t190 + (-t160 + t161) * MDP(22) + t180 * t157) * t149 + ((-t148 * t162 - t149 * t165) * MDP(19) + (-t134 * t165 + t136 * t162) * MDP(20)) * pkin(4) + t213; 0.2e1 * t187 + t160 * MDP(21) + MDP(16) + (t162 ^ 2 + t165 ^ 2) * MDP(20) * pkin(4) ^ 2 - 0.2e1 * t181 * t157; t132 * MDP(20); t130 * MDP(20) + t181 * t128; t181 * t148 + t191; 0; MDP(20); t106 * MDP(26) - t107 * MDP(27); t177; t118 * MDP(26) - t119 * MDP(27) + t179 + t192; t178; t181; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
