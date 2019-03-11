% Calculate joint inertia matrix for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRPRR3_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:07:45
% EndTime: 2019-03-08 22:07:48
% DurationCPUTime: 0.79s
% Computational Cost: add. (1024->195), mult. (2660->309), div. (0->0), fcn. (3050->14), ass. (0->100)
t157 = sin(qJ(6));
t161 = cos(qJ(6));
t172 = MDP(26) * t157 + MDP(27) * t161;
t167 = t157 * MDP(23) + t161 * MDP(24) - t172 * pkin(11);
t210 = 0.2e1 * MDP(26);
t209 = 0.2e1 * MDP(27);
t163 = cos(qJ(3));
t208 = pkin(2) * t163;
t207 = pkin(9) + qJ(4);
t155 = cos(pkin(7));
t143 = t155 * t208;
t152 = sin(pkin(7));
t159 = sin(qJ(3));
t200 = t152 * t159;
t129 = t155 * pkin(3) - t207 * t200 + t143;
t181 = t155 * t159 * pkin(2);
t199 = t152 * t163;
t132 = t207 * t199 + t181;
t151 = sin(pkin(13));
t154 = cos(pkin(13));
t118 = t151 * t129 + t154 * t132;
t116 = t155 * pkin(10) + t118;
t135 = t151 * t200 - t154 * t199;
t136 = (t151 * t163 + t154 * t159) * t152;
t141 = (-pkin(3) * t163 - pkin(2)) * t152;
t121 = t135 * pkin(4) - t136 * pkin(10) + t141;
t158 = sin(qJ(5));
t162 = cos(qJ(5));
t108 = -t158 * t116 + t162 * t121;
t105 = -t135 * pkin(5) - t108;
t206 = t105 * t157;
t205 = t105 * t161;
t145 = t151 * pkin(3) + pkin(10);
t204 = t145 * t157;
t203 = t145 * t161;
t202 = t145 * t162;
t147 = t152 ^ 2;
t201 = t147 * t159;
t153 = sin(pkin(6));
t164 = cos(qJ(2));
t198 = t153 * t164;
t197 = t155 * MDP(9);
t196 = t155 * t164;
t195 = t157 * t158;
t194 = t158 * t161;
t160 = sin(qJ(2));
t193 = t159 * t160;
t128 = t162 * t136 + t158 * t155;
t119 = t157 * t128 - t161 * t135;
t192 = t119 * MDP(24);
t120 = t161 * t128 + t157 * t135;
t191 = t120 * MDP(21);
t190 = t120 * MDP(23);
t127 = t158 * t136 - t162 * t155;
t189 = t127 * MDP(25);
t188 = t128 * MDP(14);
t187 = t128 * MDP(15);
t146 = -t154 * pkin(3) - pkin(4);
t186 = t146 * MDP(19);
t185 = t146 * MDP(20);
t184 = t157 * MDP(21);
t183 = t158 * MDP(20);
t182 = t162 * MDP(25);
t180 = t157 * t161 * MDP(22);
t117 = t154 * t129 - t151 * t132;
t179 = -t145 * MDP(19) + MDP(16);
t178 = -t145 * MDP(20) + MDP(17);
t109 = t162 * t116 + t158 * t121;
t177 = t162 * MDP(19) - t183;
t176 = MDP(23) * t161 - MDP(24) * t157;
t140 = -t162 * pkin(5) - t158 * pkin(11) + t146;
t125 = t161 * t140 - t157 * t202;
t126 = t157 * t140 + t161 * t202;
t174 = t125 * MDP(26) - t126 * MDP(27);
t173 = t161 * MDP(26) - t157 * MDP(27);
t115 = -t155 * pkin(4) - t117;
t171 = -MDP(15) + t176;
t170 = MDP(19) + t173;
t169 = (t159 * MDP(7) + t163 * MDP(8)) * t152;
t156 = cos(pkin(6));
t168 = t156 * t199 + (t163 * t196 - t193) * t153;
t106 = t135 * pkin(11) + t109;
t107 = t127 * pkin(5) - t128 * pkin(11) + t115;
t101 = -t157 * t106 + t161 * t107;
t102 = t161 * t106 + t157 * t107;
t166 = t101 * MDP(26) - t102 * MDP(27) + t189 + t190 - t192;
t150 = t161 ^ 2;
t149 = t158 ^ 2;
t148 = t157 ^ 2;
t139 = pkin(9) * t199 + t181;
t138 = -pkin(9) * t200 + t143;
t137 = -t152 * t198 + t156 * t155;
t123 = t156 * t200 + (t159 * t196 + t160 * t163) * t153;
t114 = t154 * t123 + t151 * t168;
t112 = t151 * t123 - t154 * t168;
t111 = t162 * t114 + t137 * t158;
t110 = t158 * t114 - t137 * t162;
t104 = t161 * t111 + t157 * t112;
t103 = -t157 * t111 + t161 * t112;
t1 = [MDP(1) + (t112 ^ 2 + t114 ^ 2 + t137 ^ 2) * MDP(13); MDP(3) * t198 - t153 * t160 * MDP(4) + (-t153 * t155 * t193 + ((t152 * t156 + t153 * t196) * t155 - t137 * t152) * t163) * MDP(10) + (-t123 * t155 + t137 * t200) * MDP(11) + (t112 * t136 - t114 * t135) * MDP(12) + (-t112 * t117 + t114 * t118 + t137 * t141) * MDP(13) + (-t110 * t135 + t112 * t127) * MDP(19) + (-t111 * t135 + t112 * t128) * MDP(20) + (t103 * t127 + t110 * t119) * MDP(26) + (-t104 * t127 + t110 * t120) * MDP(27); (t117 ^ 2 + t118 ^ 2 + t141 ^ 2) * MDP(13) + t135 ^ 2 * MDP(18) + MDP(2) + (MDP(5) * t159 + 0.2e1 * MDP(6) * t163) * t201 + (0.2e1 * t135 * MDP(16) + t188) * t128 + (-0.2e1 * t119 * MDP(22) + t191) * t120 + (0.2e1 * t169 + t197) * t155 + (-0.2e1 * t135 * MDP(17) - 0.2e1 * t187 + t189 + 0.2e1 * t190 - 0.2e1 * t192) * t127 + 0.2e1 * (t138 * t155 + t147 * t208) * MDP(10) + 0.2e1 * (-pkin(2) * t201 - t139 * t155) * MDP(11) + 0.2e1 * (-t117 * t136 - t118 * t135) * MDP(12) + 0.2e1 * (t108 * t135 + t115 * t127) * MDP(19) + 0.2e1 * (-t109 * t135 + t115 * t128) * MDP(20) + (t101 * t127 + t105 * t119) * t210 + (-t102 * t127 + t105 * t120) * t209; t168 * MDP(10) - t123 * MDP(11) + (-t103 * t162 + t110 * t195) * MDP(26) + (t104 * t162 + t110 * t194) * MDP(27) - t177 * t112 + (-t112 * t154 + t114 * t151) * MDP(13) * pkin(3); t128 * t185 + t138 * MDP(10) - t139 * MDP(11) + t197 + t169 + (t174 + t186) * t127 + ((-t135 * t151 - t136 * t154) * MDP(12) + (t117 * t154 + t118 * t151) * MDP(13)) * pkin(3) + (-t115 * MDP(19) + t178 * t135 - t166 + t187) * t162 + (t188 + t115 * MDP(20) + t161 * t191 + (-t119 * t161 - t120 * t157) * MDP(22) + (t119 * t145 + t206) * MDP(26) + (t120 * t145 + t205) * MDP(27) + t179 * t135 + t171 * t127) * t158; MDP(9) + (t151 ^ 2 + t154 ^ 2) * MDP(13) * pkin(3) ^ 2 + (t182 - 0.2e1 * t186) * t162 + (t150 * MDP(21) + MDP(14) - 0.2e1 * t180) * t149 + (-t125 * t162 + t149 * t204) * t210 + (t126 * t162 + t149 * t203) * t209 + 0.2e1 * (-t162 * t171 + t185) * t158; t137 * MDP(13); t141 * MDP(13) + (-t162 * t119 - t127 * t195) * MDP(26) + (-t120 * t162 - t127 * t194) * MDP(27) + t177 * t135; 0; MDP(13); -t111 * MDP(20) - t170 * t110; t128 * MDP(16) + t135 * MDP(18) + t108 * MDP(19) - t109 * MDP(20) + t120 * t184 + (-t157 * t119 + t120 * t161) * MDP(22) + (-pkin(5) * t119 - t205) * MDP(26) + (-pkin(5) * t120 + t206) * MDP(27) + (-MDP(17) + t167) * t127; (-t167 + t178) * t162 + (t161 * t184 + (-t148 + t150) * MDP(22) + (-pkin(5) * t157 - t203) * MDP(26) + (-pkin(5) * t161 + t204) * MDP(27) + t179) * t158; t170 * t162 - t183; t148 * MDP(21) + 0.2e1 * pkin(5) * t173 + MDP(18) + 0.2e1 * t180; t103 * MDP(26) - t104 * MDP(27); t166; t176 * t158 + t174 - t182; -t172 * t158; t167; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
