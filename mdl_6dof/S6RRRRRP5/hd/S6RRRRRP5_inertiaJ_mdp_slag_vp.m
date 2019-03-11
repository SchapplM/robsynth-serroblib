% Calculate joint inertia matrix for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP5_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:23:06
% EndTime: 2019-03-10 01:23:10
% DurationCPUTime: 1.03s
% Computational Cost: add. (1501->206), mult. (2951->289), div. (0->0), fcn. (3166->8), ass. (0->93)
t165 = sin(qJ(3));
t168 = cos(qJ(3));
t206 = pkin(8) + pkin(9);
t150 = t206 * t165;
t151 = t206 * t168;
t164 = sin(qJ(4));
t167 = cos(qJ(4));
t128 = -t167 * t150 - t151 * t164;
t129 = -t150 * t164 + t151 * t167;
t145 = t164 * t165 - t167 * t168;
t146 = t164 * t168 + t165 * t167;
t116 = -pkin(10) * t146 + t128;
t117 = -pkin(10) * t145 + t129;
t163 = sin(qJ(5));
t205 = cos(qJ(5));
t104 = t205 * t116 - t117 * t163;
t105 = t163 * t116 + t205 * t117;
t121 = t205 * t145 + t146 * t163;
t122 = -t163 * t145 + t205 * t146;
t177 = t122 * MDP(27) - t121 * MDP(28) + t104 * MDP(30) - t105 * MDP(31);
t170 = t146 * MDP(20) - t145 * MDP(21) + t128 * MDP(23) - t129 * MDP(24) + t177;
t218 = t170 - (t165 * MDP(16) + t168 * MDP(17)) * pkin(8) + t165 * MDP(13) + t168 * MDP(14);
t217 = (t167 * MDP(23) - MDP(24) * t164) * pkin(3);
t166 = sin(qJ(2));
t136 = t146 * t166;
t137 = t145 * t166;
t113 = t205 * t136 - t137 * t163;
t114 = -t163 * t136 - t205 * t137;
t193 = t114 * MDP(27) - t113 * MDP(28);
t216 = t137 * MDP(20) + t136 * MDP(21);
t169 = cos(qJ(2));
t149 = -pkin(2) * t169 - t166 * pkin(8) - pkin(1);
t144 = t168 * t149;
t199 = pkin(9) * t166;
t202 = pkin(7) * t165;
t123 = -t168 * t199 + t144 + (-pkin(3) - t202) * t169;
t200 = pkin(7) * t169;
t183 = t168 * t200;
t125 = t183 + (t149 - t199) * t165;
t108 = t167 * t123 - t164 * t125;
t203 = pkin(4) * t169;
t101 = t137 * pkin(10) + t108 - t203;
t109 = t164 * t123 + t167 * t125;
t106 = -pkin(10) * t136 + t109;
t96 = t205 * t101 - t163 * t106;
t182 = t205 * t106;
t97 = t163 * t101 + t182;
t215 = t96 * MDP(30) - t97 * MDP(31);
t212 = t108 * MDP(23) - t109 * MDP(24) - t216;
t211 = -2 * MDP(19);
t210 = 0.2e1 * MDP(23);
t209 = 0.2e1 * MDP(24);
t208 = -2 * MDP(26);
t207 = 2 * MDP(32);
t204 = pkin(3) * t164;
t201 = pkin(7) * t168;
t198 = t163 * pkin(4);
t197 = MDP(33) * pkin(5);
t196 = t165 * t168;
t148 = (pkin(3) * t165 + pkin(7)) * t166;
t189 = t114 * MDP(25);
t155 = pkin(3) * t167 + pkin(4);
t139 = t205 * t155 - t163 * t204;
t188 = t139 * MDP(30);
t178 = t205 * t204;
t140 = t163 * t155 + t178;
t187 = t140 * MDP(31);
t186 = t146 * MDP(18);
t184 = MDP(22) + MDP(29);
t158 = t205 * pkin(4);
t156 = -pkin(3) * t168 - pkin(2);
t181 = MDP(12) * t196;
t180 = MDP(30) * t205;
t179 = MDP(15) + t184;
t124 = pkin(4) * t136 + t148;
t176 = -t169 * MDP(29) + t193;
t133 = pkin(4) * t145 + t156;
t175 = MDP(13) * t168 - MDP(14) * t165;
t171 = -MDP(29) + t187 - t188;
t161 = t168 ^ 2;
t160 = t166 ^ 2;
t159 = t165 ^ 2;
t154 = t158 + pkin(5);
t138 = pkin(5) + t139;
t132 = t165 * t149 + t183;
t131 = -t165 * t200 + t144;
t110 = pkin(5) * t121 + t133;
t107 = pkin(5) * t113 + t124;
t99 = -t121 * qJ(6) + t105;
t98 = -qJ(6) * t122 + t104;
t95 = -t113 * qJ(6) + t97;
t94 = -pkin(5) * t169 - t114 * qJ(6) + t96;
t1 = [-0.2e1 * pkin(1) * t166 * MDP(10) + (t107 ^ 2 + t94 ^ 2 + t95 ^ 2) * MDP(33) + MDP(1) - (-t137 * MDP(18) + t136 * t211) * t137 + (t113 * t208 + t189) * t114 + t179 * t169 ^ 2 + (t161 * MDP(11) + MDP(4) - 0.2e1 * t181) * t160 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t175) * t166 + t216 - t193) * t169 + 0.2e1 * (t132 * t169 + t160 * t201) * MDP(17) + 0.2e1 * (-t131 * t169 + t160 * t202) * MDP(16) + (-t108 * t169 + t148 * t136) * t210 + (t109 * t169 - t148 * t137) * t209 + 0.2e1 * (t124 * t114 + t169 * t97) * MDP(31) + 0.2e1 * (t124 * t113 - t169 * t96) * MDP(30) + (-t113 * t95 - t114 * t94) * t207; (-t136 * t146 + t137 * t145) * MDP(19) + (t156 * t136 + t148 * t145) * MDP(23) + (-t156 * t137 + t148 * t146) * MDP(24) + (-t113 * t122 - t114 * t121) * MDP(26) + (t133 * t113 + t124 * t121) * MDP(30) + (t133 * t114 + t124 * t122) * MDP(31) + (-t113 * t99 - t114 * t98 - t121 * t95 - t122 * t94) * MDP(32) + (t107 * t110 + t94 * t98 + t95 * t99) * MDP(33) - t137 * t186 + t122 * t189 + (MDP(6) + (-t159 + t161) * MDP(12) + (-pkin(2) * t165 - t201) * MDP(16) + (-pkin(2) * t168 + t202) * MDP(17) - pkin(7) * MDP(9) + MDP(11) * t196) * t166 + (-pkin(7) * MDP(10) + MDP(7) - t218) * t169; MDP(8) + t159 * MDP(11) + 0.2e1 * t181 + t156 * t145 * t210 - t121 * t99 * t207 + (t110 ^ 2 + t98 ^ 2 + t99 ^ 2) * MDP(33) + (t145 * t211 + t156 * t209 + t186) * t146 + (MDP(25) * t122 + t121 * t208 - t98 * t207) * t122 + 0.2e1 * (t121 * MDP(30) + t122 * MDP(31)) * t133 + 0.2e1 * (t168 * MDP(16) - t165 * MDP(17)) * pkin(2); t131 * MDP(16) - t132 * MDP(17) + (-t113 * t140 - t114 * t138) * MDP(32) + (t138 * t94 + t140 * t95) * MDP(33) + t175 * t166 + (-MDP(15) - MDP(22) + t171 - t217) * t169 + t193 + t212 + t215; (-t121 * t140 - t122 * t138) * MDP(32) + (t138 * t98 + t140 * t99) * MDP(33) + t218; (t138 ^ 2 + t140 ^ 2) * MDP(33) + 0.2e1 * t217 + 0.2e1 * t188 - 0.2e1 * t187 + t179; -t169 * MDP(22) + (-t169 * t158 + t96) * MDP(30) + (-t182 + (-t101 + t203) * t163) * MDP(31) + (-t113 * t198 - t114 * t154) * MDP(32) + (t154 * t94 + t95 * t198) * MDP(33) + t176 + t212; (-t121 * t198 - t122 * t154) * MDP(32) + (t154 * t98 + t99 * t198) * MDP(33) + t170; (t139 + t158) * MDP(30) + (-t178 + (-pkin(4) - t155) * t163) * MDP(31) + (t138 * t154 + t140 * t198) * MDP(33) + t184 + t217; t154 ^ 2 * MDP(33) + (0.2e1 * t180 + (MDP(33) * t198 - 0.2e1 * MDP(31)) * t163) * pkin(4) + t184; (-t114 * MDP(32) + t94 * MDP(33)) * pkin(5) + t176 + t215; (-t122 * MDP(32) + t98 * MDP(33)) * pkin(5) + t177; t138 * t197 - t171; t154 * t197 + MDP(29) + (-MDP(31) * t163 + t180) * pkin(4); MDP(33) * pkin(5) ^ 2 + MDP(29); t107 * MDP(33); t110 * MDP(33); 0; 0; 0; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
