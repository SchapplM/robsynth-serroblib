% Calculate joint inertia matrix for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRPR3_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:19:09
% EndTime: 2019-03-09 10:19:11
% DurationCPUTime: 0.70s
% Computational Cost: add. (1409->164), mult. (2652->260), div. (0->0), fcn. (3070->10), ass. (0->80)
t164 = sin(pkin(10));
t157 = pkin(2) * t164 + pkin(8);
t168 = sin(qJ(4));
t171 = cos(qJ(4));
t177 = t168 * MDP(18) + t171 * MDP(19);
t188 = qJ(5) + t157;
t144 = t188 * t168;
t145 = t188 * t171;
t163 = sin(pkin(11));
t165 = cos(pkin(11));
t122 = -t165 * t144 - t145 * t163;
t150 = t163 * t171 + t165 * t168;
t115 = -pkin(9) * t150 + t122;
t123 = -t163 * t144 + t165 * t145;
t148 = -t163 * t168 + t165 * t171;
t116 = pkin(9) * t148 + t123;
t167 = sin(qJ(6));
t170 = cos(qJ(6));
t129 = -t170 * t148 + t150 * t167;
t130 = t148 * t167 + t150 * t170;
t179 = t130 * MDP(24) - t129 * MDP(25) + (t115 * t170 - t116 * t167) * MDP(27) - (t115 * t167 + t116 * t170) * MDP(28);
t203 = t168 * MDP(15) + t171 * MDP(16) - t157 * t177 + t179;
t178 = t171 * MDP(18) - t168 * MDP(19);
t124 = t129 * MDP(27);
t187 = -t130 * MDP(28) - t124;
t202 = t187 + t178;
t201 = -qJ(3) - pkin(7);
t166 = cos(pkin(10));
t169 = sin(qJ(2));
t196 = cos(qJ(2));
t151 = t164 * t196 + t166 * t169;
t119 = t150 * t151;
t120 = t148 * t151;
t110 = t170 * t119 + t120 * t167;
t111 = -t119 * t167 + t120 * t170;
t200 = t111 * MDP(24) - t110 * MDP(25);
t199 = 2 * MDP(20);
t198 = -2 * MDP(23);
t197 = 0.2e1 * MDP(28);
t195 = pkin(4) * t163;
t154 = t201 * t169;
t155 = t201 * t196;
t133 = t154 * t164 - t155 * t166;
t194 = t133 * t171;
t193 = t151 * t168;
t192 = t151 * t171;
t189 = t168 * t171;
t149 = t164 * t169 - t166 * t196;
t160 = -pkin(2) * t196 - pkin(1);
t128 = t149 * pkin(3) - t151 * pkin(8) + t160;
t113 = t171 * t128 - t133 * t168;
t105 = pkin(4) * t149 - qJ(5) * t192 + t113;
t109 = t194 + (-qJ(5) * t151 + t128) * t168;
t98 = t163 * t105 + t165 * t109;
t186 = t130 * MDP(22);
t158 = pkin(4) * t165 + pkin(5);
t140 = t158 * t170 - t167 * t195;
t185 = t140 * MDP(27);
t141 = t158 * t167 + t170 * t195;
t184 = t141 * MDP(28);
t183 = 0.2e1 * t196;
t182 = MDP(17) + MDP(26);
t181 = t149 * MDP(26) + t200;
t159 = -pkin(2) * t166 - pkin(3);
t180 = MDP(14) * t189;
t97 = t165 * t105 - t109 * t163;
t95 = pkin(5) * t149 - pkin(9) * t120 + t97;
t96 = -pkin(9) * t119 + t98;
t92 = -t167 * t96 + t170 * t95;
t131 = -t166 * t154 - t155 * t164;
t117 = pkin(4) * t193 + t131;
t93 = t167 * t95 + t170 * t96;
t153 = -pkin(4) * t171 + t159;
t175 = (t171 * MDP(15) - t168 * MDP(16)) * t151;
t162 = t171 ^ 2;
t161 = t168 ^ 2;
t134 = -pkin(5) * t148 + t153;
t114 = t128 * t168 + t194;
t112 = pkin(5) * t119 + t117;
t1 = [MDP(1) + pkin(1) * MDP(9) * t183 + (t131 ^ 2 + t133 ^ 2 + t160 ^ 2) * MDP(12) + (t117 ^ 2 + t97 ^ 2 + t98 ^ 2) * MDP(21) + (t162 * MDP(13) - 0.2e1 * t180) * t151 ^ 2 + t182 * t149 ^ 2 + (t111 * MDP(22) + t110 * t198) * t111 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t169 + MDP(5) * t183) * t169 + 0.2e1 * (t175 + t200) * t149 + 0.2e1 * (t131 * t151 - t133 * t149) * MDP(11) + 0.2e1 * (t113 * t149 + t131 * t193) * MDP(18) + 0.2e1 * (-t114 * t149 + t131 * t192) * MDP(19) + (-t119 * t98 - t120 * t97) * t199 + 0.2e1 * (t110 * t112 + t149 * t92) * MDP(27) + (t111 * t112 - t149 * t93) * t197; t169 * MDP(6) + t196 * MDP(7) + (-t119 * t123 - t120 * t122 + t148 * t98 - t150 * t97) * MDP(20) + (t117 * t153 + t122 * t97 + t123 * t98) * MDP(21) + t111 * t186 + (-t110 * t130 - t111 * t129) * MDP(23) + (t110 * t134 + t112 * t129) * MDP(27) + (t111 * t134 + t112 * t130) * MDP(28) - t178 * t131 + (-MDP(10) * t196 - t169 * MDP(9)) * pkin(7) + (MDP(13) * t189 + (-t161 + t162) * MDP(14) + t177 * t159) * t151 + t203 * t149 + ((-t149 * t164 - t151 * t166) * MDP(11) + (-t131 * t166 + t133 * t164) * MDP(12)) * pkin(2); MDP(8) + t161 * MDP(13) + 0.2e1 * t180 + (-t122 * t150 + t123 * t148) * t199 + (t122 ^ 2 + t123 ^ 2 + t153 ^ 2) * MDP(21) + 0.2e1 * t134 * t124 + (t164 ^ 2 + t166 ^ 2) * MDP(12) * pkin(2) ^ 2 - 0.2e1 * t178 * t159 + (t129 * t198 + t134 * t197 + t186) * t130; t160 * MDP(12) + (-t119 * t150 - t120 * t148) * MDP(20) + (t148 * t97 + t150 * t98) * MDP(21) + t202 * t149; (t122 * t148 + t123 * t150) * MDP(21); MDP(12) + (t148 ^ 2 + t150 ^ 2) * MDP(21); t149 * MDP(17) + t113 * MDP(18) - t114 * MDP(19) + (t140 * t149 + t92) * MDP(27) + (-t141 * t149 - t93) * MDP(28) + t175 + ((-t119 * t163 - t120 * t165) * MDP(20) + (t163 * t98 + t165 * t97) * MDP(21)) * pkin(4) + t181; ((t148 * t163 - t150 * t165) * MDP(20) + (t122 * t165 + t123 * t163) * MDP(21)) * pkin(4) + t203; (t148 * t165 + t150 * t163) * MDP(21) * pkin(4) + t202; (t163 ^ 2 + t165 ^ 2) * MDP(21) * pkin(4) ^ 2 + 0.2e1 * t185 - 0.2e1 * t184 + t182; MDP(21) * t117 + t110 * MDP(27) + t111 * MDP(28); MDP(21) * t153 - t187; 0; 0; MDP(21); t92 * MDP(27) - t93 * MDP(28) + t181; t179; t187; MDP(26) - t184 + t185; 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
