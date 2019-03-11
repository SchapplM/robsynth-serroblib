% Calculate joint inertia matrix for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRRPP1_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:46:55
% EndTime: 2019-03-08 22:46:57
% DurationCPUTime: 0.63s
% Computational Cost: add. (905->188), mult. (1862->275), div. (0->0), fcn. (1988->10), ass. (0->74)
t187 = MDP(19) + MDP(22);
t165 = MDP(20) + MDP(24);
t143 = cos(qJ(4));
t178 = -qJ(5) - pkin(9);
t125 = t178 * t143;
t137 = sin(pkin(11));
t139 = cos(pkin(11));
t140 = sin(qJ(4));
t161 = t178 * t140;
t110 = -t125 * t137 - t139 * t161;
t112 = -t139 * t125 + t137 * t161;
t171 = MDP(17) * t140;
t186 = t140 * MDP(14) + t143 * MDP(15) - t110 * MDP(21) + t112 * MDP(23) - (MDP(18) * t143 + t171) * pkin(9);
t185 = 2 * MDP(19);
t184 = 0.2e1 * MDP(21);
t183 = 2 * MDP(22);
t182 = 0.2e1 * MDP(23);
t181 = pkin(8) * t140;
t180 = pkin(8) * t143;
t144 = cos(qJ(3));
t179 = pkin(8) * t144;
t177 = cos(pkin(6));
t138 = sin(pkin(6));
t142 = sin(qJ(2));
t176 = t138 * t142;
t145 = cos(qJ(2));
t175 = t138 * t145;
t141 = sin(qJ(3));
t174 = t140 * t141;
t173 = t140 * t143;
t172 = t141 * t143;
t124 = -pkin(3) * t144 - pkin(9) * t141 - pkin(2);
t121 = t143 * t124;
t105 = -qJ(5) * t172 + t121 + (-pkin(4) - t181) * t144;
t163 = t143 * t179;
t108 = t163 + (-qJ(5) * t141 + t124) * t140;
t93 = t137 * t105 + t139 * t108;
t123 = pkin(4) * t174 + t141 * pkin(8);
t170 = MDP(17) * t143;
t119 = t137 * t140 - t139 * t143;
t120 = t137 * t143 + t139 * t140;
t132 = -pkin(4) * t143 - pkin(3);
t101 = pkin(5) * t119 - qJ(6) * t120 + t132;
t169 = t101 * MDP(24);
t168 = t119 * MDP(21);
t167 = t120 * MDP(23);
t166 = t141 * MDP(11);
t164 = t110 ^ 2 + t112 ^ 2;
t162 = MDP(13) * t173;
t156 = t177 * t144;
t115 = t120 * t141;
t116 = -t137 * t174 + t139 * t172;
t155 = t110 * t116 - t112 * t115;
t92 = t139 * t105 - t137 * t108;
t153 = MDP(14) * t143 - MDP(15) * t140;
t152 = -t140 * MDP(18) + t170;
t150 = t141 * t177 + t144 * t176;
t148 = MDP(20) * t132 - t167 + t168 + t169;
t147 = -t140 * t150 - t143 * t175;
t136 = t143 ^ 2;
t135 = t141 ^ 2;
t134 = t140 ^ 2;
t130 = pkin(4) * t139 + pkin(5);
t128 = pkin(4) * t137 + qJ(6);
t118 = t141 * t176 - t156;
t114 = t124 * t140 + t163;
t113 = -t140 * t179 + t121;
t107 = -t140 * t175 + t143 * t150;
t98 = pkin(5) * t115 - qJ(6) * t116 + t123;
t97 = t139 * t107 + t137 * t147;
t95 = t107 * t137 - t139 * t147;
t91 = pkin(5) * t144 - t92;
t90 = -qJ(6) * t144 + t93;
t1 = [MDP(1) + t165 * (t118 ^ 2 + t95 ^ 2 + t97 ^ 2); (t107 * t144 + t118 * t172) * MDP(18) + (t118 * t123 - t92 * t95 + t93 * t97) * MDP(20) + (t118 * t115 + t95 * t144) * MDP(21) + (-t116 * t118 - t97 * t144) * MDP(23) + (t118 * t98 + t90 * t97 + t95 * t91) * MDP(24) + (t156 + t118) * t141 * t171 + ((t144 ^ 2 * t171 - MDP(4)) * t142 + (-t166 + MDP(3) + (MDP(10) + t170) * t144) * t145) * t138 + t187 * (-t97 * t115 + t95 * t116); MDP(2) - 0.2e1 * pkin(2) * t166 + (t123 ^ 2 + t92 ^ 2 + t93 ^ 2) * MDP(20) + (t90 ^ 2 + t91 ^ 2 + t98 ^ 2) * MDP(24) + (t136 * MDP(12) + MDP(5) - 0.2e1 * t162) * t135 + (0.2e1 * pkin(2) * MDP(10) + t144 * MDP(16) + 0.2e1 * (MDP(6) - t153) * t141) * t144 + 0.2e1 * (-t113 * t144 + t135 * t181) * MDP(17) + 0.2e1 * (t114 * t144 + t135 * t180) * MDP(18) + (-t115 * t93 - t116 * t92) * t185 + (t115 * t98 + t144 * t91) * t184 + (-t115 * t90 + t116 * t91) * t183 + (-t116 * t98 - t144 * t90) * t182; -t150 * MDP(11) + (-MDP(10) + t148 - t152) * t118 + t165 * (t110 * t95 + t97 * t112) + t187 * (-t97 * t119 + t120 * t95); (-t119 * t93 - t120 * t92 + t155) * MDP(19) + (-t110 * t92 + t112 * t93 + t123 * t132) * MDP(20) + (t101 * t115 + t119 * t98) * MDP(21) + (-t119 * t90 + t120 * t91 + t155) * MDP(22) + (-t101 * t116 - t120 * t98) * MDP(23) + (t101 * t98 + t110 * t91 + t112 * t90) * MDP(24) + (-pkin(8) * MDP(11) + MDP(8) - t186) * t144 + (MDP(7) - pkin(8) * MDP(10) + MDP(12) * t173 + (-t134 + t136) * MDP(13) + (-pkin(3) * t140 - t180) * MDP(17) + (-pkin(3) * t143 + t181) * MDP(18)) * t141; MDP(9) + t134 * MDP(12) + 0.2e1 * t162 + (t132 ^ 2 + t164) * MDP(20) + t164 * MDP(24) + (-0.2e1 * t167 + 0.2e1 * t168 + t169) * t101 + 0.2e1 * t152 * pkin(3) + (t185 + t183) * (t110 * t120 - t112 * t119); t147 * MDP(17) - t107 * MDP(18) - t95 * MDP(21) + t97 * MDP(23) + (t128 * t97 - t130 * t95) * MDP(24) + (t137 * t97 - t139 * t95) * MDP(20) * pkin(4); t113 * MDP(17) - t114 * MDP(18) + t92 * MDP(21) + (-t128 * t115 - t130 * t116) * MDP(22) + t93 * MDP(23) + (t90 * t128 - t91 * t130) * MDP(24) + t153 * t141 + (-MDP(16) + (-pkin(5) - t130) * MDP(21) + (-qJ(6) - t128) * MDP(23)) * t144 + ((-t115 * t137 - t116 * t139) * MDP(19) + (t137 * t93 + t139 * t92) * MDP(20)) * pkin(4); (-t119 * t128 - t120 * t130) * MDP(22) + (-t110 * t130 + t112 * t128) * MDP(24) + ((-t119 * t137 - t120 * t139) * MDP(19) + (-t110 * t139 + t112 * t137) * MDP(20)) * pkin(4) + t186; MDP(16) + (t128 ^ 2 + t130 ^ 2) * MDP(24) + (t137 ^ 2 + t139 ^ 2) * MDP(20) * pkin(4) ^ 2 + t130 * t184 + t128 * t182; t165 * t118; t123 * MDP(20) + t115 * MDP(21) - t116 * MDP(23) + t98 * MDP(24); t148; 0; t165; t95 * MDP(24); t144 * MDP(21) + t116 * MDP(22) + t91 * MDP(24); t120 * MDP(22) + t110 * MDP(24); -MDP(24) * t130 - MDP(21); 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
