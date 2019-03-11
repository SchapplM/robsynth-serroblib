% Calculate joint inertia matrix for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR6_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:15:17
% EndTime: 2019-03-09 09:15:18
% DurationCPUTime: 0.48s
% Computational Cost: add. (696->145), mult. (1206->207), div. (0->0), fcn. (1311->8), ass. (0->69)
t129 = sin(qJ(6));
t132 = cos(qJ(6));
t171 = MDP(31) * t129 + MDP(32) * t132;
t127 = sin(pkin(10));
t128 = cos(pkin(10));
t130 = sin(qJ(5));
t133 = cos(qJ(5));
t103 = t127 * t130 - t128 * t133;
t106 = t127 * t133 + t128 * t130;
t155 = MDP(31) * t132;
t141 = -MDP(32) * t129 + t155;
t139 = MDP(24) + t141;
t170 = -t106 * MDP(25) - t139 * t103;
t131 = sin(qJ(2));
t134 = cos(qJ(2));
t114 = -t134 * pkin(2) - t131 * qJ(3) - pkin(1);
t102 = t134 * pkin(3) - t114;
t104 = t127 * t131 + t128 * t134;
t97 = pkin(4) * t104 + t102;
t169 = 0.2e1 * t97;
t168 = -pkin(2) - pkin(3);
t111 = qJ(3) * t127 - t128 * t168;
t110 = -pkin(4) - t111;
t113 = t128 * qJ(3) + t127 * t168;
t95 = -t110 * t133 + t113 * t130;
t93 = pkin(5) + t95;
t167 = pkin(5) + t93;
t166 = pkin(7) - qJ(4);
t105 = -t134 * t127 + t128 * t131;
t92 = -t104 * t130 + t105 * t133;
t165 = t132 * t92;
t164 = MDP(25) * t92;
t91 = t133 * t104 + t105 * t130;
t163 = MDP(30) * t91;
t162 = t95 * MDP(24);
t96 = t110 * t130 + t113 * t133;
t161 = t96 * MDP(25);
t115 = t166 * t134;
t148 = t166 * t131;
t100 = t128 * t115 + t127 * t148;
t124 = t131 ^ 2;
t160 = t134 ^ 2 + t124;
t159 = MDP(15) * t104;
t158 = MDP(16) * t105;
t157 = MDP(18) * t102;
t152 = t132 * MDP(27);
t123 = t129 ^ 2;
t151 = t123 * MDP(26) + MDP(23);
t150 = t129 * t152;
t149 = 0.2e1 * t150 + t151;
t147 = -pkin(2) * MDP(14) - MDP(11);
t98 = t115 * t127 - t128 * t148;
t146 = -pkin(5) * t92 - pkin(9) * t91;
t94 = -pkin(9) + t96;
t145 = -t91 * t94 + t92 * t93;
t88 = -pkin(8) * t105 - t98;
t89 = -pkin(8) * t104 + t100;
t85 = t130 * t89 - t133 * t88;
t143 = -t91 * MDP(29) + t85 * MDP(31);
t142 = MDP(28) * t132 - MDP(29) * t129;
t138 = 0.2e1 * t141;
t137 = -MDP(26) * t165 - t91 * MDP(28) - t85 * MDP(32);
t125 = t132 ^ 2;
t86 = t130 * t88 + t133 * t89;
t136 = -t91 * MDP(22) - t85 * MDP(24) - t86 * MDP(25) + (MDP(21) - (t123 - t125) * MDP(27)) * t92;
t84 = pkin(5) * t91 - pkin(9) * t92 + t97;
t83 = t129 * t84 + t132 * t86;
t82 = -t129 * t86 + t132 * t84;
t1 = [t124 * MDP(4) + (t160 * pkin(7) ^ 2 + t114 ^ 2) * MDP(14) + (t100 ^ 2 + t98 ^ 2) * MDP(18) + t164 * t169 + MDP(1) + (t157 + 0.2e1 * t158 + 0.2e1 * t159) * t102 + (MDP(26) * t125 + MDP(19) - 0.2e1 * t150) * t92 ^ 2 + (MDP(24) * t169 + t163) * t91 + 0.2e1 * (-t100 * t104 + t105 * t98) * MDP(17) + 0.2e1 * (t129 * t85 * t92 + t82 * t91) * MDP(31) + 0.2e1 * (t85 * t165 - t83 * t91) * MDP(32) + 0.2e1 * t160 * MDP(12) * pkin(7) + 0.2e1 * (-MDP(11) * t114 + MDP(9) * pkin(1)) * t134 + 0.2e1 * (-MDP(10) * pkin(1) - MDP(13) * t114 + MDP(5) * t134) * t131 + 0.2e1 * (-MDP(20) + t142) * t92 * t91; t131 * MDP(6) + t134 * MDP(7) + (-pkin(2) * t131 + qJ(3) * t134) * MDP(12) + t98 * MDP(15) + t100 * MDP(16) + (-t104 * t113 + t105 * t111) * MDP(17) + (t100 * t113 + t111 * t98) * MDP(18) + (t145 * MDP(32) + t143) * t132 + (t145 * MDP(31) + t137) * t129 + ((MDP(14) * qJ(3) - MDP(10) + MDP(13)) * t134 + (-MDP(9) + t147) * t131) * pkin(7) - t136; MDP(8) + 0.2e1 * pkin(2) * MDP(11) + 0.2e1 * qJ(3) * MDP(13) + (pkin(2) ^ 2 + qJ(3) ^ 2) * MDP(14) + (t111 ^ 2 + t113 ^ 2) * MDP(18) + t93 * t138 + 0.2e1 * MDP(15) * t111 + 0.2e1 * MDP(16) * t113 + 0.2e1 * t162 + 0.2e1 * t161 + t149; (-t104 * t127 - t105 * t128) * MDP(17) + (t100 * t127 - t128 * t98) * MDP(18) + (MDP(14) * pkin(7) + MDP(12)) * t131 + t171 * (t103 * t92 - t106 * t91); -t128 * MDP(15) + t127 * MDP(16) + (-t111 * t128 + t113 * t127) * MDP(18) + t147 - t170; MDP(14) + (t127 ^ 2 + t128 ^ 2) * MDP(18); t139 * t91 + t157 + t158 + t159 + t164; 0; 0; MDP(18); (t146 * MDP(32) - t143) * t132 + (t146 * MDP(31) - t137) * t129 + t136; -t162 - t161 - t167 * t155 + (t167 * MDP(32) - 0.2e1 * t152) * t129 - t151; t170; 0; pkin(5) * t138 + t149; t82 * MDP(31) - t83 * MDP(32) + t142 * t92 + t163; (-MDP(32) * t94 - MDP(29)) * t132 + (-MDP(31) * t94 - MDP(28)) * t129; -t171 * t106; t141; MDP(28) * t129 + MDP(29) * t132 - pkin(9) * t171; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
