% Calculate joint inertia matrix for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR3_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:14:12
% EndTime: 2019-03-08 23:14:13
% DurationCPUTime: 0.43s
% Computational Cost: add. (500->142), mult. (993->202), div. (0->0), fcn. (1062->10), ass. (0->69)
t134 = cos(qJ(6));
t152 = t134 * MDP(29);
t130 = sin(qJ(6));
t154 = t130 * MDP(28);
t171 = t152 + t154;
t151 = MDP(17) - MDP(20);
t150 = -MDP(18) + MDP(21);
t170 = -2 * MDP(20);
t169 = 2 * MDP(21);
t168 = pkin(4) + pkin(10);
t167 = -pkin(9) - pkin(8);
t166 = (MDP(22) * pkin(4));
t131 = sin(qJ(4));
t132 = sin(qJ(3));
t135 = cos(qJ(4));
t136 = cos(qJ(3));
t113 = t131 * t132 - t135 * t136;
t165 = qJ(5) * t113;
t164 = t113 * t130;
t163 = t113 * t134;
t119 = t131 * pkin(3) + qJ(5);
t162 = t119 * t113;
t128 = sin(pkin(6));
t133 = sin(qJ(2));
t161 = t128 * t133;
t137 = cos(qJ(2));
t160 = t128 * t137;
t159 = t130 * t134;
t157 = MDP(10) * t136;
t114 = t131 * t136 + t135 * t132;
t156 = t114 * MDP(26);
t122 = -t135 * pkin(3) - pkin(4);
t155 = t122 * MDP(22);
t125 = t134 * MDP(25);
t153 = t134 * MDP(28);
t149 = 0.2e1 * t114;
t123 = -t136 * pkin(3) - pkin(2);
t148 = MDP(24) * t159;
t127 = t134 ^ 2;
t147 = t127 * MDP(23) + MDP(16) - 0.2e1 * t148;
t129 = cos(pkin(6));
t107 = t129 * t136 - t132 * t161;
t108 = t129 * t132 + t136 * t161;
t93 = -t135 * t107 + t131 * t108;
t94 = t131 * t107 + t135 * t108;
t146 = -t151 * t93 + (t150 + t171) * t94;
t145 = t114 * t168 + t165;
t118 = -pkin(10) + t122;
t144 = -t114 * t118 + t162;
t115 = t167 * t132;
t116 = t167 * t136;
t102 = -t135 * t115 - t131 * t116;
t103 = t131 * t115 - t135 * t116;
t143 = MDP(25) * t130 + MDP(26) * t134;
t142 = -t130 * MDP(29) + t153;
t141 = -t114 * qJ(5) + t123;
t140 = t169 + 0.2e1 * t152 + 0.2e1 * t154;
t126 = t130 ^ 2;
t92 = -t113 * pkin(5) + t103;
t139 = t92 * t152 + (t125 + MDP(14)) * t114 + t150 * t103 - t151 * t102 + (MDP(23) * t159 - MDP(15) + (-t126 + t127) * MDP(24)) * t113;
t97 = t113 * pkin(4) + t141;
t91 = t114 * pkin(5) + t102;
t85 = t92 * t130;
t82 = t168 * t113 + t141;
t81 = -t130 * t93 + t134 * t160;
t80 = t130 * t160 + t134 * t93;
t79 = t130 * t91 + t134 * t82;
t78 = -t130 * t82 + t134 * t91;
t1 = [MDP(1) + (t128 ^ 2 * t137 ^ 2 + t93 ^ 2 + t94 ^ 2) * MDP(22); (-t94 * t113 + t93 * t114) * MDP(19) + (t93 * t102 + t94 * t103) * MDP(22) + (t80 * t114 - t94 * t163) * MDP(28) + (t81 * t114 + t94 * t164) * MDP(29) + (-t133 * MDP(4) + (-MDP(11) * t132 - MDP(22) * t97 - t151 * t113 + t150 * t114 + MDP(3) + t157) * t137) * t128; MDP(2) + 0.2e1 * pkin(2) * t157 + (t102 ^ 2 + t103 ^ 2 + t97 ^ 2) * MDP(22) + (t123 * MDP(18) - t97 * MDP(21)) * t149 + (MDP(12) + MDP(27)) * t114 ^ 2 + (t126 * MDP(23) + 0.2e1 * t148) * t113 ^ 2 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t132 + 0.2e1 * t136 * MDP(6)) * t132 + (0.2e1 * t123 * MDP(17) + t97 * t170 + (-MDP(13) + t143) * t149) * t113 + 0.2e1 * (t102 * t114 - t103 * t113) * MDP(19) + 0.2e1 * (t78 * t114 - t92 * t163) * MDP(28) + 0.2e1 * (-t79 * t114 + t92 * t164) * MDP(29); t107 * MDP(10) - t108 * MDP(11) + (t94 * t119 + t93 * t122) * MDP(22) + t146; (-t144 * t134 + t85) * MDP(28) + t136 * MDP(8) + t132 * MDP(7) + (-t132 * MDP(10) - t136 * MDP(11)) * pkin(8) + (t144 * MDP(29) - t156) * t130 + t139 + (t122 * t114 - t162) * MDP(19) + (t102 * t122 + t103 * t119) * MDP(22); MDP(9) + ((2 * MDP(20)) + t155) * t122 + 0.2e1 * (t135 * MDP(17) - t131 * MDP(18)) * pkin(3) + (t119 * MDP(22) + t140) * t119 + t147; (-t93 * pkin(4) + t94 * qJ(5)) * MDP(22) + t146; (t145 * MDP(29) - t156) * t130 + (-t145 * t134 + t85) * MDP(28) + t139 + (-t102 * pkin(4) + t103 * qJ(5)) * MDP(22) + (-pkin(4) * t114 - t165) * MDP(19); pkin(4) * t170 + qJ(5) * t169 + (-t122 * pkin(4) + t119 * qJ(5)) * MDP(22) + (t150 * t131 + t151 * t135) * pkin(3) + t147 + t171 * (qJ(5) + t119); (t170 + t166) * pkin(4) + (MDP(22) * qJ(5) + t140) * qJ(5) + t147; t93 * MDP(22); t102 * MDP(22) + (MDP(19) + t142) * t114; MDP(20) + t155; MDP(20) - t166; MDP(22); t80 * MDP(28) + t81 * MDP(29); t114 * MDP(27) + t78 * MDP(28) - t79 * MDP(29) + t143 * t113; -t130 * MDP(26) + t142 * t118 + t125; -t168 * t153 + t125 + (MDP(29) * t168 - MDP(26)) * t130; t142; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
