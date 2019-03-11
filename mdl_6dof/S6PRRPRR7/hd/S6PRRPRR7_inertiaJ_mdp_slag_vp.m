% Calculate joint inertia matrix for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR7_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:34:12
% EndTime: 2019-03-08 22:34:14
% DurationCPUTime: 0.46s
% Computational Cost: add. (409->143), mult. (822->199), div. (0->0), fcn. (828->10), ass. (0->80)
t183 = pkin(4) + pkin(8);
t132 = sin(qJ(6));
t136 = cos(qJ(6));
t137 = cos(qJ(5));
t138 = cos(qJ(3));
t166 = t137 * t138;
t133 = sin(qJ(5));
t167 = t133 * t138;
t102 = t132 * t167 - t136 * t166;
t113 = t132 * t137 + t136 * t133;
t103 = t113 * t138;
t182 = -t103 * MDP(25) + t102 * MDP(26);
t140 = -pkin(3) - pkin(9);
t114 = -t132 * t133 + t136 * t137;
t176 = -pkin(10) + t140;
t116 = t176 * t133;
t117 = t176 * t137;
t150 = t114 * MDP(25) - t113 * MDP(26) + (-t132 * t116 + t136 * t117) * MDP(28) - (t136 * t116 + t132 * t117) * MDP(29);
t181 = (MDP(21) * t140 + MDP(18)) * t137 + (-MDP(22) * t140 - MDP(19)) * t133 + t150;
t180 = pkin(8) * MDP(15) + MDP(12);
t179 = 0.2e1 * t138;
t178 = 0.2e1 * MDP(29);
t134 = sin(qJ(3));
t177 = t134 * pkin(5);
t131 = cos(pkin(6));
t130 = sin(pkin(6));
t135 = sin(qJ(2));
t171 = t130 * t135;
t105 = -t131 * t138 + t134 * t171;
t139 = cos(qJ(2));
t170 = t130 * t139;
t95 = t105 * t137 + t133 * t170;
t96 = -t105 * t133 + t137 * t170;
t86 = t132 * t96 + t136 * t95;
t87 = t132 * t95 - t136 * t96;
t175 = t86 * MDP(28) - t87 * MDP(29);
t174 = MDP(15) * pkin(3);
t153 = -t134 * qJ(4) - pkin(2);
t112 = t140 * t138 + t153;
t152 = pkin(10) * t138 - t112;
t119 = t183 * t134;
t169 = t133 * t119;
t90 = -t152 * t137 + t169;
t172 = t136 * t90;
t168 = t133 * t137;
t165 = t114 * MDP(28) - t113 * MDP(29);
t120 = t183 * t138;
t118 = -t138 * pkin(3) + t153;
t164 = MDP(15) * t118;
t163 = qJ(4) * MDP(15);
t162 = t103 * MDP(23);
t161 = t113 * MDP(28);
t160 = t133 * MDP(21);
t159 = t137 * MDP(22);
t158 = pkin(8) ^ 2 * MDP(15);
t157 = -MDP(11) + MDP(14);
t156 = MDP(20) + MDP(27);
t155 = t134 * MDP(27) + t182;
t154 = MDP(17) * t168;
t115 = t137 * t119;
t89 = t152 * t133 + t115 + t177;
t82 = -t132 * t90 + t136 * t89;
t151 = MDP(13) - t174;
t149 = -MDP(10) + t151;
t147 = -t133 * MDP(18) - t137 * MDP(19);
t146 = t137 * MDP(21) - t133 * MDP(22);
t145 = t159 + t160;
t144 = (MDP(28) * t136 - MDP(29) * t132) * pkin(5);
t143 = t146 + t180;
t129 = t138 ^ 2;
t128 = t137 ^ 2;
t127 = t134 ^ 2;
t126 = t133 ^ 2;
t122 = t133 * pkin(5) + qJ(4);
t106 = t131 * t134 + t138 * t171;
t104 = pkin(5) * t166 + t120;
t92 = t137 * t112 + t169;
t91 = -t133 * t112 + t115;
t83 = t132 * t89 + t172;
t1 = [MDP(1) + (t130 ^ 2 * t139 ^ 2 + t105 ^ 2 + t106 ^ 2) * MDP(15); (t106 * t166 + t95 * t134) * MDP(21) + (-t106 * t167 + t96 * t134) * MDP(22) + (-t106 * t102 + t86 * t134) * MDP(28) + (-t106 * t103 - t87 * t134) * MDP(29) + (-t135 * MDP(4) + (-t164 + MDP(3) + (MDP(10) - MDP(13)) * t138 + t157 * t134) * t139) * t130 + t180 * (t105 * t134 + t106 * t138); pkin(2) * MDP(10) * t179 + MDP(2) + (MDP(13) * t179 + t164) * t118 + (t126 * MDP(16) + 0.2e1 * t154 + t158) * t129 - (0.2e1 * t102 * MDP(24) - t162) * t103 + (MDP(5) + t156 + t158) * t127 + 0.2e1 * (-pkin(2) * MDP(11) - t118 * MDP(14) + (MDP(6) + t147) * t138 + t182) * t134 + 0.2e1 * (t120 * t166 + t91 * t134) * MDP(21) + 0.2e1 * (-t120 * t167 - t92 * t134) * MDP(22) + 0.2e1 * (-t104 * t102 + t82 * t134) * MDP(28) + (-t104 * t103 - t83 * t134) * t178 + 0.2e1 * (t127 + t129) * MDP(12) * pkin(8); t149 * t105 + (t114 * MDP(29) + t145 + t157 + t161 + t163) * t106; -t114 * t162 + (t114 * t102 + t103 * t113) * MDP(24) + (-t122 * t102 + t104 * t113) * MDP(28) + (-t122 * t103 + t104 * t114) * MDP(29) + t145 * t120 + (MDP(8) - MDP(16) * t168 + (t126 - t128) * MDP(17) + t157 * pkin(8) + t143 * qJ(4)) * t138 + (-pkin(3) * MDP(12) + t149 * pkin(8) + MDP(7) + t181) * t134; -0.2e1 * t154 + 0.2e1 * t122 * t161 + t128 * MDP(16) + MDP(9) + (-0.2e1 * MDP(13) + t174) * pkin(3) + (MDP(23) * t114 - 0.2e1 * t113 * MDP(24) + t122 * t178) * t114 + (0.2e1 * MDP(14) + 0.2e1 * t159 + 0.2e1 * t160 + t163) * qJ(4); t105 * MDP(15); (t143 + t165) * t134; t151; MDP(15); t95 * MDP(21) + t96 * MDP(22) + t175; t134 * MDP(20) + t91 * MDP(21) - t92 * MDP(22) + (t136 * t177 + t82) * MDP(28) + (-t172 + (-t89 - t177) * t132) * MDP(29) + t147 * t138 + t155; t181; t146 + t165; 0.2e1 * t144 + t156; t175; t82 * MDP(28) - t83 * MDP(29) + t155; t150; t165; MDP(27) + t144; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
