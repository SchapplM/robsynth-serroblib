% Calculate joint inertia matrix for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRPRRR3_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:34:13
% EndTime: 2019-03-08 20:34:15
% DurationCPUTime: 0.42s
% Computational Cost: add. (627->114), mult. (1291->170), div. (0->0), fcn. (1546->12), ass. (0->72)
t128 = sin(qJ(6));
t132 = cos(qJ(6));
t142 = t132 * MDP(28) - t128 * MDP(29);
t172 = MDP(21) + t142;
t157 = t128 * MDP(25) + t132 * MDP(26);
t171 = MDP(8) * qJ(3);
t124 = sin(pkin(12));
t126 = cos(pkin(12));
t130 = sin(qJ(4));
t134 = cos(qJ(4));
t109 = t130 * t124 - t134 * t126;
t114 = -t126 * pkin(3) - pkin(2);
t104 = t109 * pkin(4) + t114;
t170 = 0.2e1 * t104;
t169 = 0.2e1 * t114;
t168 = pkin(2) * MDP(8);
t166 = pkin(8) + qJ(3);
t129 = sin(qJ(5));
t133 = cos(qJ(5));
t110 = t134 * t124 + t130 * t126;
t111 = t166 * t124;
t112 = t166 * t126;
t148 = -t134 * t111 - t130 * t112;
t92 = -t110 * pkin(9) + t148;
t144 = t130 * t111 - t134 * t112;
t93 = -t109 * pkin(9) - t144;
t84 = t129 * t93 - t133 * t92;
t165 = t84 * t132;
t103 = -t129 * t109 + t133 * t110;
t164 = t103 * t128;
t163 = t103 * t132;
t162 = t124 * MDP(6);
t125 = sin(pkin(6));
t131 = sin(qJ(2));
t161 = t125 * t131;
t135 = cos(qJ(2));
t160 = t125 * t135;
t159 = t126 * MDP(5);
t158 = t128 * t132;
t102 = t133 * t109 + t129 * t110;
t155 = t102 * MDP(27);
t154 = t103 * MDP(22);
t153 = t109 * MDP(14);
t150 = MDP(24) * t158;
t122 = t128 ^ 2;
t149 = t122 * MDP(23) + MDP(20) + 0.2e1 * t150;
t147 = -pkin(5) * t103 - pkin(10) * t102;
t115 = t129 * pkin(4) + pkin(10);
t116 = -t133 * pkin(4) - pkin(5);
t146 = -t102 * t115 + t103 * t116;
t143 = MDP(25) * t132 - MDP(26) * t128;
t141 = -MDP(28) * t128 - MDP(29) * t132;
t127 = cos(pkin(6));
t105 = -t124 * t161 + t127 * t126;
t106 = t127 * t124 + t126 * t161;
t95 = t134 * t105 - t130 * t106;
t96 = t130 * t105 + t134 * t106;
t89 = t129 * t96 - t133 * t95;
t90 = t129 * t95 + t133 * t96;
t140 = -t90 * MDP(22) - t172 * t89;
t139 = (t133 * MDP(21) - t129 * MDP(22)) * pkin(4);
t123 = t132 ^ 2;
t85 = t129 * t92 + t133 * t93;
t138 = -t84 * MDP(21) - t85 * MDP(22) + (MDP(18) + (-t122 + t123) * MDP(24) + MDP(23) * t158) * t103 + (-MDP(19) + t157) * t102;
t137 = t110 * MDP(15) + t153 + t154 - t159 + t162 - t168;
t86 = t102 * pkin(5) - t103 * pkin(10) + t104;
t80 = t84 * t128;
t79 = -t128 * t160 + t132 * t90;
t78 = -t128 * t90 - t132 * t160;
t77 = t128 * t86 + t132 * t85;
t76 = -t128 * t85 + t132 * t86;
t1 = [MDP(1) + (t125 ^ 2 * t135 ^ 2 + t105 ^ 2 + t106 ^ 2) * MDP(8); (t78 * t102 + t89 * t164) * MDP(28) + (-t79 * t102 + t89 * t163) * MDP(29) + (-t131 * MDP(4) + (-t102 * MDP(21) + MDP(3) - t137) * t135) * t125 + (MDP(7) + t171) * (-t105 * t124 + t106 * t126); t153 * t169 + t154 * t170 + MDP(2) + (0.2e1 * t159 - 0.2e1 * t162 + t168) * pkin(2) + (-0.2e1 * t109 * MDP(10) + MDP(15) * t169 + MDP(9) * t110) * t110 + (t123 * MDP(23) + MDP(16) - 0.2e1 * t150) * t103 ^ 2 + (MDP(21) * t170 + t155 + 0.2e1 * (-MDP(17) + t143) * t103) * t102 + 0.2e1 * (t76 * t102 + t84 * t164) * MDP(28) + 0.2e1 * (-t77 * t102 + t84 * t163) * MDP(29) + (0.2e1 * MDP(7) + t171) * (t124 ^ 2 + t126 ^ 2) * qJ(3); -MDP(8) * t160; t172 * t102 + t137; MDP(8); t95 * MDP(14) - t96 * MDP(15) + t140; t110 * MDP(11) - t109 * MDP(12) + t148 * MDP(14) + t144 * MDP(15) + (t146 * t128 - t165) * MDP(28) + (t146 * t132 + t80) * MDP(29) + t138; 0; -0.2e1 * t116 * t142 + MDP(13) + 0.2e1 * t139 + t149; t140; (t147 * t128 - t165) * MDP(28) + (t147 * t132 + t80) * MDP(29) + t138; 0; t139 + t149 + t142 * (pkin(5) - t116); 0.2e1 * pkin(5) * t142 + t149; t78 * MDP(28) - t79 * MDP(29); t76 * MDP(28) - t77 * MDP(29) + t143 * t103 + t155; t142; t141 * t115 + t157; t141 * pkin(10) + t157; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
