% Calculate joint inertia matrix for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRRR2_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:29:57
% EndTime: 2019-03-08 20:29:58
% DurationCPUTime: 0.43s
% Computational Cost: add. (404->122), mult. (901->190), div. (0->0), fcn. (971->12), ass. (0->67)
t124 = sin(qJ(5));
t128 = cos(qJ(5));
t135 = MDP(18) * t124 + MDP(19) * t128;
t123 = sin(qJ(6));
t127 = cos(qJ(6));
t106 = t123 * t124 - t127 * t128;
t107 = t123 * t128 + t127 * t124;
t156 = pkin(9) + pkin(10);
t109 = t156 * t124;
t110 = t156 * t128;
t139 = t107 * MDP(22) - t106 * MDP(23) + (-t127 * t109 - t123 * t110) * MDP(25) - (-t123 * t109 + t127 * t110) * MDP(26);
t161 = t124 * MDP(15) + t128 * MDP(16) - pkin(9) * t135 + t139;
t133 = (MDP(25) * t127 - MDP(26) * t123) * pkin(5);
t125 = sin(qJ(4));
t100 = t106 * t125;
t99 = t107 * t125;
t153 = -t100 * MDP(22) - t99 * MDP(23);
t121 = cos(pkin(12));
t112 = -t121 * pkin(2) - pkin(3);
t129 = cos(qJ(4));
t105 = -t129 * pkin(4) - t125 * pkin(9) + t112;
t101 = t128 * t105;
t148 = t125 * t128;
t119 = sin(pkin(12));
t111 = t119 * pkin(2) + pkin(8);
t152 = t111 * t124;
t82 = -pkin(10) * t148 + t101 + (-pkin(5) - t152) * t129;
t150 = t111 * t129;
t141 = t128 * t150;
t83 = t141 + (-pkin(10) * t125 + t105) * t124;
t77 = -t123 * t83 + t127 * t82;
t78 = t123 * t82 + t127 * t83;
t160 = t77 * MDP(25) - t78 * MDP(26) + t153;
t158 = -2 * MDP(21);
t157 = 0.2e1 * MDP(26);
t122 = cos(pkin(6));
t120 = sin(pkin(6));
t126 = sin(qJ(2));
t130 = cos(qJ(2));
t98 = (t119 * t130 + t121 * t126) * t120;
t87 = t122 * t125 + t98 * t129;
t96 = (t119 * t126 - t121 * t130) * t120;
t79 = -t87 * t124 + t96 * t128;
t80 = t96 * t124 + t87 * t128;
t75 = -t123 * t80 + t127 * t79;
t76 = t123 * t79 + t127 * t80;
t155 = t75 * MDP(25) - t76 * MDP(26);
t154 = -t99 * MDP(25) + t100 * MDP(26);
t151 = t111 * t128;
t149 = t124 * t128;
t145 = t106 * MDP(25);
t144 = t107 * MDP(20);
t143 = t125 * MDP(12);
t142 = MDP(17) + MDP(24);
t140 = MDP(14) * t149;
t138 = t128 * MDP(15) - t124 * MDP(16);
t136 = t128 * MDP(18) - t124 * MDP(19);
t132 = -t107 * MDP(26) + MDP(11) + t136 - t145;
t117 = t128 ^ 2;
t116 = t125 ^ 2;
t115 = t124 ^ 2;
t114 = -t128 * pkin(5) - pkin(4);
t102 = (pkin(5) * t124 + t111) * t125;
t86 = -t122 * t129 + t98 * t125;
t85 = t124 * t105 + t141;
t84 = -t124 * t150 + t101;
t1 = [MDP(1) + (t122 ^ 2 + t96 ^ 2 + t98 ^ 2) * MDP(5); (t86 * t124 * t125 - t79 * t129) * MDP(18) + (t80 * t129 + t148 * t86) * MDP(19) + (-t75 * t129 + t86 * t99) * MDP(25) + (-t86 * t100 + t76 * t129) * MDP(26) + (-t129 * MDP(11) + t143) * t96 + (t130 * MDP(3) - t126 * MDP(4)) * t120 + (t119 * t98 - t121 * t96) * MDP(5) * pkin(2); 0.2e1 * t112 * t143 + MDP(2) + (t119 ^ 2 + t121 ^ 2) * MDP(5) * pkin(2) ^ 2 + t142 * t129 ^ 2 - (-t100 * MDP(20) + t158 * t99) * t100 + (t117 * MDP(13) + MDP(6) - 0.2e1 * t140) * t116 + 0.2e1 * (-t112 * MDP(11) + (MDP(7) - t138) * t125 - t153) * t129 + 0.2e1 * (t116 * t152 - t84 * t129) * MDP(18) + 0.2e1 * (t116 * t151 + t85 * t129) * MDP(19) + 0.2e1 * (t102 * t99 - t77 * t129) * MDP(25) + (-t102 * t100 + t78 * t129) * t157; t122 * MDP(5); 0; MDP(5); -t87 * MDP(12) - t132 * t86; -t100 * t144 + (t100 * t106 - t107 * t99) * MDP(21) + (t102 * t106 + t114 * t99) * MDP(25) + (-t114 * t100 + t102 * t107) * MDP(26) + (-t111 * MDP(12) + MDP(9) - t161) * t129 + (MDP(8) - t111 * MDP(11) + MDP(13) * t149 + (-t115 + t117) * MDP(14) + (-pkin(4) * t124 - t151) * MDP(18) + (-pkin(4) * t128 + t152) * MDP(19)) * t125; t129 * t132 - t143; 0.2e1 * t140 + 0.2e1 * t114 * t145 + t115 * MDP(13) + MDP(10) + 0.2e1 * t136 * pkin(4) + (t106 * t158 + t114 * t157 + t144) * t107; t79 * MDP(18) - t80 * MDP(19) + t155; t84 * MDP(18) - t85 * MDP(19) + (-t142 - t133) * t129 + t138 * t125 + t160; -t125 * t135 + t154; t161; 0.2e1 * t133 + t142; t155; -t129 * MDP(24) + t160; t154; t139; MDP(24) + t133; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
