% Calculate joint inertia matrix for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPPR1_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:39:36
% EndTime: 2019-03-09 02:39:37
% DurationCPUTime: 0.40s
% Computational Cost: add. (735->123), mult. (1316->188), div. (0->0), fcn. (1431->10), ass. (0->69)
t115 = sin(pkin(10));
t106 = t115 * pkin(3) + qJ(5);
t114 = sin(pkin(11));
t117 = cos(pkin(11));
t141 = t114 ^ 2 + t117 ^ 2;
t134 = t141 * MDP(17);
t156 = t106 * t134;
t116 = sin(pkin(9));
t109 = t116 * pkin(1) + pkin(7);
t155 = qJ(4) + t109;
t154 = t141 * MDP(16);
t118 = cos(pkin(10));
t121 = sin(qJ(3));
t135 = t155 * t121;
t151 = cos(qJ(3));
t95 = t155 * t151;
t88 = t115 * t95 + t118 * t135;
t153 = t88 ^ 2;
t98 = t115 * t121 - t118 * t151;
t96 = t98 ^ 2;
t110 = -t118 * pkin(3) - pkin(4);
t103 = -t117 * pkin(5) + t110;
t152 = 0.2e1 * t103;
t150 = t88 * t98;
t149 = pkin(8) + t106;
t101 = t115 * t151 + t118 * t121;
t119 = cos(pkin(9));
t111 = -t119 * pkin(1) - pkin(2);
t104 = -t151 * pkin(3) + t111;
t83 = t98 * pkin(4) - t101 * qJ(5) + t104;
t90 = -t115 * t135 + t118 * t95;
t78 = t114 * t83 + t117 * t90;
t148 = MDP(13) * pkin(3);
t147 = t101 * t114;
t120 = sin(qJ(6));
t122 = cos(qJ(6));
t102 = t122 * t114 + t120 * t117;
t84 = t102 * t101;
t146 = t84 * MDP(19);
t100 = t120 * t114 - t122 * t117;
t85 = t100 * t101;
t145 = t85 * MDP(18);
t144 = t101 * t154;
t143 = t98 * MDP(20);
t142 = t98 * MDP(21);
t140 = t100 * MDP(23);
t139 = t110 * MDP(17);
t138 = t114 * MDP(15);
t137 = t117 * MDP(14);
t136 = t151 * MDP(10);
t77 = -t114 * t90 + t117 * t83;
t133 = t78 * t114 + t77 * t117;
t132 = -t77 * t114 + t78 * t117;
t131 = t101 * t110 - t106 * t98;
t75 = -t117 * t101 * pkin(8) + t98 * pkin(5) + t77;
t76 = -pkin(8) * t147 + t78;
t130 = (-t120 * t76 + t122 * t75) * MDP(23) - (t120 * t75 + t122 * t76) * MDP(24);
t129 = t84 * MDP(23) - t85 * MDP(24);
t128 = t114 * MDP(14) + t117 * MDP(15);
t127 = -t102 * MDP(24) - t140;
t126 = t127 + t137 - t138;
t125 = -t126 + t139;
t97 = t101 ^ 2;
t94 = t149 * t117;
t93 = t149 * t114;
t87 = -t120 * t93 + t122 * t94;
t86 = -t120 * t94 - t122 * t93;
t79 = pkin(5) * t147 + t88;
t1 = [MDP(1) - 0.2e1 * t111 * t136 + (t104 ^ 2 + t90 ^ 2 + t153) * MDP(13) + (t77 ^ 2 + t78 ^ 2 + t153) * MDP(17) - 0.2e1 * t84 * t142 + t96 * MDP(22) + (t116 ^ 2 + t119 ^ 2) * MDP(4) * pkin(1) ^ 2 - (0.2e1 * t143 - t145 - 0.2e1 * t146) * t85 + (0.2e1 * t111 * MDP(11) + MDP(5) * t121 + 0.2e1 * MDP(6) * t151) * t121 + 0.2e1 * t129 * t79 + 0.2e1 * (-t90 * MDP(12) + t77 * MDP(14) - t78 * MDP(15) + t130) * t98 + 0.2e1 * (-t133 * MDP(16) + (MDP(12) + t128) * t88) * t101; (t90 * t101 + t150) * MDP(13) + (t101 * t132 + t150) * MDP(17); MDP(4) + (t97 + t96) * MDP(13) + (t141 * t97 + t96) * MDP(17); t121 * MDP(7) + t151 * MDP(8) + (t114 * t131 - t88 * t117) * MDP(14) + (t88 * t114 + t117 * t131) * MDP(15) + t132 * MDP(16) + (t106 * t132 + t88 * t110) * MDP(17) + (t103 * t84 + t86 * t98) * MDP(23) + (-t103 * t85 - t87 * t98) * MDP(24) + (-t121 * MDP(10) - MDP(11) * t151) * t109 + (t85 * MDP(19) + t79 * MDP(23) - t142) * t100 + (t79 * MDP(24) + t143 - t145 - t146) * t102 + ((-t101 * t118 - t115 * t98) * MDP(12) + (t115 * t90 - t118 * t88) * MDP(13)) * pkin(3); t136 - t121 * MDP(11) + t144 + (t115 * t148 + t156) * t101 + (-t118 * t148 + t125) * t98; t140 * t152 + MDP(9) + (t115 ^ 2 + t118 ^ 2) * MDP(13) * pkin(3) ^ 2 + (-0.2e1 * t137 + 0.2e1 * t138 + t139) * t110 + (MDP(18) * t102 - 0.2e1 * t100 * MDP(19) + MDP(24) * t152) * t102 + (0.2e1 * t154 + t156) * t106; t104 * MDP(13) + MDP(17) * t133 + t126 * t98 - t144; 0; 0; MDP(13) + t134; t88 * MDP(17) + t101 * t128 + t129; t98 * MDP(17); t125; 0; MDP(17); -t85 * MDP(20) - t84 * MDP(21) + t98 * MDP(22) + t130; -t129; t102 * MDP(20) - t100 * MDP(21) + t86 * MDP(23) - t87 * MDP(24); t127; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
