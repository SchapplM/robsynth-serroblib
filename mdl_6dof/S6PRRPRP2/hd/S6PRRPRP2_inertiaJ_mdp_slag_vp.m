% Calculate joint inertia matrix for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:57
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP2_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:54:51
% EndTime: 2021-01-16 02:54:55
% DurationCPUTime: 0.58s
% Computational Cost: add. (779->163), mult. (1560->240), div. (0->0), fcn. (1739->10), ass. (0->69)
t112 = sin(pkin(11));
t114 = cos(pkin(11));
t116 = sin(qJ(3));
t155 = cos(qJ(3));
t102 = t112 * t155 + t114 * t116;
t161 = 0.2e1 * t102;
t160 = qJ(4) + pkin(8);
t115 = sin(qJ(5));
t110 = t115 ^ 2;
t117 = cos(qJ(5));
t111 = t117 ^ 2;
t145 = t110 + t111;
t159 = MDP(24) * t145;
t107 = pkin(3) * t112 + pkin(9);
t132 = t107 * MDP(26) + MDP(24);
t138 = MDP(22) - MDP(25);
t139 = MDP(21) + MDP(23);
t158 = -t115 * t138 + t117 * t139 + MDP(12);
t113 = sin(pkin(6));
t154 = sin(qJ(2));
t137 = t113 * t154;
t151 = cos(pkin(6));
t120 = -t116 * t137 + t151 * t155;
t98 = t116 * t151 + t137 * t155;
t87 = t112 * t98 - t114 * t120;
t157 = t87 ^ 2;
t109 = -pkin(3) * t155 - pkin(2);
t156 = 0.2e1 * t109;
t101 = t112 * t116 - t114 * t155;
t153 = pkin(5) * t101;
t92 = t101 * pkin(4) - t102 * pkin(9) + t109;
t104 = t160 * t155;
t135 = t160 * t116;
t96 = t114 * t104 - t112 * t135;
t82 = t115 * t92 + t117 * t96;
t152 = MDP(15) * pkin(3);
t150 = qJ(6) * t101;
t149 = t101 * t107;
t148 = t102 * t117;
t118 = cos(qJ(2));
t147 = t113 * t118;
t108 = -pkin(3) * t114 - pkin(4);
t130 = -pkin(5) * t117 - qJ(6) * t115;
t99 = t108 + t130;
t146 = t99 * MDP(26);
t144 = MDP(17) * t117;
t143 = t101 * MDP(20);
t142 = t102 * MDP(13);
t140 = t109 * MDP(15);
t136 = MDP(10) * t155;
t134 = t115 * t96 - t117 * t92;
t133 = -MDP(26) * pkin(5) - MDP(23);
t94 = t104 * t112 + t114 * t135;
t131 = MDP(21) - t133;
t129 = pkin(5) * t115 - qJ(6) * t117;
t79 = t150 + t82;
t80 = t134 - t153;
t128 = t115 * t79 - t117 * t80;
t89 = t112 * t120 + t114 * t98;
t85 = t115 * t89 + t117 * t147;
t86 = -t115 * t147 + t117 * t89;
t127 = t115 * t86 - t117 * t85;
t125 = -t102 * t99 + t149;
t124 = MDP(26) * qJ(6) - t138;
t123 = -MDP(21) * t134 - t82 * MDP(22);
t122 = t102 * t108 - t149;
t121 = t117 * MDP(18) - t115 * MDP(19);
t83 = t102 * t129 + t94;
t1 = [MDP(1) + (t113 ^ 2 * t118 ^ 2 + t89 ^ 2 + t157) * MDP(15) + (t85 ^ 2 + t86 ^ 2 + t157) * MDP(26); (t87 * t94 + t89 * t96) * MDP(15) + (t79 * t86 + t80 * t85 + t83 * t87) * MDP(26) + (-t89 * MDP(14) - t138 * t86 - t139 * t85) * t101 + (-t154 * MDP(4) + (-MDP(11) * t116 - t101 * MDP(12) + MDP(3) + t136 - t140 - t142) * t118) * t113 + (-t127 * MDP(24) + (t115 * t139 + t117 * t138 + MDP(14)) * t87) * t102; MDP(2) + 0.2e1 * pkin(2) * t136 + t142 * t156 + (t109 ^ 2 + t94 ^ 2 + t96 ^ 2) * MDP(15) + (t79 ^ 2 + t80 ^ 2 + t83 ^ 2) * MDP(26) + (t111 * MDP(16) - 0.2e1 * t115 * t144) * t102 ^ 2 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t116 + 0.2e1 * MDP(6) * t155) * t116 + (MDP(12) * t156 + t121 * t161 + t143) * t101 + 0.2e1 * (-t96 * MDP(14) - t80 * MDP(23) + t79 * MDP(25) + t123) * t101 + (-t128 * MDP(24) + (t115 * MDP(23) - t117 * MDP(25)) * t83 + (t115 * MDP(21) + t117 * MDP(22) + MDP(14)) * t94) * t161; t120 * MDP(10) - t98 * MDP(11) + (t112 * t152 - MDP(13)) * t89 + (-t114 * t152 + t146 - t158) * t87 + t132 * (t115 * t85 + t117 * t86); t83 * t146 - t94 * MDP(12) - t96 * MDP(13) + t116 * MDP(7) + t155 * MDP(8) + (-t110 + t111) * MDP(17) * t102 + (-t116 * MDP(10) - MDP(11) * t155) * pkin(8) + (t101 * MDP(19) - t94 * MDP(21) + MDP(22) * t122 - t83 * MDP(23) + MDP(25) * t125 + t132 * t79) * t117 + (MDP(16) * t148 + t101 * MDP(18) + MDP(21) * t122 + t94 * MDP(22) - MDP(23) * t125 - t83 * MDP(25) + t132 * t80) * t115 + ((-t101 * t112 - t102 * t114) * MDP(14) + (t112 * t96 - t114 * t94) * MDP(15)) * pkin(3); MDP(9) + t110 * MDP(16) + (t107 ^ 2 * t145 + t99 ^ 2) * MDP(26) + (t112 ^ 2 + t114 ^ 2) * MDP(15) * pkin(3) ^ 2 + 0.2e1 * t107 * t159 + 0.2e1 * (-t108 * MDP(21) - t99 * MDP(23)) * t117 + 0.2e1 * (t114 * MDP(12) - t112 * MDP(13)) * pkin(3) + 0.2e1 * (MDP(22) * t108 - MDP(25) * t99 + t144) * t115; -MDP(15) * t147 + MDP(26) * t127; t140 + t128 * MDP(26) + (MDP(13) - t159) * t102 + t158 * t101; 0; MDP(26) * t145 + MDP(15); t124 * t86 - t131 * t85; t143 + (-t134 + 0.2e1 * t153) * MDP(23) + (0.2e1 * t150 + t82) * MDP(25) + (-pkin(5) * t80 + qJ(6) * t79) * MDP(26) + (MDP(24) * t130 + t121) * t102 + t123; t115 * MDP(18) + t117 * MDP(19) - t129 * MDP(24) + (-t115 * t131 + t117 * t124) * t107; t115 * t124 + t117 * t131; MDP(20) + 0.2e1 * pkin(5) * MDP(23) + 0.2e1 * qJ(6) * MDP(25) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(26); t85 * MDP(26); -t101 * MDP(23) + MDP(24) * t148 + t80 * MDP(26); t132 * t115; -t117 * MDP(26); t133; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
