% Calculate joint inertia matrix for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP5_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:49:17
% EndTime: 2019-03-08 21:49:19
% DurationCPUTime: 0.51s
% Computational Cost: add. (444->157), mult. (833->210), div. (0->0), fcn. (773->8), ass. (0->69)
t162 = pkin(4) + pkin(8);
t116 = sin(qJ(5));
t109 = t116 ^ 2;
t119 = cos(qJ(5));
t111 = t119 ^ 2;
t106 = t109 + t111;
t161 = t106 * MDP(26);
t122 = -pkin(3) - pkin(9);
t145 = MDP(26) * t122;
t160 = MDP(24) - t145;
t140 = MDP(21) + MDP(23);
t115 = cos(pkin(6));
t117 = sin(qJ(3));
t120 = cos(qJ(3));
t114 = sin(pkin(6));
t118 = sin(qJ(2));
t152 = t114 * t118;
t94 = t115 * t117 + t120 * t152;
t159 = t94 ^ 2;
t157 = 2 * MDP(23);
t156 = t117 * pkin(5);
t103 = t162 * t117;
t137 = -t117 * qJ(4) - pkin(2);
t96 = t122 * t120 + t137;
t84 = t116 * t103 + t119 * t96;
t155 = MDP(15) * pkin(8);
t154 = pkin(3) * MDP(15);
t93 = -t115 * t120 + t117 * t152;
t153 = t93 * t117;
t121 = cos(qJ(2));
t151 = t114 * t121;
t150 = t116 * t120;
t149 = t117 * qJ(6);
t148 = t117 * t122;
t147 = t119 * t120;
t104 = t162 * t120;
t101 = -t120 * pkin(3) + t137;
t146 = MDP(15) * t101;
t144 = qJ(4) * MDP(15);
t129 = t116 * pkin(5) - t119 * qJ(6);
t100 = qJ(4) + t129;
t143 = t100 * MDP(26);
t142 = pkin(8) ^ 2 * MDP(15);
t141 = -MDP(11) + MDP(14);
t139 = MDP(22) - MDP(25);
t138 = t119 * t116 * MDP(17);
t136 = -t119 * t103 + t116 * t96;
t134 = MDP(13) - t154;
t133 = -MDP(26) * pkin(5) - MDP(23);
t132 = t140 * t119;
t131 = -MDP(10) + t134;
t130 = MDP(21) - t133;
t102 = pkin(5) * t119 + t116 * qJ(6);
t85 = t116 * t151 + t93 * t119;
t87 = -t93 * t116 + t119 * t151;
t79 = -t87 * t116 + t85 * t119;
t128 = t141 + t144;
t127 = MDP(26) * qJ(6) - t139;
t126 = -t136 * MDP(21) - t84 * MDP(22);
t125 = -t116 * MDP(18) - t119 * MDP(19);
t124 = -t139 * t116 + t132;
t112 = t120 ^ 2;
t110 = t117 ^ 2;
t99 = t106 * t122;
t89 = t102 * t120 + t104;
t82 = t136 - t156;
t81 = t149 + t84;
t78 = t81 * t116 - t82 * t119;
t1 = [MDP(1) + (t114 ^ 2 * t121 ^ 2 + t93 ^ 2 + t159) * MDP(15) + (t85 ^ 2 + t87 ^ 2 + t159) * MDP(26); MDP(12) * t153 + (-t87 * t81 - t85 * t82 + t94 * t89) * MDP(26) - t139 * (-t87 * t117 + t94 * t150) + (t94 * MDP(12) + (t116 * t85 + t119 * t87) * MDP(24)) * t120 + (t94 * t120 + t153) * t155 + (-t118 * MDP(4) + (-t146 + MDP(3) + (MDP(10) - MDP(13)) * t120 + t141 * t117) * t121) * t114 + t140 * (t85 * t117 + t94 * t147); MDP(2) + (t81 ^ 2 + t82 ^ 2 + t89 ^ 2) * MDP(26) + t146 * t101 + (t109 * MDP(16) + 0.2e1 * t138 + t142) * t112 + (MDP(20) + MDP(5) + t142) * t110 + 0.2e1 * (t110 + t112) * MDP(12) * pkin(8) + 0.2e1 * (pkin(2) * MDP(10) + t101 * MDP(13) + (-t116 * t82 - t119 * t81) * MDP(24) + (t119 * MDP(23) + t116 * MDP(25)) * t89 + (t119 * MDP(21) - t116 * MDP(22)) * t104) * t120 + 0.2e1 * (-pkin(2) * MDP(11) - t101 * MDP(14) + (MDP(6) + t125) * t120 - t82 * MDP(23) + t81 * MDP(25) + t126) * t117; t131 * t93 + (t140 * t116 + t139 * t119 + t128 + t143) * t94 - t160 * t79; t89 * t143 - t78 * MDP(24) + (-pkin(3) * MDP(12) + MDP(7)) * t117 + t132 * t148 + (MDP(8) + qJ(4) * MDP(12) + (t109 - t111) * MDP(17)) * t120 + (-t82 * t145 + t117 * MDP(18) + t104 * MDP(22) - t89 * MDP(25) + (qJ(4) * MDP(21) + t100 * MDP(23)) * t120) * t119 + (-MDP(16) * t147 - t117 * MDP(19) + t104 * MDP(21) + (-qJ(4) * t120 - t148) * MDP(22) + t89 * MDP(23) + (t100 * t120 + t148) * MDP(25) + t81 * t145) * t116 + (t131 * t117 + t128 * t120) * pkin(8); -0.2e1 * t138 - 0.2e1 * t99 * MDP(24) + t111 * MDP(16) + MDP(9) + t122 ^ 2 * t161 + (-0.2e1 * t119 * MDP(25) + t116 * t157 + t143) * t100 + (-0.2e1 * MDP(13) + t154) * pkin(3) + (0.2e1 * t116 * MDP(21) + 0.2e1 * t119 * MDP(22) + 0.2e1 * MDP(14) + t144) * qJ(4); t93 * MDP(15) + t79 * MDP(26); t78 * MDP(26) + (MDP(12) + t124 + t155) * t117; -t106 * MDP(24) + t99 * MDP(26) + t134; MDP(15) + t161; -t127 * t87 + t130 * t85; t117 * MDP(20) + (-t136 + 0.2e1 * t156) * MDP(23) + (0.2e1 * t149 + t84) * MDP(25) + (-t82 * pkin(5) + t81 * qJ(6)) * MDP(26) + (t129 * MDP(24) + t125) * t120 + t126; t119 * MDP(18) - t116 * MDP(19) - t102 * MDP(24) + (t127 * t116 + t130 * t119) * t122; t102 * MDP(26) + t124; MDP(20) + pkin(5) * t157 + 0.2e1 * qJ(6) * MDP(25) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(26); -t85 * MDP(26); -t117 * MDP(23) - MDP(24) * t150 + t82 * MDP(26); t160 * t119; -t119 * MDP(26); t133; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
