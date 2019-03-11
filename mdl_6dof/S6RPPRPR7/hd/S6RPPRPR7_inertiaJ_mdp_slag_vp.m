% Calculate joint inertia matrix for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRPR7_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:53:45
% EndTime: 2019-03-09 01:53:46
% DurationCPUTime: 0.40s
% Computational Cost: add. (657->121), mult. (1138->180), div. (0->0), fcn. (1248->8), ass. (0->66)
t118 = sin(qJ(6));
t120 = cos(qJ(6));
t115 = cos(pkin(10));
t114 = sin(pkin(9));
t116 = cos(pkin(9));
t119 = sin(qJ(4));
t149 = cos(qJ(4));
t94 = t114 * t119 - t116 * t149;
t144 = t115 * t94;
t113 = sin(pkin(10));
t103 = pkin(3) * t114 + qJ(2);
t95 = t114 * t149 + t116 * t119;
t85 = pkin(4) * t95 + qJ(5) * t94 + t103;
t117 = -pkin(1) - qJ(3);
t148 = -pkin(7) + t117;
t98 = t148 * t114;
t99 = t148 * t116;
t88 = t119 * t99 + t149 * t98;
t76 = -t113 * t88 + t115 * t85;
t74 = pkin(5) * t95 + pkin(8) * t144 + t76;
t145 = t113 * t94;
t77 = t113 * t85 + t115 * t88;
t75 = pkin(8) * t145 + t77;
t96 = t113 * t120 + t115 * t118;
t82 = t96 * t94;
t93 = t113 * t118 - t115 * t120;
t84 = t93 * t94;
t156 = (-t118 * t75 + t120 * t74) * MDP(27) - (t118 * t74 + t120 * t75) * MDP(28) + t84 * MDP(24) + t82 * MDP(25);
t152 = t94 ^ 2;
t91 = t95 ^ 2;
t154 = -t91 - t152;
t139 = t113 ^ 2 + t115 ^ 2;
t135 = t139 * MDP(21);
t153 = qJ(5) * t135;
t106 = -pkin(5) * t115 - pkin(4);
t151 = 0.2e1 * t106;
t150 = -2 * MDP(23);
t147 = pkin(8) + qJ(5);
t146 = pkin(4) * MDP(21);
t143 = MDP(21) * t95;
t142 = t84 * MDP(22);
t141 = t93 * MDP(27);
t140 = t94 * MDP(21);
t102 = t114 ^ 2 + t116 ^ 2;
t138 = t113 * MDP(19);
t137 = t115 * MDP(18);
t136 = t139 * MDP(20);
t134 = t113 * t77 + t115 * t76;
t133 = -t113 * t76 + t115 * t77;
t130 = -t82 * MDP(27) + t84 * MDP(28);
t129 = -MDP(28) * t96 - t141;
t128 = t114 * MDP(7) + t116 * MDP(8);
t127 = -MDP(17) + t136;
t126 = -MDP(18) * t113 - MDP(19) * t115;
t87 = t119 * t98 - t149 * t99;
t100 = t147 * t113;
t101 = t147 * t115;
t125 = t96 * MDP(24) - t93 * MDP(25) + (-t100 * t120 - t101 * t118) * MDP(27) - (-t100 * t118 + t101 * t120) * MDP(28);
t124 = t129 + t137 - t138;
t123 = -t124 - t146;
t122 = qJ(2) ^ 2;
t92 = t102 * t117;
t83 = t93 * t95;
t81 = t96 * t95;
t78 = -pkin(5) * t145 + t87;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2) + t122) * MDP(6) + (t102 * t117 ^ 2 + t122) * MDP(10) + (t76 ^ 2 + t77 ^ 2 + t87 ^ 2) * MDP(21) + t91 * MDP(26) + (MDP(11) * t94 - 0.2e1 * t103 * MDP(17)) * t94 + (-t150 * t82 + t142) * t84 - 0.2e1 * t92 * MDP(9) + 0.2e1 * t130 * t78 + 0.2e1 * (MDP(20) * t134 + t126 * t87) * t94 + 0.2e1 * (MDP(5) + t128) * qJ(2) + 0.2e1 * (MDP(12) * t94 + MDP(16) * t103 + t76 * MDP(18) - t77 * MDP(19) + t156) * t95; MDP(4) - pkin(1) * MDP(6) - t102 * MDP(9) + t92 * MDP(10) + t87 * t140 + (-t81 * t95 - t82 * t94) * MDP(27) + (t83 * t95 + t84 * t94) * MDP(28) + (MDP(19) * t154 + t77 * t143) * t115 + (MDP(18) * t154 - t76 * t143) * t113; MDP(6) + t102 * MDP(10) + (t139 * t91 + t152) * MDP(21); qJ(2) * MDP(10) + t134 * MDP(21) + t127 * t94 + (MDP(16) + t124) * t95 + t128; 0; MDP(10) + t135; -t94 * MDP(13) - t87 * MDP(16) - t88 * MDP(17) + (pkin(4) * t145 - t87 * t115) * MDP(18) + (pkin(4) * t144 + t87 * t113) * MDP(19) + t133 * MDP(20) + (-t87 * pkin(4) + qJ(5) * t133) * MDP(21) + t96 * t142 + (t82 * t96 - t84 * t93) * MDP(23) + (-t106 * t82 + t78 * t93) * MDP(27) + (t106 * t84 + t78 * t96) * MDP(28) + (qJ(5) * t126 - MDP(14) + t125) * t95; (t127 + t153) * t95 + (-MDP(16) + t123) * t94; 0; t141 * t151 + MDP(15) + (0.2e1 * t137 - 0.2e1 * t138 + t146) * pkin(4) + (MDP(22) * t96 + MDP(28) * t151 + t150 * t93) * t96 + (0.2e1 * t136 + t153) * qJ(5); t87 * MDP(21) + t126 * t94 + t130; t140; 0; t123; MDP(21); MDP(26) * t95 + t156; -t81 * MDP(27) + t83 * MDP(28); t129; t125; 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
