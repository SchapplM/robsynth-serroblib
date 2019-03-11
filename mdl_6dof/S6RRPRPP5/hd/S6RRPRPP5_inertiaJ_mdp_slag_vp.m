% Calculate joint inertia matrix for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPRPP5_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:06:18
% EndTime: 2019-03-09 10:06:20
% DurationCPUTime: 0.61s
% Computational Cost: add. (595->181), mult. (914->219), div. (0->0), fcn. (722->4), ass. (0->66)
t156 = pkin(3) + pkin(7);
t117 = sin(qJ(4));
t155 = 0.2e1 * t117;
t142 = MDP(22) + MDP(26);
t119 = cos(qJ(4));
t121 = pkin(4) + pkin(5);
t150 = qJ(5) * t117;
t97 = t119 * t121 + t150;
t154 = (2 * pkin(4) * MDP(22)) + 0.2e1 * t121 * MDP(26) + MDP(19);
t123 = -pkin(2) - pkin(8);
t118 = sin(qJ(2));
t103 = t156 * t118;
t120 = cos(qJ(2));
t137 = -qJ(3) * t118 - pkin(1);
t92 = t120 * t123 + t137;
t152 = -t119 * t103 + t117 * t92;
t85 = t117 * t103 + t119 * t92;
t109 = t118 * qJ(5);
t82 = t109 + t85;
t151 = t117 * t82;
t149 = t117 * t120;
t148 = qJ(6) + t123;
t104 = t156 * t120;
t112 = t117 ^ 2;
t114 = t119 ^ 2;
t106 = t112 + t114;
t147 = MDP(16) * t119;
t146 = MDP(25) * t123;
t145 = t117 * MDP(26);
t144 = pkin(7) ^ 2 * MDP(14);
t143 = -MDP(21) + MDP(24);
t141 = -MDP(23) + MDP(28);
t140 = MDP(24) + MDP(27);
t139 = MDP(25) + MDP(29);
t138 = 0.2e1 * t109 + t85;
t136 = pkin(7) * MDP(14) + MDP(11);
t135 = -pkin(2) * MDP(14) + MDP(12);
t133 = pkin(4) * MDP(25) + MDP(22);
t132 = qJ(5) * t119 - qJ(3);
t102 = pkin(4) * t119 + t150;
t131 = -t152 * MDP(20) - t85 * MDP(21);
t100 = pkin(4) * t117 - t132;
t91 = -t117 * t121 + t132;
t130 = t100 * MDP(22) - t91 * MDP(26);
t129 = -t100 * MDP(24) + t91 * MDP(27);
t96 = t148 * t117;
t98 = t148 * t119;
t128 = t119 * MDP(17) + t98 * MDP(26) + t96 * MDP(27);
t79 = qJ(6) * t149 - t118 * t121 + t152;
t86 = -t120 * t97 - t104;
t88 = t102 * t120 + t104;
t127 = -t104 * MDP(21) + t88 * MDP(24) - t86 * MDP(27) + t79 * MDP(28);
t126 = (MDP(20) + t142) * t119 + (-MDP(21) + t140) * t117;
t124 = qJ(5) ^ 2;
t115 = t120 ^ 2;
t113 = t118 ^ 2;
t107 = t119 * t120 * qJ(6);
t105 = t119 * t123 * t118;
t101 = -pkin(2) * t120 + t137;
t95 = t106 * t123;
t87 = t117 * t96 + t119 * t98;
t83 = -pkin(4) * t118 + t152;
t81 = t107 + t82;
t80 = t81 * t117;
t78 = -t119 * t83 + t151;
t1 = [(t82 ^ 2 + t83 ^ 2 + t88 ^ 2) * MDP(25) + (t79 ^ 2 + t81 ^ 2 + t86 ^ 2) * MDP(29) + MDP(1) + t101 ^ 2 * MDP(14) + (t112 * MDP(15) + t147 * t155 + t144) * t115 + (MDP(19) + MDP(4) + t144) * t113 + 0.2e1 * (t113 + t115) * MDP(11) * pkin(7) + 0.2e1 * (pkin(1) * MDP(9) + t101 * MDP(12) + (t104 * MDP(20) + t88 * MDP(22) - t82 * MDP(23) - t86 * MDP(26) + t81 * MDP(28)) * t119 + (-t83 * MDP(23) + t127) * t117) * t120 + 0.2e1 * (-pkin(1) * MDP(10) - t101 * MDP(13) + (-t117 * MDP(17) - t119 * MDP(18) + MDP(5)) * t120 - t83 * MDP(22) + t82 * MDP(24) - t79 * MDP(26) + t81 * MDP(27) + t131) * t118; (t104 * t117 + t105) * MDP(20) + (t117 * t88 + t105) * MDP(22) - t78 * MDP(23) + (t100 * t88 + t123 * t151) * MDP(25) - t86 * t145 + t80 * MDP(28) + (-t79 * t98 + t81 * t96 + t86 * t91) * MDP(29) + (-t83 * t146 - t127) * t119 + (-pkin(2) * MDP(11) + MDP(6) + (t123 * t143 - MDP(18)) * t117 + (-MDP(9) + t135) * pkin(7) + t128) * t118 + (MDP(7) + (t112 - t114) * MDP(16) + (-MDP(10) + MDP(13)) * pkin(7) + (t96 * MDP(28) + t130) * t119 + (-t119 * MDP(15) - t98 * MDP(28) - t129) * t117 + (t119 * MDP(20) - t117 * MDP(21) + t136) * qJ(3)) * t120; MDP(8) - 0.2e1 * pkin(2) * MDP(12) + 0.2e1 * qJ(3) * MDP(13) + (pkin(2) ^ 2 + qJ(3) ^ 2) * MDP(14) + t114 * MDP(15) + (t106 * t123 ^ 2 + t100 ^ 2) * MDP(25) + (t91 ^ 2 + t96 ^ 2 + t98 ^ 2) * MDP(29) - 0.2e1 * MDP(23) * t95 + 0.2e1 * MDP(28) * t87 + 0.2e1 * (qJ(3) * MDP(21) + t129) * t119 + (qJ(3) * MDP(20) + t130 - t147) * t155; t78 * MDP(25) + (-t119 * t79 + t80) * MDP(29) + (t126 + t136) * t118; MDP(25) * t95 + MDP(29) * t87 + t106 * t141 + t135; t106 * t139 + MDP(14); t138 * MDP(24) + (-pkin(4) * t83 + qJ(5) * t82) * MDP(25) + (t107 + t138) * MDP(27) + (qJ(5) * t81 - t121 * t79) * MDP(29) + t154 * t118 + ((qJ(5) * t141 - MDP(18)) * t119 + (MDP(23) * pkin(4) - MDP(26) * qJ(6) - MDP(28) * t121 - MDP(17)) * t117) * t120 + t131 - t142 * t152; -t117 * MDP(18) - t102 * MDP(23) + t97 * MDP(28) + (qJ(5) * t96 + t121 * t98) * MDP(29) + ((MDP(20) + t133) * t119 + (MDP(25) * qJ(5) + t143) * t117) * t123 + t128; MDP(25) * t102 + MDP(29) * t97 + t126; ((pkin(4) ^ 2) + t124) * MDP(25) + (t121 ^ 2 + t124) * MDP(29) + 0.2e1 * t140 * qJ(5) + t154; MDP(25) * t83 + t79 * MDP(29) - t118 * t142 + t141 * t149; -t98 * MDP(29) + (-t141 - t146) * t119; -t139 * t119; -MDP(29) * t121 - MDP(26) - t133; t139; t86 * MDP(29) + (-t119 * MDP(26) - t117 * MDP(27)) * t120; t119 * MDP(27) + t91 * MDP(29) - t145; 0; 0; 0; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
