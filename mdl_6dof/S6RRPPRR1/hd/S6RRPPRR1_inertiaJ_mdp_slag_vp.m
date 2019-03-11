% Calculate joint inertia matrix for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPPRR1_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:48:04
% EndTime: 2019-03-09 08:48:05
% DurationCPUTime: 0.43s
% Computational Cost: add. (663->129), mult. (1183->179), div. (0->0), fcn. (1310->8), ass. (0->65)
t119 = sin(qJ(6));
t122 = cos(qJ(6));
t161 = MDP(29) * t119 + MDP(30) * t122;
t120 = sin(qJ(5));
t123 = cos(qJ(5));
t143 = t122 * MDP(29);
t131 = -t119 * MDP(30) + t143;
t129 = MDP(22) + t131;
t160 = -t120 * MDP(23) + t123 * t129;
t117 = sin(pkin(10));
t118 = cos(pkin(10));
t121 = sin(qJ(2));
t124 = cos(qJ(2));
t100 = t117 * t121 - t118 * t124;
t101 = t117 * t124 + t118 * t121;
t112 = -t124 * pkin(2) - pkin(1);
t85 = t100 * pkin(3) - t101 * qJ(4) + t112;
t81 = -t100 * pkin(4) - t85;
t159 = 0.2e1 * t81;
t158 = 0.2e1 * t124;
t110 = t118 * pkin(2) + pkin(3);
t107 = -pkin(4) - t110;
t108 = t117 * pkin(2) + qJ(4);
t95 = -t123 * t107 + t120 * t108;
t93 = pkin(5) + t95;
t157 = pkin(5) + t93;
t156 = -qJ(3) - pkin(7);
t88 = t120 * t100 + t123 * t101;
t155 = t122 * t88;
t154 = t85 * MDP(16);
t87 = -t123 * t100 + t120 * t101;
t153 = t87 * MDP(28);
t152 = t88 * MDP(23);
t151 = t95 * MDP(22);
t96 = t120 * t107 + t123 * t108;
t150 = t96 * MDP(23);
t105 = t156 * t124;
t138 = t156 * t121;
t92 = -t118 * t105 + t117 * t138;
t147 = t100 * MDP(13);
t146 = t101 * MDP(15);
t144 = t122 * MDP(25);
t115 = t119 ^ 2;
t142 = t115 * MDP(24) + MDP(21);
t90 = -t117 * t105 - t118 * t138;
t141 = t90 ^ 2 + t92 ^ 2;
t140 = t119 * t144;
t139 = 0.2e1 * t140 + t142;
t136 = -pkin(5) * t88 - pkin(9) * t87;
t94 = -pkin(9) + t96;
t135 = -t87 * t94 + t88 * t93;
t82 = -t101 * pkin(8) + t90;
t83 = t100 * pkin(8) + t92;
t78 = t120 * t83 - t123 * t82;
t133 = -t87 * MDP(27) + t78 * MDP(29);
t132 = MDP(26) * t122 - MDP(27) * t119;
t128 = 0.2e1 * t131;
t127 = -MDP(24) * t155 - t87 * MDP(26) - t78 * MDP(30);
t116 = t122 ^ 2;
t79 = t120 * t82 + t123 * t83;
t126 = -t87 * MDP(20) - t78 * MDP(22) - t79 * MDP(23) + (MDP(19) - (t115 - t116) * MDP(25)) * t88;
t77 = t87 * pkin(5) - t88 * pkin(9) + t81;
t76 = t119 * t77 + t122 * t79;
t75 = -t119 * t79 + t122 * t77;
t1 = [MDP(1) + pkin(1) * MDP(9) * t158 + (t112 ^ 2 + t141) * MDP(12) + t141 * MDP(16) + t152 * t159 + (-0.2e1 * t146 + 0.2e1 * t147 + t154) * t85 + (t116 * MDP(24) + MDP(17) - 0.2e1 * t140) * t88 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t121 + MDP(5) * t158) * t121 + (MDP(22) * t159 + t153 + 0.2e1 * (-MDP(18) + t132) * t88) * t87 + 0.2e1 * (t78 * t119 * t88 + t75 * t87) * MDP(29) + 0.2e1 * (t78 * t155 - t76 * t87) * MDP(30) + 0.2e1 * (MDP(11) + MDP(14)) * (-t92 * t100 + t90 * t101); t121 * MDP(6) + t124 * MDP(7) - t90 * MDP(13) + (-t108 * t100 - t110 * t101) * MDP(14) + t92 * MDP(15) + (t92 * t108 - t90 * t110) * MDP(16) + (-t124 * MDP(10) - t121 * MDP(9)) * pkin(7) + (t135 * MDP(30) + t133) * t122 + (t135 * MDP(29) + t127) * t119 + ((-t100 * t117 - t101 * t118) * MDP(11) + (t117 * t92 - t118 * t90) * MDP(12)) * pkin(2) - t126; MDP(8) + (t108 ^ 2 + t110 ^ 2) * MDP(16) + t93 * t128 + (t117 ^ 2 + t118 ^ 2) * MDP(12) * pkin(2) ^ 2 + 0.2e1 * t110 * MDP(13) + 0.2e1 * t108 * MDP(15) + 0.2e1 * t151 + 0.2e1 * t150 + t139; t112 * MDP(12) - t129 * t87 - t146 + t147 - t152 + t154; 0; MDP(12) + MDP(16); t101 * MDP(14) + t90 * MDP(16) + t161 * (-t120 * t87 - t123 * t88); -t110 * MDP(16) - MDP(13) - t160; 0; MDP(16); (t136 * MDP(30) - t133) * t122 + (t136 * MDP(29) - t127) * t119 + t126; -t151 - t150 - t157 * t143 + (t157 * MDP(30) - 0.2e1 * t144) * t119 - t142; 0; t160; pkin(5) * t128 + t139; t75 * MDP(29) - t76 * MDP(30) + t132 * t88 + t153; (-t94 * MDP(30) - MDP(27)) * t122 + (-t94 * MDP(29) - MDP(26)) * t119; -t131; -t161 * t120; t119 * MDP(26) + t122 * MDP(27) - pkin(9) * t161; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
