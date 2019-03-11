% Calculate joint inertia matrix for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPPR5_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:51:39
% EndTime: 2019-03-09 02:51:41
% DurationCPUTime: 0.43s
% Computational Cost: add. (719->127), mult. (1288->165), div. (0->0), fcn. (1414->8), ass. (0->65)
t124 = sin(qJ(6));
t126 = cos(qJ(6));
t120 = sin(pkin(9));
t122 = cos(pkin(9));
t125 = sin(qJ(3));
t127 = cos(qJ(3));
t104 = t120 * t127 + t122 * t125;
t119 = sin(pkin(10));
t102 = t120 * t125 - t127 * t122;
t158 = pkin(8) * t102;
t113 = -t122 * pkin(2) - pkin(1);
t133 = -t104 * qJ(4) + t113;
t157 = pkin(3) + qJ(5);
t85 = t157 * t102 + t133;
t121 = cos(pkin(10));
t156 = pkin(7) + qJ(2);
t107 = t156 * t120;
t108 = t156 * t122;
t95 = t107 * t127 + t108 * t125;
t88 = pkin(4) * t104 + t95;
t87 = t121 * t88;
t79 = t104 * pkin(5) + t87 + (-t85 - t158) * t119;
t82 = t119 * t88 + t121 * t85;
t80 = t121 * t158 + t82;
t103 = -t119 * t124 + t121 * t126;
t90 = t103 * t102;
t101 = t119 * t126 + t121 * t124;
t91 = t101 * t102;
t165 = MDP(28) * (-t124 * t80 + t126 * t79) - MDP(29) * (t124 * t79 + t126 * t80) + MDP(25) * t91 + MDP(26) * t90;
t149 = t119 ^ 2 + t121 ^ 2;
t164 = t149 * MDP(21);
t111 = pkin(5) * t119 + qJ(4);
t161 = 0.2e1 * t111;
t160 = -2 * MDP(24);
t159 = pkin(1) * MDP(7);
t155 = -pkin(8) - t157;
t153 = t120 * MDP(5);
t151 = t122 * MDP(4);
t81 = -t119 * t85 + t87;
t78 = t119 * t82 + t121 * t81;
t150 = t78 * MDP(22);
t147 = MDP(23) * t103;
t146 = MDP(28) * t101;
t145 = qJ(4) * MDP(22);
t144 = MDP(22) * t149 + MDP(18);
t143 = MDP(17) - MDP(14);
t142 = -pkin(3) * MDP(18) + MDP(16);
t141 = -t119 * t81 + t121 * t82;
t138 = -t90 * MDP(28) + t91 * MDP(29);
t96 = -t125 * t107 + t127 * t108;
t137 = MDP(19) * t121 - MDP(20) * t119;
t136 = MDP(19) * t119 + MDP(20) * t121;
t135 = MDP(28) * t103 - MDP(29) * t101;
t134 = -MDP(29) * t103 - t146;
t132 = -MDP(15) - t137;
t105 = t155 * t119;
t106 = t155 * t121;
t131 = t103 * MDP(25) - t101 * MDP(26) + (-t105 * t124 + t106 * t126) * MDP(28) - (t105 * t126 + t106 * t124) * MDP(29);
t130 = -t134 + t136;
t128 = qJ(4) ^ 2;
t100 = t149 * t157;
t92 = pkin(3) * t102 + t133;
t89 = -pkin(4) * t102 + t96;
t83 = (-pkin(5) * t121 - pkin(4)) * t102 + t96;
t1 = [(t92 ^ 2 + t95 ^ 2 + t96 ^ 2) * MDP(18) + (t81 ^ 2 + t82 ^ 2 + t89 ^ 2) * MDP(22) + MDP(1) + (MDP(8) + MDP(27)) * t104 ^ 2 + (MDP(23) * t91 - t90 * t160) * t91 + (0.2e1 * t151 - 0.2e1 * t153 + t159) * pkin(1) + 0.2e1 * t138 * t83 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t120 ^ 2 + t122 ^ 2) * qJ(2) + 0.2e1 * (t113 * MDP(13) - t96 * MDP(15) - t92 * MDP(16) + t141 * MDP(21) - t137 * t89) * t102 + 0.2e1 * (MDP(14) * t113 + MDP(15) * t95 - MDP(17) * t92 + MDP(19) * t81 - MDP(20) * t82 - MDP(9) * t102 + t165) * t104; -t151 + t153 - t159 + t92 * MDP(18) + t141 * MDP(22) + (MDP(13) - MDP(16) + t164) * t102 + (-t130 - t143) * t104; MDP(7) + t144; -t78 * MDP(21) + t91 * t147 + (-t101 * t91 + t103 * t90) * MDP(24) + (t101 * t83 - t111 * t90) * MDP(28) + (t103 * t83 + t111 * t91) * MDP(29) + (MDP(18) * qJ(4) + t143) * t96 + (-MDP(13) + t142) * t95 - t157 * t150 + (t136 + t145) * t89 + (t132 * qJ(4) - MDP(11)) * t102 + (-pkin(3) * MDP(15) - t137 * t157 + MDP(10) + t131) * t104; 0; MDP(12) - 0.2e1 * pkin(3) * MDP(16) + (pkin(3) ^ 2 + t128) * MDP(18) + 0.2e1 * t100 * MDP(21) + (t149 * t157 ^ 2 + t128) * MDP(22) + t146 * t161 + (MDP(29) * t161 + t101 * t160 + t147) * t103 + 0.2e1 * (MDP(17) + t136) * qJ(4); t95 * MDP(18) + t150 + (-t132 + t135) * t104; 0; -MDP(22) * t100 + t142 - t164; t144; t89 * MDP(22) - t137 * t102 + t138; 0; t130 + t145; 0; MDP(22); MDP(27) * t104 + t165; t134; t131; t135; 0; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
