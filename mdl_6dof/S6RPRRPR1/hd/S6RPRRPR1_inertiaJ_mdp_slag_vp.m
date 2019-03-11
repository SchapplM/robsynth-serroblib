% Calculate joint inertia matrix for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRPR1_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:58:46
% EndTime: 2019-03-09 04:58:48
% DurationCPUTime: 0.41s
% Computational Cost: add. (769->116), mult. (1373->179), div. (0->0), fcn. (1543->10), ass. (0->64)
t120 = sin(qJ(6));
t123 = cos(qJ(6));
t136 = t123 * MDP(26) - t120 * MDP(27);
t122 = sin(qJ(3));
t121 = sin(qJ(4));
t124 = cos(qJ(3));
t149 = t121 * t124;
t157 = cos(qJ(4));
t101 = t157 * t122 + t149;
t116 = sin(pkin(11));
t118 = cos(pkin(11));
t150 = t121 * t122;
t134 = -t157 * t124 + t150;
t88 = t116 * t101 + t118 * t134;
t160 = t136 * t88;
t147 = t120 * MDP(23) + t123 * MDP(24);
t159 = t134 * MDP(17);
t119 = cos(pkin(10));
t144 = -t119 * pkin(1) - pkin(2);
t102 = -t124 * pkin(3) + t144;
t158 = 0.2e1 * t102;
t156 = pkin(3) * t121;
t117 = sin(pkin(10));
t107 = t117 * pkin(1) + pkin(7);
t155 = pkin(8) + t107;
t154 = MDP(20) * pkin(4);
t141 = t155 * t122;
t128 = -t157 * t141 - t155 * t149;
t127 = -t101 * qJ(5) + t128;
t140 = t157 * t155;
t80 = (t157 * qJ(5) + t140) * t124 + (-qJ(5) - t155) * t150;
t74 = t116 * t80 - t118 * t127;
t73 = t74 * t120;
t153 = t74 * t123;
t151 = t120 * t123;
t148 = t88 * MDP(25);
t110 = t157 * pkin(3) + pkin(4);
t95 = t116 * t110 + t118 * t156;
t143 = MDP(22) * t151;
t114 = t120 ^ 2;
t142 = t114 * MDP(21) + MDP(16) + 0.2e1 * t143;
t90 = t118 * t101 - t116 * t134;
t94 = t118 * t110 - t116 * t156;
t92 = -pkin(5) - t94;
t93 = pkin(9) + t95;
t139 = -t88 * t93 + t90 * t92;
t106 = t116 * pkin(4) + pkin(9);
t108 = -t118 * pkin(4) - pkin(5);
t138 = -t106 * t88 + t108 * t90;
t137 = -t124 * MDP(10) + t122 * MDP(11);
t135 = -MDP(26) * t120 - MDP(27) * t123;
t133 = -t101 * MDP(18) - t159 - t160;
t132 = (MDP(23) * t123 - MDP(24) * t120) * t90;
t115 = t123 ^ 2;
t131 = t128 * MDP(17) + (t121 * t141 - t124 * t140) * MDP(18) - t134 * MDP(15) + t101 * MDP(14) + ((-t114 + t115) * MDP(22) + MDP(21) * t151) * t90 + t147 * t88;
t130 = 0.2e1 * t136;
t129 = (t157 * MDP(17) - t121 * MDP(18)) * pkin(3);
t91 = t134 * pkin(4) + t102;
t87 = t90 ^ 2;
t77 = t88 * pkin(5) - t90 * pkin(9) + t91;
t76 = t116 * t127 + t118 * t80;
t72 = t120 * t77 + t123 * t76;
t71 = -t120 * t76 + t123 * t77;
t1 = [MDP(1) + t158 * t159 + (t74 ^ 2 + t76 ^ 2 + t91 ^ 2) * MDP(20) + (t115 * MDP(21) - 0.2e1 * t143) * t87 + (t117 ^ 2 + t119 ^ 2) * MDP(4) * pkin(1) ^ 2 + t148 * t88 + (MDP(12) * t101 - 0.2e1 * t134 * MDP(13) + MDP(18) * t158) * t101 + 0.2e1 * (t74 * t90 - t76 * t88) * MDP(19) + 0.2e1 * (t71 * t88 + t90 * t73) * MDP(26) + 0.2e1 * (t90 * t153 - t72 * t88) * MDP(27) + 0.2e1 * t132 * t88 + 0.2e1 * t137 * t144 + (MDP(5) * t122 + 0.2e1 * t124 * MDP(6)) * t122; (t74 * t88 + t76 * t90) * MDP(20); MDP(4) + (t88 ^ 2 + t87) * MDP(20); t122 * MDP(7) + t124 * MDP(8) + (-t95 * t88 - t94 * t90) * MDP(19) + (-t74 * t94 + t76 * t95) * MDP(20) + (t139 * t120 - t153) * MDP(26) + (t139 * t123 + t73) * MDP(27) + (-t122 * MDP(10) - t124 * MDP(11)) * t107 + t131; (-t88 * t94 + t90 * t95) * MDP(20) + t133 - t137; MDP(9) + (t94 ^ 2 + t95 ^ 2) * MDP(20) - t92 * t130 + 0.2e1 * t129 + t142; (t138 * t120 - t153) * MDP(26) + (t138 * t123 + t73) * MDP(27) + ((-t116 * t88 - t118 * t90) * MDP(19) + (t116 * t76 - t118 * t74) * MDP(20)) * pkin(4) + t131; (t116 * t90 - t118 * t88) * t154 + t133; (t116 * t95 + t118 * t94) * t154 + t129 + t142 - t136 * (t108 + t92); (t116 ^ 2 + t118 ^ 2) * MDP(20) * pkin(4) ^ 2 - t108 * t130 + t142; t91 * MDP(20) + t160; 0; 0; 0; MDP(20); t71 * MDP(26) - t72 * MDP(27) + t132 + t148; t135 * t90; t135 * t93 + t147; t135 * t106 + t147; t136; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
