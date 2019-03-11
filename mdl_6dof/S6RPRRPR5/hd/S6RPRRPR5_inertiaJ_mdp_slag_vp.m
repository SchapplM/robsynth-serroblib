% Calculate joint inertia matrix for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRPR5_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:13:35
% EndTime: 2019-03-09 05:13:37
% DurationCPUTime: 0.42s
% Computational Cost: add. (808->133), mult. (1477->175), div. (0->0), fcn. (1683->8), ass. (0->68)
t145 = MDP(20) - MDP(23);
t144 = -MDP(21) + MDP(24);
t126 = cos(qJ(6));
t147 = MDP(32) * t126;
t123 = sin(qJ(6));
t149 = MDP(31) * t123;
t166 = t147 + t149;
t122 = cos(pkin(10));
t113 = -pkin(2) * t122 - pkin(1);
t165 = 0.2e1 * t113;
t164 = -2 * MDP(23);
t163 = 2 * MDP(24);
t162 = pkin(4) + pkin(9);
t161 = pkin(1) * MDP(7);
t160 = pkin(7) + qJ(2);
t159 = (MDP(25) * pkin(4));
t121 = sin(pkin(10));
t125 = sin(qJ(3));
t128 = cos(qJ(3));
t104 = t121 * t125 - t128 * t122;
t105 = t121 * t128 + t122 * t125;
t124 = sin(qJ(4));
t127 = cos(qJ(4));
t99 = t127 * t104 + t105 * t124;
t158 = qJ(5) * t99;
t111 = pkin(3) * t124 + qJ(5);
t157 = t111 * t99;
t106 = t160 * t121;
t107 = t160 * t122;
t140 = -t128 * t106 - t107 * t125;
t91 = -pkin(8) * t105 + t140;
t137 = t106 * t125 - t107 * t128;
t92 = -pkin(8) * t104 - t137;
t88 = t124 * t91 + t127 * t92;
t81 = -pkin(5) * t99 + t88;
t77 = t81 * t123;
t78 = t81 * t126;
t156 = MDP(4) * t122;
t155 = MDP(5) * t121;
t154 = t123 * t126;
t151 = MDP(13) * t104;
t114 = -pkin(3) * t127 - pkin(4);
t150 = MDP(25) * t114;
t148 = MDP(31) * t126;
t100 = -t104 * t124 + t105 * t127;
t146 = t100 * MDP(29);
t116 = t126 * MDP(28);
t143 = 0.2e1 * t99;
t142 = MDP(27) * t154;
t120 = t126 ^ 2;
t141 = t120 * MDP(26) + MDP(19) - 0.2e1 * t142;
t87 = t124 * t92 - t127 * t91;
t139 = t100 * t162 + t158;
t109 = -pkin(9) + t114;
t138 = t100 * t109 - t157;
t101 = pkin(3) * t104 + t113;
t136 = MDP(28) * t123 + MDP(29) * t126;
t135 = -t123 * MDP(32) + t148;
t133 = t163 + 0.2e1 * t147 + 0.2e1 * t149;
t132 = -qJ(5) * t100 + t101;
t119 = t123 ^ 2;
t131 = t78 * MDP(32) + t144 * t88 - t145 * t87 + (t116 + MDP(17)) * t100 + ((-t119 + t120) * MDP(27) + MDP(26) * t154 - MDP(18)) * t99;
t86 = pkin(4) * t99 + t132;
t80 = pkin(5) * t100 + t87;
t79 = t162 * t99 + t132;
t76 = t123 * t80 + t126 * t79;
t75 = -t123 * t79 + t126 * t80;
t1 = [(t86 ^ 2 + t87 ^ 2 + t88 ^ 2) * MDP(25) + t151 * t165 + MDP(1) + (MDP(20) * t101 - MDP(23) * t86) * t143 + (MDP(15) + MDP(30)) * t100 ^ 2 + (MDP(26) * t119 + 0.2e1 * t142) * t99 ^ 2 + (-0.2e1 * t155 + 0.2e1 * t156 + t161) * pkin(1) + (MDP(14) * t165 + MDP(8) * t105 - 0.2e1 * MDP(9) * t104) * t105 + (0.2e1 * t101 * MDP(21) - 0.2e1 * t86 * MDP(24) + (-MDP(16) + t136) * t143) * t100 + 0.2e1 * (t100 * t87 - t88 * t99) * MDP(22) + 0.2e1 * (t100 * t75 - t99 * t78) * MDP(31) + 0.2e1 * (-t100 * t76 + t99 * t77) * MDP(32) + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t121 ^ 2 + t122 ^ 2) * qJ(2); t151 + MDP(14) * t105 + MDP(25) * t86 - t156 + t155 - t161 + t145 * t99 + (-t166 - t144) * t100; MDP(7) + MDP(25); t140 * MDP(13) + t137 * MDP(14) + (t111 * t88 + t114 * t87) * MDP(25) + (t100 * t114 - t157) * MDP(22) - t104 * MDP(11) + t105 * MDP(10) + t131 + (-t138 * MDP(32) - t146) * t123 + (t138 * t126 + t77) * MDP(31); 0; MDP(12) + ((2 * MDP(23)) + t150) * t114 + 0.2e1 * (t127 * MDP(20) - t124 * MDP(21)) * pkin(3) + (MDP(25) * t111 + t133) * t111 + t141; (t139 * MDP(32) - t146) * t123 + (-pkin(4) * t100 - t158) * MDP(22) + (-pkin(4) * t87 + qJ(5) * t88) * MDP(25) + t131 + (-t139 * t126 + t77) * MDP(31); 0; pkin(4) * t164 + qJ(5) * t163 + (-pkin(4) * t114 + qJ(5) * t111) * MDP(25) + (t144 * t124 + t145 * t127) * pkin(3) + t141 + t166 * (qJ(5) + t111); (t164 + t159) * pkin(4) + (MDP(25) * qJ(5) + t133) * qJ(5) + t141; MDP(25) * t87 + (MDP(22) + t135) * t100; 0; MDP(23) + t150; MDP(23) - t159; MDP(25); MDP(30) * t100 + MDP(31) * t75 - MDP(32) * t76 + t136 * t99; -t166; -MDP(29) * t123 + t135 * t109 + t116; -t162 * t148 + t116 + (MDP(32) * t162 - MDP(29)) * t123; t135; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
