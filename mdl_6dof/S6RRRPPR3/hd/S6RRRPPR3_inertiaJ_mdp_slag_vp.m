% Calculate joint inertia matrix for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPPR3_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:30:09
% EndTime: 2019-03-09 15:30:11
% DurationCPUTime: 0.44s
% Computational Cost: add. (559->151), mult. (921->189), div. (0->0), fcn. (880->6), ass. (0->63)
t129 = (-pkin(3) - pkin(4));
t123 = sin(qJ(6));
t120 = t123 ^ 2;
t126 = cos(qJ(6));
t143 = t123 * t126 * MDP(27);
t142 = t120 * MDP(26) + MDP(15) + 0.2e1 * t143;
t160 = 2 * MDP(23);
t161 = 2 * MDP(18);
t165 = (pkin(3) * t161) + (t129 * t160) + t142;
t146 = (MDP(20) + MDP(22));
t164 = 2 * t146;
t150 = MDP(31) * t126;
t136 = -t123 * MDP(32) + t150;
t124 = sin(qJ(3));
t157 = pkin(2) * t124;
t110 = qJ(4) + t157;
t163 = t110 ^ 2;
t128 = cos(qJ(2));
t162 = 0.2e1 * t128;
t159 = (-pkin(4) - pkin(9));
t158 = -pkin(8) - pkin(7);
t116 = -t128 * pkin(2) - pkin(1);
t125 = sin(qJ(2));
t127 = cos(qJ(3));
t101 = t124 * t125 - t127 * t128;
t104 = t158 * t125;
t105 = t158 * t128;
t94 = t104 * t124 - t105 * t127;
t85 = qJ(5) * t101 + t94;
t156 = t101 * t85;
t155 = qJ(4) * t101;
t154 = qJ(4) * t110;
t153 = t101 * t110;
t152 = t101 * t126;
t102 = t124 * t128 + t125 * t127;
t149 = t102 * MDP(29);
t147 = t124 * MDP(17);
t145 = MDP(23) - MDP(18);
t144 = 0.2e1 * t102;
t114 = pkin(2) * t127 + pkin(3);
t141 = -t101 * pkin(3) + t102 * qJ(4) - t116;
t93 = -t127 * t104 - t105 * t124;
t109 = -pkin(4) - t114;
t107 = -pkin(9) + t109;
t108 = pkin(5) + t110;
t140 = t101 * t108 - t102 * t107;
t119 = -pkin(3) + t159;
t122 = qJ(4) + pkin(5);
t139 = t101 * t122 - t102 * t119;
t138 = MDP(28) * t126 - MDP(29) * t123;
t137 = -MDP(28) * t123 - MDP(29) * t126;
t135 = -MDP(31) * t123 - MDP(32) * t126;
t133 = 0.2e1 * t136;
t132 = -MDP(26) * t152 - t102 * MDP(28) - t85 * MDP(32);
t121 = t126 ^ 2;
t84 = -qJ(5) * t102 + t93;
t131 = t102 * MDP(13) + t84 * MDP(23) + (MDP(20) - MDP(17)) * t94 + (-MDP(18) - MDP(16)) * t93 + (t150 + MDP(22)) * t85 + ((t120 - t121) * MDP(27) - MDP(14)) * t101;
t130 = qJ(4) ^ 2;
t80 = -pkin(4) * t101 + t141;
t79 = pkin(5) * t102 + t159 * t101 + t141;
t78 = t123 * t79 + t126 * t84;
t77 = -t123 * t84 + t126 * t79;
t1 = [pkin(1) * MDP(9) * t162 + (t141 ^ 2 + t93 ^ 2 + t94 ^ 2) * MDP(21) + (t80 ^ 2 + t84 ^ 2 + t85 ^ 2) * MDP(25) + MDP(1) + (t121 * MDP(26) - 0.2e1 * t143) * t101 ^ 2 + (MDP(11) + MDP(30)) * t102 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t125 + MDP(5) * t162) * t125 + (t116 * MDP(17) + MDP(20) * t141 + t80 * MDP(22)) * t144 + (0.2e1 * t116 * MDP(16) - t141 * t161 + t80 * t160 + (-MDP(12) + t138) * t144) * t101 + 0.2e1 * (-t101 * t94 + t102 * t93) * MDP(19) + 0.2e1 * (-t102 * t84 + t156) * MDP(24) + 0.2e1 * (t102 * t77 + t123 * t156) * MDP(31) + 0.2e1 * (-t78 * t102 + t85 * t152) * MDP(32); t131 + (-MDP(10) * t128 - MDP(9) * t125) * pkin(7) + t125 * MDP(6) + t128 * MDP(7) + (-t102 * t109 + t153) * MDP(24) + (t109 * t84 + t110 * t85) * MDP(25) + (-t102 * t114 - t153) * MDP(19) + (t110 * t94 - t114 * t93) * MDP(21) + (t140 * MDP(31) + t132) * t123 + (t140 * MDP(32) - t149) * t126; MDP(8) + (t114 ^ 2 + t163) * MDP(21) + (t109 ^ 2 + t163) * MDP(25) + t108 * t133 + 0.2e1 * (MDP(16) * t127 - t147) * pkin(2) + t114 * t161 + t109 * t160 + t110 * t164 + t142; t131 + (t139 * MDP(31) + t132) * t123 + (t139 * MDP(32) - t149) * t126 + (-t102 * t129 + t155) * MDP(24) + (qJ(4) * t85 + t129 * t84) * MDP(25) + (-pkin(3) * t102 - t155) * MDP(19) + (-pkin(3) * t93 + qJ(4) * t94) * MDP(21); (pkin(3) * t114 + t154) * MDP(21) + (t109 * t129 + t154) * MDP(25) + t146 * (0.2e1 * qJ(4) + t157) + (-t147 + (MDP(16) - t145) * t127) * pkin(2) + t136 * (t108 + t122) + t165; ((pkin(3) ^ 2) + t130) * MDP(21) + ((t129 ^ 2) + t130) * MDP(25) + t122 * t133 + qJ(4) * t164 + t165; t93 * MDP(21) + MDP(25) * t84 + (MDP(19) - MDP(24) + t135) * t102; -MDP(21) * t114 + MDP(25) * t109 + t145; -pkin(3) * MDP(21) + MDP(25) * t129 + t145; MDP(21) + MDP(25); t101 * MDP(23) + MDP(25) * t80 + (MDP(22) + t136) * t102; 0; 0; 0; MDP(25); MDP(30) * t102 + t77 * MDP(31) - t78 * MDP(32) + t138 * t101; t135 * t107 + t137; t135 * t119 + t137; t135; t136; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
