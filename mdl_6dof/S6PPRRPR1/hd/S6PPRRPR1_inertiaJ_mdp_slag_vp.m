% Calculate joint inertia matrix for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRPR1_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:47:29
% EndTime: 2019-03-08 18:47:31
% DurationCPUTime: 0.45s
% Computational Cost: add. (570->142), mult. (1452->232), div. (0->0), fcn. (1716->14), ass. (0->73)
t146 = (MDP(16) * qJ(5));
t161 = MDP(15) + t146;
t160 = 2 * MDP(15);
t159 = -2 * MDP(18);
t158 = 2 * MDP(23);
t113 = sin(pkin(13));
t157 = pkin(9) * t113;
t125 = cos(qJ(4));
t156 = pkin(9) * t125;
t122 = sin(qJ(4));
t155 = pkin(10) * t122;
t154 = pkin(10) + qJ(5);
t153 = pkin(4) * MDP(16);
t103 = -pkin(4) * t125 - qJ(5) * t122 - pkin(3);
t117 = cos(pkin(13));
t92 = t103 * t113 + t117 * t156;
t115 = sin(pkin(7));
t123 = sin(qJ(3));
t152 = t115 * t123;
t126 = cos(qJ(3));
t151 = t115 * t126;
t118 = cos(pkin(12));
t119 = cos(pkin(7));
t150 = t118 * t119;
t121 = sin(qJ(6));
t124 = cos(qJ(6));
t101 = t113 * t124 + t117 * t121;
t93 = t101 * t122;
t149 = t93 * MDP(20);
t100 = t113 * t121 - t117 * t124;
t94 = t100 * t122;
t148 = t94 * MDP(19);
t145 = MDP(17) * t101;
t144 = t100 * MDP(22);
t143 = t113 * MDP(14);
t142 = t117 * MDP(13);
t141 = t122 * MDP(12);
t140 = t125 * MDP(21);
t114 = sin(pkin(12));
t116 = sin(pkin(6));
t120 = cos(pkin(6));
t85 = t120 * t152 + (t114 * t126 + t123 * t150) * t116;
t95 = -t115 * t116 * t118 + t119 * t120;
t82 = t122 * t95 + t125 * t85;
t84 = -t120 * t151 + (t114 * t123 - t126 * t150) * t116;
t75 = -t113 * t82 + t117 * t84;
t76 = t113 * t84 + t117 * t82;
t136 = (-t121 * t76 + t124 * t75) * MDP(22) - (t121 * t75 + t124 * t76) * MDP(23);
t97 = t119 * t122 + t125 * t152;
t86 = -t113 * t97 - t117 * t151;
t87 = -t113 * t151 + t117 * t97;
t135 = (-t121 * t87 + t124 * t86) * MDP(22) - (t121 * t86 + t124 * t87) * MDP(23);
t134 = t93 * MDP(22) - t94 * MDP(23);
t133 = MDP(13) * t113 + MDP(14) * t117;
t132 = -t142 + t143 - t153;
t131 = MDP(16) * pkin(9) + t133;
t104 = t154 * t113;
t105 = t154 * t117;
t130 = t101 * MDP(19) - t100 * MDP(20) + (-t104 * t124 - t105 * t121) * MDP(22) - (-t104 * t121 + t105 * t124) * MDP(23);
t129 = MDP(23) * t101 + t132 + t144;
t128 = -MDP(11) + t129;
t112 = t122 ^ 2;
t108 = -pkin(5) * t117 - pkin(4);
t102 = (pkin(5) * t113 + pkin(9)) * t122;
t99 = t117 * t103;
t96 = -t119 * t125 + t122 * t152;
t91 = -t113 * t156 + t99;
t88 = -t113 * t155 + t92;
t83 = -t117 * t155 + t99 + (-pkin(5) - t157) * t125;
t81 = t122 * t85 - t125 * t95;
t78 = t121 * t83 + t124 * t88;
t77 = -t121 * t88 + t124 * t83;
t1 = [MDP(1) + (t120 ^ 2 + (t114 ^ 2 + t118 ^ 2) * t116 ^ 2) * MDP(2) + (t75 ^ 2 + t76 ^ 2 + t81 ^ 2) * MDP(16); t120 * MDP(2) + (t75 * t86 + t76 * t87 + t81 * t96) * MDP(16); MDP(2) + (t86 ^ 2 + t87 ^ 2 + t96 ^ 2) * MDP(16); -t84 * MDP(4) - t85 * MDP(5) + (t75 * t91 + t76 * t92) * MDP(16) + t134 * t81 + (-MDP(11) * t84 - MDP(13) * t75 + MDP(14) * t76 - t136) * t125 + (t84 * MDP(12) + (-t113 * t76 - t117 * t75) * MDP(15) + t131 * t81) * t122; (t86 * t91 + t87 * t92) * MDP(16) + t134 * t96 + (-MDP(13) * t86 + MDP(14) * t87 - t135) * t125 + ((-t113 * t87 - t117 * t86) * MDP(15) + t131 * t96) * t122 + (-t123 * MDP(5) + (MDP(11) * t125 + MDP(4) - t141) * t126) * t115; MDP(3) + t112 * MDP(6) - 0.2e1 * pkin(3) * t141 + (pkin(9) ^ 2 * t112 + t91 ^ 2 + t92 ^ 2) * MDP(16) - (-MDP(17) * t94 + t159 * t93) * t94 + (0.2e1 * MDP(11) * pkin(3) + 0.2e1 * MDP(7) * t122 + t140 + 0.2e1 * t148 + 0.2e1 * t149) * t125 + 0.2e1 * (t112 * t157 - t125 * t91) * MDP(13) + 0.2e1 * (pkin(9) * t112 * t117 + t125 * t92) * MDP(14) + 0.2e1 * (t102 * t93 - t125 * t77) * MDP(22) + (-t102 * t94 + t125 * t78) * t158 + (-t113 * t92 - t117 * t91) * t122 * t160; -t82 * MDP(12) + t128 * t81 + t161 * (-t113 * t75 + t117 * t76); -t97 * MDP(12) + t128 * t96 + t161 * (-t113 * t86 + t117 * t87); -t94 * t145 + (t100 * t94 - t101 * t93) * MDP(18) + (t100 * t102 + t108 * t93) * MDP(22) + (t101 * t102 - t108 * t94) * MDP(23) + (-pkin(9) * MDP(12) + qJ(5) * t133 + MDP(9) - t130) * t125 + (MDP(8) - t133 * pkin(4) + (-MDP(11) + t132) * pkin(9)) * t122 + t161 * (-t113 * t91 + t117 * t92); 0.2e1 * t108 * t144 + MDP(10) + (0.2e1 * t142 - 0.2e1 * t143 + t153) * pkin(4) + (t100 * t159 + t108 * t158 + t145) * t101 + (t160 + t146) * (t113 ^ 2 + t117 ^ 2) * qJ(5); t81 * MDP(16); t96 * MDP(16); t122 * t131 + t134; t129; MDP(16); t136; t135; MDP(22) * t77 - MDP(23) * t78 - t140 - t148 - t149; t130; 0; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
