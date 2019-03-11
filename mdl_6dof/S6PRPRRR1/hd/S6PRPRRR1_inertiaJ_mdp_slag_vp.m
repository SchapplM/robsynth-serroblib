% Calculate joint inertia matrix for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRRR1_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:25:07
% EndTime: 2019-03-08 20:25:08
% DurationCPUTime: 0.33s
% Computational Cost: add. (380->101), mult. (809->156), div. (0->0), fcn. (899->12), ass. (0->63)
t119 = sin(qJ(6));
t123 = cos(qJ(6));
t143 = t119 * MDP(22) + t123 * MDP(23);
t132 = t123 * MDP(25) - t119 * MDP(26);
t120 = sin(qJ(5));
t121 = sin(qJ(4));
t124 = cos(qJ(4));
t152 = cos(qJ(5));
t100 = t120 * t121 - t152 * t124;
t101 = t120 * t124 + t152 * t121;
t96 = t101 * MDP(19);
t154 = t100 * MDP(18) + t96;
t117 = cos(pkin(12));
t150 = t117 * pkin(2);
t106 = -pkin(3) - t150;
t102 = -pkin(4) * t124 + t106;
t153 = 0.2e1 * t102;
t115 = sin(pkin(12));
t151 = pkin(2) * t115;
t105 = pkin(8) + t151;
t148 = pkin(9) + t105;
t97 = t148 * t121;
t98 = t148 * t124;
t81 = t120 * t98 + t152 * t97;
t147 = t123 * t81;
t146 = t101 * t119;
t145 = t101 * t123;
t144 = t119 * t123;
t142 = t100 * MDP(24);
t139 = t124 * MDP(11);
t138 = MDP(21) * t144;
t113 = t119 ^ 2;
t137 = t113 * MDP(20) + MDP(17) + 0.2e1 * t138;
t136 = -pkin(5) * t101 - pkin(10) * t100;
t108 = pkin(4) * t120 + pkin(10);
t109 = -t152 * pkin(4) - pkin(5);
t135 = -t100 * t108 + t101 * t109;
t134 = t121 * MDP(12) - t139;
t133 = MDP(22) * t123 - MDP(23) * t119;
t131 = -MDP(25) * t119 - MDP(26) * t123;
t118 = cos(pkin(6));
t116 = sin(pkin(6));
t122 = sin(qJ(2));
t125 = cos(qJ(2));
t89 = (t115 * t125 + t117 * t122) * t116;
t84 = t118 * t124 - t89 * t121;
t85 = t118 * t121 + t124 * t89;
t75 = t120 * t85 - t152 * t84;
t76 = t120 * t84 + t152 * t85;
t130 = -t76 * MDP(19) + (-MDP(18) - t132) * t75;
t129 = -t100 * t132 - t154;
t114 = t123 ^ 2;
t82 = -t120 * t97 + t152 * t98;
t128 = -t81 * MDP(18) - t82 * MDP(19) + ((-t113 + t114) * MDP(21) + MDP(20) * t144 + MDP(15)) * t101 + (-MDP(16) + t143) * t100;
t127 = (t152 * MDP(18) - t120 * MDP(19)) * pkin(4);
t87 = (t115 * t122 - t117 * t125) * t116;
t78 = t100 * pkin(5) - t101 * pkin(10) + t102;
t77 = t81 * t119;
t71 = t119 * t78 + t123 * t82;
t70 = -t119 * t82 + t123 * t78;
t69 = t87 * t119 + t123 * t76;
t68 = -t119 * t76 + t123 * t87;
t1 = [MDP(1) + (t118 ^ 2 + t87 ^ 2 + t89 ^ 2) * MDP(5); t89 * MDP(5) * t151 + (t100 * t68 + t75 * t146) * MDP(25) + (-t69 * t100 + t75 * t145) * MDP(26) + (t125 * MDP(3) - t122 * MDP(4)) * t116 + (-MDP(5) * t150 + t134 + t154) * t87; -0.2e1 * t106 * t139 + t96 * t153 + MDP(2) + (t115 ^ 2 + t117 ^ 2) * MDP(5) * pkin(2) ^ 2 + (MDP(18) * t153 + t142 + 0.2e1 * (-MDP(14) + t133) * t101) * t100 + 0.2e1 * (t100 * t70 + t81 * t146) * MDP(25) + 0.2e1 * (-t71 * t100 + t81 * t145) * MDP(26) + (t114 * MDP(20) + MDP(13) - 0.2e1 * t138) * t101 ^ 2 + (0.2e1 * t106 * MDP(12) + MDP(6) * t121 + 0.2e1 * t124 * MDP(7)) * t121; t118 * MDP(5); 0; MDP(5); t84 * MDP(11) - t85 * MDP(12) + t130; t121 * MDP(8) + t124 * MDP(9) + (t135 * t119 - t147) * MDP(25) + (t135 * t123 + t77) * MDP(26) + (-t121 * MDP(11) - t124 * MDP(12)) * t105 + t128; t129 - t134; -0.2e1 * t109 * t132 + MDP(10) + 0.2e1 * t127 + t137; t130; (t136 * t119 - t147) * MDP(25) + (t136 * t123 + t77) * MDP(26) + t128; t129; t127 + t137 + t132 * (pkin(5) - t109); 0.2e1 * pkin(5) * t132 + t137; t68 * MDP(25) - t69 * MDP(26); t70 * MDP(25) - t71 * MDP(26) + t133 * t101 + t142; t131 * t101; t131 * t108 + t143; t131 * pkin(10) + t143; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
