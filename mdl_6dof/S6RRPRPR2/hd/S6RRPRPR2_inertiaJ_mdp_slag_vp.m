% Calculate joint inertia matrix for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR2_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:13:56
% EndTime: 2019-03-09 10:13:58
% DurationCPUTime: 0.43s
% Computational Cost: add. (869->138), mult. (1579->185), div. (0->0), fcn. (1775->8), ass. (0->67)
t171 = MDP(18) - MDP(21);
t170 = -MDP(19) + MDP(22);
t130 = sin(qJ(6));
t133 = cos(qJ(6));
t168 = MDP(29) * t130 + MDP(30) * t133;
t169 = (2 * MDP(22)) + 0.2e1 * t168;
t129 = cos(pkin(10));
t122 = pkin(2) * t129 + pkin(3);
t131 = sin(qJ(4));
t134 = cos(qJ(4));
t128 = sin(pkin(10));
t162 = pkin(2) * t128;
t110 = t122 * t134 - t131 * t162;
t109 = -pkin(4) - t110;
t167 = MDP(23) * t109;
t135 = cos(qJ(2));
t166 = 0.2e1 * t135;
t165 = -2 * MDP(21);
t163 = pkin(4) + pkin(9);
t161 = -qJ(3) - pkin(7);
t160 = MDP(23) * pkin(4);
t132 = sin(qJ(2));
t113 = -t128 * t132 + t129 * t135;
t114 = t128 * t135 + t129 * t132;
t101 = -t134 * t113 + t114 * t131;
t118 = t161 * t132;
t119 = t161 * t135;
t103 = t129 * t118 + t119 * t128;
t93 = -pkin(8) * t114 + t103;
t104 = t128 * t118 - t129 * t119;
t94 = pkin(8) * t113 + t104;
t90 = t131 * t93 + t134 * t94;
t83 = -pkin(5) * t101 + t90;
t79 = t83 * t130;
t80 = t83 * t133;
t159 = qJ(5) * t101;
t140 = t122 * t131 + t134 * t162;
t107 = qJ(5) + t140;
t158 = t101 * t107;
t157 = t130 * t133;
t154 = MDP(29) * t133;
t102 = t113 * t131 + t114 * t134;
t152 = t102 * MDP(27);
t151 = t110 * MDP(18);
t150 = t140 * MDP(19);
t125 = t133 * MDP(26);
t149 = 0.2e1 * t102;
t123 = -pkin(2) * t135 - pkin(1);
t148 = MDP(25) * t157;
t127 = t133 ^ 2;
t147 = t127 * MDP(24) + MDP(17) - 0.2e1 * t148;
t89 = t131 * t94 - t134 * t93;
t146 = t102 * t163 + t159;
t106 = -pkin(9) + t109;
t145 = -t102 * t106 + t158;
t105 = -pkin(3) * t113 + t123;
t144 = MDP(26) * t130 + MDP(27) * t133;
t143 = -t130 * MDP(30) + t154;
t139 = -qJ(5) * t102 + t105;
t126 = t130 ^ 2;
t138 = t80 * MDP(30) + t170 * t90 - t171 * t89 + (t125 + MDP(15)) * t102 + ((-t126 + t127) * MDP(25) + MDP(24) * t157 - MDP(16)) * t101;
t88 = pkin(4) * t101 + t139;
t82 = pkin(5) * t102 + t89;
t81 = t163 * t101 + t139;
t78 = t130 * t82 + t133 * t81;
t77 = -t130 * t81 + t133 * t82;
t1 = [MDP(1) + pkin(1) * MDP(9) * t166 + (t103 ^ 2 + t104 ^ 2 + t123 ^ 2) * MDP(12) + (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) * MDP(23) + (MDP(24) * t126 + 0.2e1 * t148) * t101 ^ 2 + (MDP(19) * t105 - MDP(22) * t88) * t149 + (MDP(13) + MDP(28)) * t102 ^ 2 + (-0.2e1 * MDP(10) * pkin(1) + MDP(4) * t132 + MDP(5) * t166) * t132 + (0.2e1 * t105 * MDP(18) + t88 * t165 + (-MDP(14) + t144) * t149) * t101 + 0.2e1 * (-t103 * t114 + t104 * t113) * MDP(11) + 0.2e1 * (-t101 * t90 + t102 * t89) * MDP(20) + 0.2e1 * (-t101 * t80 + t102 * t77) * MDP(29) + 0.2e1 * (t101 * t79 - t78 * t102) * MDP(30); t132 * MDP(6) + t135 * MDP(7) + (t107 * t90 + t109 * t89) * MDP(23) + ((t113 * t128 - t114 * t129) * MDP(11) + (t103 * t129 + t104 * t128) * MDP(12)) * pkin(2) + (-MDP(10) * t135 - MDP(9) * t132) * pkin(7) + (t102 * t109 - t158) * MDP(20) + (t145 * MDP(30) - t152) * t130 + (-t145 * t133 + t79) * MDP(29) + t138; MDP(8) + (t128 ^ 2 + t129 ^ 2) * MDP(12) * pkin(2) ^ 2 + 0.2e1 * t151 - 0.2e1 * t150 + t147 + ((2 * MDP(21)) + t167) * t109 + (MDP(23) * t107 + t169) * t107; MDP(12) * t123 + MDP(23) * t88 + t171 * t101 + (-t168 - t170) * t102; 0; MDP(12) + MDP(23); (t146 * MDP(30) - t152) * t130 + (-pkin(4) * t102 - t159) * MDP(20) + (-pkin(4) * t89 + qJ(5) * t90) * MDP(23) + (-t146 * t133 + t79) * MDP(29) + t138; t151 - t150 + (-0.2e1 * pkin(4) - t110) * MDP(21) + (0.2e1 * qJ(5) + t140) * MDP(22) + (-pkin(4) * t109 + qJ(5) * t107) * MDP(23) + t147 + t168 * (qJ(5) + t107); 0; (t165 + t160) * pkin(4) + (MDP(23) * qJ(5) + t169) * qJ(5) + t147; MDP(23) * t89 + (MDP(20) + t143) * t102; MDP(21) + t167; 0; MDP(21) - t160; MDP(23); MDP(28) * t102 + MDP(29) * t77 - MDP(30) * t78 + t144 * t101; -MDP(27) * t130 + t143 * t106 + t125; -t168; -t163 * t154 + t125 + (MDP(30) * t163 - MDP(27)) * t130; t143; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
