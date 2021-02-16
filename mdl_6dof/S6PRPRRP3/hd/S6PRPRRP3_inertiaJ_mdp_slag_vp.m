% Calculate joint inertia matrix for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRRP3_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:38:34
% EndTime: 2021-01-16 01:38:37
% DurationCPUTime: 0.56s
% Computational Cost: add. (672->149), mult. (1377->214), div. (0->0), fcn. (1550->10), ass. (0->68)
t109 = cos(qJ(5));
t106 = sin(qJ(5));
t123 = MDP(21) + MDP(23);
t119 = t123 * t106;
t124 = MDP(20) + MDP(22);
t148 = t124 * t109 + MDP(13) - t119;
t102 = sin(pkin(11));
t104 = cos(pkin(11));
t107 = sin(qJ(4));
t143 = cos(qJ(4));
t89 = t143 * t102 + t107 * t104;
t147 = 0.2e1 * t89;
t96 = -t104 * pkin(3) - pkin(2);
t145 = 0.2e1 * t96;
t88 = t107 * t102 - t143 * t104;
t144 = t88 * pkin(5);
t142 = qJ(3) + pkin(8);
t141 = -qJ(6) - pkin(9);
t140 = t102 ^ 2 + t104 ^ 2;
t139 = MDP(25) * pkin(5);
t138 = t106 * t89;
t90 = t142 * t102;
t91 = t142 * t104;
t84 = -t107 * t90 + t143 * t91;
t137 = t109 * t84;
t136 = t109 * t89;
t103 = sin(pkin(6));
t108 = sin(qJ(2));
t135 = t103 * t108;
t110 = cos(qJ(2));
t134 = t103 * t110;
t133 = t104 * MDP(5);
t132 = t88 * MDP(19);
t131 = t89 * MDP(14);
t93 = t141 * t109;
t130 = t93 * MDP(23);
t97 = -t109 * pkin(5) - pkin(4);
t129 = t97 * MDP(25);
t100 = t106 ^ 2;
t101 = t109 ^ 2;
t128 = t100 + t101;
t127 = t106 * MDP(18);
t126 = t106 * MDP(23);
t125 = t109 * MDP(22);
t122 = t106 * t109 * MDP(16);
t82 = t88 * pkin(4) - t89 * pkin(9) + t96;
t74 = -t106 * t84 + t109 * t82;
t121 = -MDP(24) * pkin(5) + MDP(17);
t120 = MDP(22) + t139;
t118 = MDP(20) + t120;
t117 = -pkin(2) * MDP(7) - t133;
t105 = cos(pkin(6));
t85 = -t102 * t135 + t104 * t105;
t86 = t105 * t102 + t104 * t135;
t80 = t107 * t85 + t143 * t86;
t76 = -t106 * t80 - t109 * t134;
t77 = -t106 * t134 + t109 * t80;
t115 = t77 * t106 + t76 * t109;
t114 = t109 * MDP(20) - t106 * MDP(21);
t83 = t107 * t91 + t143 * t90;
t73 = t137 + (-qJ(6) * t89 + t82) * t106;
t113 = t74 * MDP(20) - (t106 * t82 + t137) * MDP(21) - t73 * MDP(23);
t112 = -t125 + t126 + t129;
t92 = t141 * t106;
t79 = t107 * t86 - t143 * t85;
t78 = pkin(5) * t138 + t83;
t72 = -qJ(6) * t136 + t144 + t74;
t1 = [MDP(1) + (t103 ^ 2 * t110 ^ 2 + t85 ^ 2 + t86 ^ 2) * MDP(7) + (t76 ^ 2 + t77 ^ 2 + t79 ^ 2) * MDP(25); (t76 * t72 + t77 * t73 + t79 * t78) * MDP(25) - t115 * MDP(24) * t89 + t123 * (t79 * t136 - t77 * t88) + t124 * (t79 * t138 + t76 * t88) + (-t108 * MDP(4) + (-t88 * MDP(13) + MDP(3) - t117 - t131) * t110) * t103 + (MDP(7) * qJ(3) + MDP(6)) * (-t85 * t102 + t86 * t104); MDP(2) + 0.2e1 * pkin(2) * t133 + (t140 * qJ(3) ^ 2 + pkin(2) ^ 2) * MDP(7) + t131 * t145 + (t72 ^ 2 + t73 ^ 2 + t78 ^ 2) * MDP(25) + (t101 * MDP(15) + MDP(8) - 0.2e1 * t122) * t89 ^ 2 + (MDP(13) * t145 + t132 + (t109 * MDP(17) - MDP(9) - t127) * t147) * t88 + 0.2e1 * (t72 * MDP(22) + t113) * t88 + ((t83 * MDP(21) + t78 * MDP(23) - t72 * MDP(24)) * t109 + (t83 * MDP(20) + t78 * MDP(22) - t73 * MDP(24)) * t106) * t147 + 0.2e1 * t140 * MDP(6) * qJ(3); t115 * MDP(25) - MDP(7) * t134; (t73 * t106 + t72 * t109) * MDP(25) + (-t128 * MDP(24) + MDP(14)) * t89 + t148 * t88 + t117; t128 * MDP(25) + MDP(7); -t80 * MDP(14) + (-t76 * t106 + t77 * t109) * MDP(24) + (t76 * t92 - t77 * t93) * MDP(25) + (t129 - t148) * t79; -t84 * MDP(14) + (-t72 * t106 + t73 * t109) * MDP(24) + (t72 * t92 - t73 * t93) * MDP(25) + t112 * t78 + (-MDP(13) - t114) * t83 + (t106 * MDP(17) + t109 * MDP(18) + t92 * MDP(22) + t130 - MDP(11) + (-MDP(20) * t106 - MDP(21) * t109) * pkin(9)) * t88 + (MDP(10) + (-t100 + t101) * MDP(16) + (-pkin(4) * MDP(21) + t97 * MDP(23) - t92 * MDP(24)) * t109 + (t109 * MDP(15) - pkin(4) * MDP(20) + t97 * MDP(22) + t93 * MDP(24)) * t106) * t89; (-t106 * t93 + t109 * t92) * MDP(25); MDP(12) + t100 * MDP(15) + 0.2e1 * t122 + 0.2e1 * (-t92 * t106 - t93 * t109) * MDP(24) + (t92 ^ 2 + t93 ^ 2) * MDP(25) + (-0.2e1 * t125 + 0.2e1 * t126 + t129) * t97 + 0.2e1 * t114 * pkin(4); t118 * t76 - t123 * t77; t132 + (t74 + 0.2e1 * t144) * MDP(22) + t72 * t139 + (-t127 + (-MDP(22) * qJ(6) + t121) * t109) * t89 + t113; t118 * t109 - t119; t130 + t120 * t92 + (-MDP(21) * pkin(9) + MDP(18)) * t109 + (-MDP(20) * pkin(9) + t121) * t106; MDP(19) + (0.2e1 * MDP(22) + t139) * pkin(5); t79 * MDP(25); t78 * MDP(25) + (t106 * MDP(22) + t109 * MDP(23)) * t89; 0; t112; 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
