% Calculate joint inertia matrix for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRRR12_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:21
% EndTime: 2024-09-28 18:08:21
% DurationCPUTime: 0.25s
% Computational Cost: add. (354->92), mult. (881->126), div. (0->0), fcn. (885->12), ass. (0->62)
t111 = cos(pkin(6));
t146 = t111 * MDP(15);
t109 = sin(pkin(6));
t117 = cos(qJ(5));
t133 = t109 * t117;
t113 = sin(qJ(5));
t134 = t109 * t113;
t145 = MDP(13) * t134 + MDP(14) * t133;
t108 = t109 ^ 2;
t136 = t108 * t113;
t143 = 2 * MDP(16);
t142 = 2 * MDP(17);
t115 = sin(qJ(3));
t141 = pkin(2) * t115;
t135 = t108 * t117;
t131 = t111 * t117;
t106 = t109 * pkin(10);
t119 = cos(qJ(3));
t103 = pkin(2) * t119 + pkin(3);
t114 = sin(qJ(4));
t118 = cos(qJ(4));
t89 = -t103 * t114 - t118 * t141;
t81 = t106 - t89;
t97 = t118 * t103;
t88 = -t114 * t141 + t97;
t87 = pkin(4) + t88;
t71 = -t113 * t81 + t87 * t131;
t140 = t71 * t111 + t87 * t135;
t107 = t118 * pkin(3);
t102 = t107 + pkin(4);
t96 = pkin(3) * t114 + t106;
t78 = t102 * t131 - t113 * t96;
t139 = t102 * t135 + t78 * t111;
t138 = t88 * MDP(9);
t90 = pkin(4) * t131 - pkin(10) * t134;
t137 = pkin(4) * t135 + t90 * t111;
t132 = t111 * t113;
t130 = t119 * MDP(6);
t129 = t89 * MDP(10);
t128 = MDP(10) * t114;
t127 = t145 + t146;
t112 = cos(pkin(5));
t110 = sin(pkin(5));
t116 = sin(qJ(2));
t120 = cos(qJ(2));
t85 = (-t115 * t116 + t119 * t120) * t110;
t86 = (t115 * t120 + t116 * t119) * t110;
t75 = -t114 * t86 + t118 * t85;
t125 = t109 * t112 + t111 * t75;
t76 = t114 * t85 + t118 * t86;
t67 = -t113 * t76 + t125 * t117;
t68 = t125 * t113 + t117 * t76;
t69 = -t109 * t75 + t111 * t112;
t126 = t75 * MDP(9) - t76 * MDP(10) + (t111 * t67 - t69 * t133) * MDP(16) + (-t111 * t68 + t69 * t134) * MDP(17);
t124 = t85 * MDP(6) - t86 * MDP(7) + t126;
t123 = MDP(8) + (MDP(11) * t136 + 0.2e1 * MDP(12) * t135) * t113 + (0.2e1 * t145 + t146) * t111;
t122 = MDP(5) + t123;
t121 = (MDP(9) * t118 - t128) * pkin(3);
t91 = pkin(4) * t132 + pkin(10) * t133;
t79 = t102 * t132 + t117 * t96;
t72 = t117 * t81 + t87 * t132;
t1 = [MDP(1); (MDP(3) * t120 - MDP(4) * t116) * t110 + t124; MDP(2) + 0.2e1 * (-MDP(7) * t115 + t130) * pkin(2) + 0.2e1 * t138 + 0.2e1 * t129 + t140 * t143 + (-t111 * t72 - t87 * t136) * t142 + t122; t124; (t107 + t97) * MDP(9) + (t139 + t140) * MDP(16) + ((-t72 - t79) * t111 + (-t102 - t87) * t136) * MDP(17) + (-pkin(3) - t103) * t128 + (t130 + (-MDP(10) * t118 - MDP(9) * t114 - MDP(7)) * t115) * pkin(2) + t122; 0.2e1 * t121 + t139 * t143 + (-t102 * t136 - t111 * t79) * t142 + t122; t126; t138 + t129 + (t137 + t140) * MDP(16) + ((-t72 - t91) * t111 + (-pkin(4) - t87) * t136) * MDP(17) + t123; (t137 + t139) * MDP(16) + ((-t79 - t91) * t111 + (-pkin(4) - t102) * t136) * MDP(17) + t121 + t123; t137 * t143 + (-pkin(4) * t136 - t111 * t91) * t142 + t123; MDP(16) * t67 - MDP(17) * t68; MDP(16) * t71 - MDP(17) * t72 + t127; MDP(16) * t78 - MDP(17) * t79 + t127; MDP(16) * t90 - MDP(17) * t91 + t127; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
