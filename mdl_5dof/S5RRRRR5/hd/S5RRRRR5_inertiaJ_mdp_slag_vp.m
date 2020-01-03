% Calculate joint inertia matrix for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRRR5_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:24
% EndTime: 2020-01-03 12:13:25
% DurationCPUTime: 0.21s
% Computational Cost: add. (317->98), mult. (536->114), div. (0->0), fcn. (478->8), ass. (0->57)
t110 = sin(qJ(4));
t114 = cos(qJ(4));
t126 = MDP(15) * t114;
t117 = -0.2e1 * MDP(16) * t110 + 0.2e1 * t126;
t109 = sin(qJ(5));
t113 = cos(qJ(5));
t116 = cos(qJ(2));
t102 = t116 * pkin(1) + pkin(2);
t111 = sin(qJ(3));
t115 = cos(qJ(3));
t112 = sin(qJ(2));
t135 = pkin(1) * t112;
t84 = -t111 * t102 - t115 * t135;
t82 = pkin(8) - t84;
t73 = (-pkin(9) - t82) * t110;
t107 = t114 * pkin(9);
t74 = t114 * t82 + t107;
t138 = (-t109 * t74 + t113 * t73) * MDP(22) + (-t109 * t73 - t113 * t74) * MDP(23);
t100 = t111 * pkin(2) + pkin(8);
t88 = (-pkin(9) - t100) * t110;
t89 = t114 * t100 + t107;
t137 = (-t109 * t89 + t113 * t88) * MDP(22) + (-t109 * t88 - t113 * t89) * MDP(23);
t93 = (-pkin(8) - pkin(9)) * t110;
t94 = t114 * pkin(8) + t107;
t136 = (-t109 * t94 + t113 * t93) * MDP(22) + (-t109 * t93 - t113 * t94) * MDP(23);
t134 = pkin(3) * t110;
t133 = t115 * pkin(2);
t90 = t109 * t110 - t113 * t114;
t91 = t109 * t114 + t113 * t110;
t132 = t91 * MDP(19) - t90 * MDP(20);
t131 = -t115 * t102 + t111 * t135;
t130 = t131 * MDP(8);
t129 = t84 * MDP(9);
t128 = t111 * MDP(9);
t127 = t116 * MDP(5);
t103 = -t114 * pkin(4) - pkin(3);
t125 = t110 * MDP(12) + t114 * MDP(13) + t132;
t124 = MDP(7) + (MDP(17) * t91 - 0.2e1 * MDP(18) * t90) * t91 + (MDP(10) * t110 + 0.2e1 * MDP(11) * t114) * t110;
t122 = -MDP(15) * t110 - MDP(16) * t114;
t121 = MDP(4) + t124;
t120 = (t115 * MDP(8) - t128) * pkin(2);
t119 = (MDP(22) * t113 - MDP(23) * t109) * pkin(4);
t118 = 0.2e1 * MDP(22) * t90 + 0.2e1 * MDP(23) * t91;
t108 = pkin(3) * t114;
t101 = -pkin(3) - t133;
t96 = t101 * t110;
t92 = t103 - t133;
t81 = -pkin(3) + t131;
t80 = t103 * t91;
t79 = t103 * t90;
t78 = t81 * t110;
t77 = t103 + t131;
t76 = t92 * t91;
t75 = t92 * t90;
t69 = t77 * t91;
t68 = t77 * t90;
t1 = [MDP(1) - t81 * t117 + t77 * t118 + 0.2e1 * (-t112 * MDP(6) + t127) * pkin(1) - 0.2e1 * t130 + 0.2e1 * t129 + t121; (-t131 + t133) * MDP(8) + (t96 + t78) * MDP(16) + (t75 + t68) * MDP(22) + (t76 + t69) * MDP(23) + (-t101 - t81) * t126 + (-pkin(2) - t102) * t128 + (t127 + (-MDP(9) * t115 - MDP(6)) * t112) * pkin(1) + t121; -t101 * t117 + t92 * t118 + 0.2e1 * t120 + t121; -t130 + t129 + (-t81 * t114 + t108) * MDP(15) + (t78 - t134) * MDP(16) + (t79 + t68) * MDP(22) + (t80 + t69) * MDP(23) + t124; (-t101 * t114 + t108) * MDP(15) + (t96 - t134) * MDP(16) + (t79 + t75) * MDP(22) + (t80 + t76) * MDP(23) + t120 + t124; pkin(3) * t117 + t103 * t118 + t124; t122 * t82 + t125 + t138; t122 * t100 + t125 + t137; t122 * pkin(8) + t125 + t136; MDP(14) + MDP(21) + 0.2e1 * t119; t132 + t138; t132 + t137; t132 + t136; MDP(21) + t119; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
