% Calculate joint inertia matrix for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR6_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:52
% EndTime: 2019-12-31 17:30:54
% DurationCPUTime: 0.37s
% Computational Cost: add. (350->123), mult. (862->185), div. (0->0), fcn. (864->8), ass. (0->67)
t86 = sin(pkin(4));
t129 = 0.2e1 * t86;
t88 = sin(qJ(4));
t91 = cos(qJ(4));
t95 = -(MDP(23) * t88 + MDP(24) * t91) * pkin(8) + t88 * MDP(20) + t91 * MDP(21);
t128 = 0.2e1 * MDP(23);
t127 = 0.2e1 * MDP(24);
t90 = sin(qJ(2));
t126 = pkin(1) * t90;
t93 = cos(qJ(2));
t125 = pkin(1) * t93;
t124 = pkin(7) * t88;
t123 = pkin(7) * t91;
t92 = cos(qJ(3));
t122 = pkin(7) * t92;
t118 = t86 * t93;
t104 = pkin(6) * t118;
t87 = cos(pkin(4));
t73 = t104 + (pkin(7) + t126) * t87;
t74 = (-pkin(2) * t93 - pkin(7) * t90 - pkin(1)) * t86;
t89 = sin(qJ(3));
t66 = -t89 * t73 + t92 * t74;
t64 = pkin(3) * t118 - t66;
t121 = t64 * t88;
t120 = t64 * t91;
t119 = t86 * t90;
t117 = pkin(2) * MDP(16);
t116 = pkin(2) * MDP(17);
t115 = t87 * MDP(8);
t114 = t90 * MDP(6);
t113 = MDP(15) * t93;
t76 = t92 * t119 + t87 * t89;
t68 = t91 * t118 + t76 * t88;
t112 = t68 * MDP(21);
t69 = -t88 * t118 + t76 * t91;
t111 = t69 * MDP(18);
t110 = t69 * MDP(20);
t75 = t89 * t119 - t87 * t92;
t109 = t75 * MDP(22);
t108 = t76 * MDP(12);
t107 = t76 * MDP(13);
t106 = t91 * MDP(18);
t105 = t92 * MDP(22);
t103 = t88 * t91 * MDP(19);
t102 = pkin(7) * MDP(16) - MDP(13);
t101 = pkin(7) * MDP(17) - MDP(14);
t67 = t92 * t73 + t89 * t74;
t100 = MDP(20) * t91 - MDP(21) * t88;
t79 = -t92 * pkin(3) - t89 * pkin(8) - pkin(2);
t70 = -t88 * t122 + t91 * t79;
t71 = t91 * t122 + t88 * t79;
t98 = t70 * MDP(23) - t71 * MDP(24);
t80 = pkin(6) * t119;
t72 = t80 + (-pkin(2) - t125) * t87;
t96 = -MDP(12) + t100;
t63 = t75 * pkin(3) - t76 * pkin(8) + t72;
t65 = -pkin(8) * t118 + t67;
t61 = t91 * t63 - t88 * t65;
t62 = t88 * t63 + t91 * t65;
t94 = t61 * MDP(23) - t62 * MDP(24) + t109 + t110 - t112;
t85 = t91 ^ 2;
t84 = t89 ^ 2;
t83 = t88 ^ 2;
t82 = t86 ^ 2;
t78 = t87 * t126 + t104;
t77 = t87 * t125 - t80;
t1 = [t82 * t90 ^ 2 * MDP(4) + t76 ^ 2 * MDP(11) + MDP(1) + (t114 * t129 + t115) * t87 + (-0.2e1 * t68 * MDP(19) + t111) * t69 + ((MDP(7) * t87 - t107) * t129 + (0.2e1 * MDP(5) * t90 + t113) * t82) * t93 + (0.2e1 * MDP(14) * t118 - 0.2e1 * t108 + t109 + 0.2e1 * t110 - 0.2e1 * t112) * t75 + 0.2e1 * (t82 * t125 + t77 * t87) * MDP(9) + 0.2e1 * (-t82 * t126 - t78 * t87) * MDP(10) + 0.2e1 * (-t66 * t118 + t72 * t75) * MDP(16) + 0.2e1 * (t67 * t118 + t72 * t76) * MDP(17) + (t61 * t75 + t64 * t68) * t128 + (-t62 * t75 + t64 * t69) * t127; -t76 * t116 - t78 * MDP(10) + t115 + t77 * MDP(9) + (t93 * MDP(7) + t114) * t86 + (t98 - t117) * t75 + (-t72 * MDP(16) + t101 * t118 + t108 - t94) * t92 + (t76 * MDP(11) + t72 * MDP(17) + t69 * t106 + (-t68 * t91 - t69 * t88) * MDP(19) + (pkin(7) * t68 + t121) * MDP(23) + (pkin(7) * t69 + t120) * MDP(24) + t102 * t118 + t96 * t75) * t89; MDP(8) + (t105 + 0.2e1 * t117) * t92 + (t85 * MDP(18) + MDP(11) - 0.2e1 * t103) * t84 + (t84 * t124 - t70 * t92) * t128 + (t84 * t123 + t71 * t92) * t127 + 0.2e1 * (-t92 * t96 - t116) * t89; t107 - t86 * t113 + t66 * MDP(16) - t67 * MDP(17) + t88 * t111 + (-t88 * t68 + t69 * t91) * MDP(19) + (-pkin(3) * t68 - t120) * MDP(23) + (-pkin(3) * t69 + t121) * MDP(24) + (-MDP(14) + t95) * t75; (-t101 - t95) * t92 + (t88 * t106 + (-t83 + t85) * MDP(19) + (-pkin(3) * t88 - t123) * MDP(23) + (-pkin(3) * t91 + t124) * MDP(24) - t102) * t89; 0.2e1 * t103 + t83 * MDP(18) + MDP(15) + 0.2e1 * (t91 * MDP(23) - t88 * MDP(24)) * pkin(3); t94; t100 * t89 - t105 + t98; t95; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
