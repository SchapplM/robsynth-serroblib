% Calculate joint inertia matrix for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRPP6_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:02:00
% EndTime: 2019-12-31 21:02:01
% DurationCPUTime: 0.43s
% Computational Cost: add. (604->143), mult. (1137->208), div. (0->0), fcn. (1076->6), ass. (0->54)
t105 = sin(qJ(3));
t107 = cos(qJ(3));
t103 = sin(pkin(8));
t104 = cos(pkin(8));
t125 = -qJ(4) - pkin(7);
t113 = t125 * t105;
t91 = t125 * t107;
t78 = -t103 * t91 - t104 * t113;
t80 = t103 * t113 - t104 * t91;
t133 = t105 * MDP(13) + t107 * MDP(14) - t78 * MDP(20) + t80 * MDP(22) - (t105 * MDP(16) + t107 * MDP(17)) * pkin(7);
t132 = 2 * MDP(18);
t131 = 0.2e1 * MDP(20);
t130 = 2 * MDP(21);
t129 = 0.2e1 * MDP(22);
t128 = pkin(6) * t105;
t127 = pkin(6) * t107;
t108 = cos(qJ(2));
t126 = pkin(6) * t108;
t106 = sin(qJ(2));
t122 = t106 * t107;
t90 = -t108 * pkin(2) - t106 * pkin(7) - pkin(1);
t87 = t107 * t90;
t74 = -qJ(4) * t122 + t87 + (-pkin(3) - t128) * t108;
t117 = t107 * t126;
t76 = t117 + (-qJ(4) * t106 + t90) * t105;
t67 = t103 * t74 + t104 * t76;
t124 = t105 * t106;
t89 = pkin(3) * t124 + t106 * pkin(6);
t123 = t105 * t107;
t85 = t103 * t105 - t104 * t107;
t86 = t103 * t107 + t104 * t105;
t98 = -t107 * pkin(3) - pkin(2);
t71 = t85 * pkin(4) - t86 * qJ(5) + t98;
t121 = t71 * MDP(23);
t120 = t85 * MDP(20);
t119 = t86 * MDP(22);
t118 = t78 ^ 2 + t80 ^ 2;
t116 = MDP(12) * t123;
t83 = t86 * t106;
t84 = -t103 * t124 + t104 * t122;
t115 = t78 * t84 - t80 * t83;
t66 = -t103 * t76 + t104 * t74;
t112 = t107 * MDP(13) - t105 * MDP(14);
t102 = t107 ^ 2;
t101 = t106 ^ 2;
t100 = t105 ^ 2;
t96 = t104 * pkin(3) + pkin(4);
t94 = t103 * pkin(3) + qJ(5);
t82 = t105 * t90 + t117;
t81 = -t105 * t126 + t87;
t68 = t83 * pkin(4) - t84 * qJ(5) + t89;
t65 = t108 * pkin(4) - t66;
t64 = -t108 * qJ(5) + t67;
t1 = [MDP(1) - 0.2e1 * pkin(1) * t106 * MDP(10) + (t66 ^ 2 + t67 ^ 2 + t89 ^ 2) * MDP(19) + (t64 ^ 2 + t65 ^ 2 + t68 ^ 2) * MDP(23) + (t102 * MDP(11) + MDP(4) - 0.2e1 * t116) * t101 + (t108 * MDP(15) + 0.2e1 * pkin(1) * MDP(9) + 0.2e1 * (MDP(5) - t112) * t106) * t108 + 0.2e1 * (t101 * t128 - t81 * t108) * MDP(16) + 0.2e1 * (t101 * t127 + t82 * t108) * MDP(17) + (-t66 * t84 - t67 * t83) * t132 + (t65 * t108 + t68 * t83) * t131 + (-t64 * t83 + t65 * t84) * t130 + (-t64 * t108 - t68 * t84) * t129; (-t66 * t86 - t67 * t85 + t115) * MDP(18) + (-t66 * t78 + t67 * t80 + t89 * t98) * MDP(19) + (t68 * t85 + t71 * t83) * MDP(20) + (-t64 * t85 + t65 * t86 + t115) * MDP(21) + (-t68 * t86 - t71 * t84) * MDP(22) + (t64 * t80 + t65 * t78 + t68 * t71) * MDP(23) + (-pkin(6) * MDP(10) + MDP(7) - t133) * t108 + (MDP(6) - pkin(6) * MDP(9) + MDP(11) * t123 + (-t100 + t102) * MDP(12) + (-pkin(2) * t105 - t127) * MDP(16) + (-pkin(2) * t107 + t128) * MDP(17)) * t106; MDP(8) + t100 * MDP(11) + 0.2e1 * t116 + (t98 ^ 2 + t118) * MDP(19) + t118 * MDP(23) + (-0.2e1 * t119 + 0.2e1 * t120 + t121) * t71 + 0.2e1 * (t107 * MDP(16) - t105 * MDP(17)) * pkin(2) + (t132 + t130) * (t78 * t86 - t80 * t85); t81 * MDP(16) - t82 * MDP(17) + t66 * MDP(20) + (-t94 * t83 - t96 * t84) * MDP(21) + t67 * MDP(22) + (t64 * t94 - t65 * t96) * MDP(23) + t112 * t106 + (-MDP(15) + (-pkin(4) - t96) * MDP(20) + (-qJ(5) - t94) * MDP(22)) * t108 + ((-t103 * t83 - t104 * t84) * MDP(18) + (t103 * t67 + t104 * t66) * MDP(19)) * pkin(3); (-t94 * t85 - t96 * t86) * MDP(21) + (-t78 * t96 + t80 * t94) * MDP(23) + ((-t103 * t85 - t104 * t86) * MDP(18) + (t103 * t80 - t104 * t78) * MDP(19)) * pkin(3) + t133; MDP(15) + (t94 ^ 2 + t96 ^ 2) * MDP(23) + (t103 ^ 2 + t104 ^ 2) * MDP(19) * pkin(3) ^ 2 + t96 * t131 + t94 * t129; t89 * MDP(19) + t83 * MDP(20) - t84 * MDP(22) + t68 * MDP(23); t98 * MDP(19) - t119 + t120 + t121; 0; MDP(19) + MDP(23); t108 * MDP(20) + t84 * MDP(21) + t65 * MDP(23); t86 * MDP(21) + t78 * MDP(23); -t96 * MDP(23) - MDP(20); 0; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
