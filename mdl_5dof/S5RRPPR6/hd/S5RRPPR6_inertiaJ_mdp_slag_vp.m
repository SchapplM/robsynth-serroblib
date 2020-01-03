% Calculate joint inertia matrix for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPPR6_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:06
% EndTime: 2019-12-31 19:33:08
% DurationCPUTime: 0.33s
% Computational Cost: add. (507->102), mult. (961->158), div. (0->0), fcn. (1030->8), ass. (0->58)
t95 = sin(pkin(9));
t97 = cos(pkin(9));
t122 = t95 ^ 2 + t97 ^ 2;
t132 = t122 * MDP(15);
t100 = sin(qJ(2));
t123 = -qJ(3) - pkin(6);
t114 = t123 * t100;
t102 = cos(qJ(2));
t86 = t123 * t102;
t96 = sin(pkin(8));
t98 = cos(pkin(8));
t75 = -t98 * t114 - t96 * t86;
t131 = t75 ^ 2;
t80 = t96 * t100 - t98 * t102;
t130 = 0.2e1 * t80;
t91 = -t98 * pkin(2) - pkin(3);
t85 = -t97 * pkin(4) + t91;
t129 = 0.2e1 * t85;
t128 = 0.2e1 * t102;
t127 = -2 * MDP(18);
t88 = t96 * pkin(2) + qJ(4);
t126 = pkin(7) + t88;
t82 = t98 * t100 + t96 * t102;
t125 = t82 * t95;
t124 = t82 * t97;
t92 = -t102 * pkin(2) - pkin(1);
t73 = t80 * pkin(3) - t82 * qJ(4) + t92;
t77 = t96 * t114 - t98 * t86;
t65 = t95 * t73 + t97 * t77;
t101 = cos(qJ(5));
t99 = sin(qJ(5));
t83 = t101 * t95 + t99 * t97;
t67 = t83 * t82;
t121 = t67 * MDP(20);
t110 = t101 * t97 - t99 * t95;
t68 = t110 * t82;
t120 = t68 * MDP(17);
t119 = t80 * MDP(21);
t118 = t110 * MDP(22);
t117 = t91 * MDP(16);
t116 = t95 * MDP(14);
t115 = t97 * MDP(13);
t64 = t97 * t73 - t95 * t77;
t113 = t122 * MDP(16);
t112 = t64 * t97 + t65 * t95;
t111 = -t64 * t95 + t65 * t97;
t109 = MDP(13) * t95 + MDP(14) * t97;
t62 = t80 * pkin(4) - pkin(7) * t124 + t64;
t63 = -pkin(7) * t125 + t65;
t108 = (t101 * t62 - t99 * t63) * MDP(22) - (t101 * t63 + t99 * t62) * MDP(23);
t107 = t67 * MDP(22) + t68 * MDP(23);
t106 = -t83 * MDP(23) + t118;
t105 = t106 + t115 - t116;
t78 = t126 * t95;
t79 = t126 * t97;
t104 = t83 * MDP(19) + t110 * MDP(20) + (-t101 * t78 - t99 * t79) * MDP(22) - (t101 * t79 - t99 * t78) * MDP(23);
t66 = pkin(4) * t125 + t75;
t1 = [MDP(1) + pkin(1) * MDP(9) * t128 + (t77 ^ 2 + t92 ^ 2 + t131) * MDP(12) + (t64 ^ 2 + t65 ^ 2 + t131) * MDP(16) + (t119 - 0.2e1 * t121) * t80 + (MDP(19) * t130 + t67 * t127 + t120) * t68 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t100 + MDP(5) * t128) * t100 + 0.2e1 * t107 * t66 + (-t77 * MDP(11) + t64 * MDP(13) - t65 * MDP(14) + t108) * t130 + 0.2e1 * (-t112 * MDP(15) + (MDP(11) + t109) * t75) * t82; t100 * MDP(6) + t102 * MDP(7) + (t91 * t125 - t75 * t97) * MDP(13) + (t91 * t124 + t75 * t95) * MDP(14) + t111 * MDP(15) + (t111 * t88 + t75 * t91) * MDP(16) + t83 * t120 + (t110 * t68 - t83 * t67) * MDP(18) + (-t110 * t66 + t85 * t67) * MDP(22) + (t66 * t83 + t85 * t68) * MDP(23) + (-t109 * t88 + t104) * t80 + (-t102 * MDP(10) - t100 * MDP(9)) * pkin(6) + ((-t80 * t96 - t82 * t98) * MDP(11) + (-t75 * t98 + t77 * t96) * MDP(12)) * pkin(2); -t118 * t129 + MDP(8) + (-0.2e1 * t115 + 0.2e1 * t116 + t117) * t91 + (t96 ^ 2 + t98 ^ 2) * MDP(12) * pkin(2) ^ 2 + (MDP(17) * t83 + MDP(23) * t129 - t110 * t127) * t83 + (t113 * t88 + 0.2e1 * t132) * t88; t92 * MDP(12) + t112 * MDP(16) + t105 * t80 - t82 * t132; 0; MDP(12) + t113; t75 * MDP(16) + t109 * t82 + t107; -t105 + t117; 0; MDP(16); t68 * MDP(19) + t108 + t119 - t121; t104; t106; 0; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
