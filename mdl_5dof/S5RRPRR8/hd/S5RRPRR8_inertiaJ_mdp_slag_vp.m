% Calculate joint inertia matrix for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:36
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR8_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:35:14
% EndTime: 2021-01-15 21:35:16
% DurationCPUTime: 0.32s
% Computational Cost: add. (519->101), mult. (988->147), div. (0->0), fcn. (1095->8), ass. (0->56)
t100 = sin(qJ(5));
t103 = cos(qJ(5));
t126 = t100 * MDP(24) + t103 * MDP(25);
t109 = t103 * MDP(27) - t100 * MDP(28);
t102 = sin(qJ(2));
t105 = cos(qJ(2));
t98 = sin(pkin(9));
t99 = cos(pkin(9));
t85 = t98 * t102 - t99 * t105;
t92 = -t105 * pkin(2) - pkin(1);
t79 = t85 * pkin(3) + t92;
t131 = 0.2e1 * t79;
t130 = 0.2e1 * t105;
t129 = pkin(2) * t98;
t127 = -qJ(3) - pkin(6);
t101 = sin(qJ(4));
t104 = cos(qJ(4));
t88 = t127 * t102;
t89 = t127 * t105;
t77 = t99 * t88 + t98 * t89;
t86 = t99 * t102 + t98 * t105;
t67 = -t86 * pkin(7) + t77;
t78 = t98 * t88 - t99 * t89;
t68 = -t85 * pkin(7) + t78;
t64 = t101 * t68 - t104 * t67;
t60 = t64 * t100;
t125 = t64 * t103;
t124 = t100 * t103;
t75 = t101 * t86 + t104 * t85;
t123 = t75 * MDP(26);
t76 = -t101 * t85 + t104 * t86;
t122 = t76 * MDP(21);
t91 = t99 * pkin(2) + pkin(3);
t82 = -t101 * t129 + t104 * t91;
t121 = t82 * MDP(20);
t83 = -t101 * t91 - t104 * t129;
t120 = t83 * MDP(21);
t119 = t85 * MDP(11);
t118 = t86 * MDP(12);
t117 = t92 * MDP(14);
t113 = MDP(23) * t124;
t96 = t100 ^ 2;
t114 = t96 * MDP(22) + MDP(19) + 0.2e1 * t113;
t112 = -pkin(4) * t76 - pkin(8) * t75;
t80 = -pkin(4) - t82;
t81 = pkin(8) - t83;
t111 = -t75 * t81 + t76 * t80;
t110 = MDP(24) * t103 - MDP(25) * t100;
t108 = -MDP(27) * t100 - MDP(28) * t103;
t65 = t101 * t67 + t104 * t68;
t97 = t103 ^ 2;
t107 = -t64 * MDP(20) - t65 * MDP(21) + ((-t96 + t97) * MDP(23) + MDP(22) * t124 + MDP(17)) * t76 + (-MDP(18) + t126) * t75;
t63 = t75 * pkin(4) - t76 * pkin(8) + t79;
t59 = t100 * t63 + t103 * t65;
t58 = -t100 * t65 + t103 * t63;
t1 = [MDP(1) + pkin(1) * MDP(9) * t130 + (t77 ^ 2 + t78 ^ 2) * MDP(14) + t122 * t131 + (t117 + 0.2e1 * t118 + 0.2e1 * t119) * t92 + (t97 * MDP(22) + MDP(15) - 0.2e1 * t113) * t76 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t102 + MDP(5) * t130) * t102 + (MDP(20) * t131 + t123 + 0.2e1 * (-MDP(16) + t110) * t76) * t75 + 0.2e1 * (-t77 * t86 - t78 * t85) * MDP(13) + 0.2e1 * (t58 * t75 + t76 * t60) * MDP(27) + 0.2e1 * (t76 * t125 - t59 * t75) * MDP(28); (t111 * t100 - t125) * MDP(27) + (t111 * t103 + t60) * MDP(28) + ((-t85 * t98 - t86 * t99) * MDP(13) + (t77 * t99 + t78 * t98) * MDP(14)) * pkin(2) + t102 * MDP(6) + t105 * MDP(7) + t77 * MDP(11) - t78 * MDP(12) + (-t105 * MDP(10) - t102 * MDP(9)) * pkin(6) + t107; MDP(8) + (t98 ^ 2 + t99 ^ 2) * MDP(14) * pkin(2) ^ 2 + 0.2e1 * t121 + 0.2e1 * t120 + t114 - 0.2e1 * t109 * t80 + 0.2e1 * (t99 * MDP(11) - t98 * MDP(12)) * pkin(2); t119 + t118 + t117 + t122 + (MDP(20) + t109) * t75; 0; MDP(14); (t112 * t100 - t125) * MDP(27) + (t112 * t103 + t60) * MDP(28) + t107; t114 + t120 + t121 + t109 * (pkin(4) - t80); 0; 0.2e1 * pkin(4) * t109 + t114; t58 * MDP(27) - t59 * MDP(28) + t110 * t76 + t123; t108 * t81 + t126; t109; t108 * pkin(8) + t126; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
