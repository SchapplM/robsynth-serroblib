% Calculate joint inertia matrix for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR9_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:08:01
% EndTime: 2019-12-31 19:08:02
% DurationCPUTime: 0.30s
% Computational Cost: add. (464->88), mult. (889->124), div. (0->0), fcn. (1010->8), ass. (0->52)
t96 = sin(qJ(5));
t99 = cos(qJ(5));
t121 = t96 * MDP(24) + t99 * MDP(25);
t105 = t99 * MDP(27) - t96 * MDP(28);
t100 = cos(qJ(3));
t94 = sin(pkin(9));
t95 = cos(pkin(9));
t98 = sin(qJ(3));
t79 = -t100 * t95 + t98 * t94;
t84 = -t95 * pkin(2) - pkin(1);
t76 = t79 * pkin(3) + t84;
t129 = 0.2e1 * t76;
t128 = 0.2e1 * t84;
t126 = cos(qJ(4));
t125 = pkin(1) * MDP(7);
t122 = pkin(6) + qJ(2);
t81 = t122 * t94;
t82 = t122 * t95;
t110 = -t100 * t81 - t98 * t82;
t80 = t100 * t94 + t98 * t95;
t66 = -t80 * pkin(7) + t110;
t107 = -t100 * t82 + t98 * t81;
t67 = -t79 * pkin(7) - t107;
t97 = sin(qJ(4));
t62 = -t126 * t66 + t97 * t67;
t59 = t62 * t96;
t124 = t62 * t99;
t123 = t96 * t99;
t119 = t94 * MDP(5);
t118 = t95 * MDP(4);
t74 = t126 * t79 + t97 * t80;
t117 = t74 * MDP(26);
t75 = t126 * t80 - t97 * t79;
t116 = t75 * MDP(21);
t115 = t79 * MDP(13);
t112 = MDP(23) * t123;
t92 = t96 ^ 2;
t111 = t92 * MDP(22) + MDP(19) + 0.2e1 * t112;
t109 = -pkin(4) * t75 - pkin(8) * t74;
t85 = t97 * pkin(3) + pkin(8);
t86 = -t126 * pkin(3) - pkin(4);
t108 = -t74 * t85 + t75 * t86;
t106 = MDP(24) * t99 - MDP(25) * t96;
t104 = -MDP(27) * t96 - MDP(28) * t99;
t63 = t126 * t67 + t97 * t66;
t93 = t99 ^ 2;
t103 = -t62 * MDP(20) - t63 * MDP(21) + ((-t92 + t93) * MDP(23) + MDP(22) * t123 + MDP(17)) * t75 + (-MDP(18) + t121) * t74;
t102 = (t126 * MDP(20) - t97 * MDP(21)) * pkin(3);
t64 = t74 * pkin(4) - t75 * pkin(8) + t76;
t58 = t99 * t63 + t96 * t64;
t57 = -t96 * t63 + t99 * t64;
t1 = [t115 * t128 + t116 * t129 + MDP(1) + (0.2e1 * t118 - 0.2e1 * t119 + t125) * pkin(1) + (MDP(14) * t128 + MDP(8) * t80 - 0.2e1 * t79 * MDP(9)) * t80 + (t93 * MDP(22) + MDP(15) - 0.2e1 * t112) * t75 ^ 2 + (MDP(20) * t129 + t117 + 0.2e1 * (-MDP(16) + t106) * t75) * t74 + 0.2e1 * (t57 * t74 + t75 * t59) * MDP(27) + 0.2e1 * (t75 * t124 - t58 * t74) * MDP(28) + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t94 ^ 2 + t95 ^ 2) * qJ(2); t115 + t80 * MDP(14) + t116 - t118 + t119 - t125 + (MDP(20) + t105) * t74; MDP(7); t80 * MDP(10) - t79 * MDP(11) + t110 * MDP(13) + t107 * MDP(14) + (t108 * t96 - t124) * MDP(27) + (t108 * t99 + t59) * MDP(28) + t103; 0; -0.2e1 * t105 * t86 + MDP(12) + 0.2e1 * t102 + t111; (t109 * t96 - t124) * MDP(27) + (t109 * t99 + t59) * MDP(28) + t103; 0; t102 + t111 + t105 * (pkin(4) - t86); 0.2e1 * pkin(4) * t105 + t111; t57 * MDP(27) - t58 * MDP(28) + t106 * t75 + t117; t105; t104 * t85 + t121; t104 * pkin(8) + t121; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
