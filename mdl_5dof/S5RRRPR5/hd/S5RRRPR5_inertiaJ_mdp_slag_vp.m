% Calculate joint inertia matrix for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:12
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR5_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:10:46
% EndTime: 2021-01-15 23:10:49
% DurationCPUTime: 0.42s
% Computational Cost: add. (637->116), mult. (1194->176), div. (0->0), fcn. (1299->8), ass. (0->57)
t108 = sin(qJ(5));
t111 = cos(qJ(5));
t135 = t108 * MDP(24) + t111 * MDP(25);
t121 = t111 * MDP(27) - t108 * MDP(28);
t109 = sin(qJ(3));
t110 = sin(qJ(2));
t112 = cos(qJ(3));
t113 = cos(qJ(2));
t123 = t109 * t110 - t112 * t113;
t100 = -t113 * pkin(2) - pkin(1);
t146 = 0.2e1 * t100;
t145 = 2 * MDP(18);
t144 = pkin(6) + pkin(7);
t143 = pkin(2) * t109;
t106 = sin(pkin(9));
t107 = cos(pkin(9));
t128 = t144 * t113;
t129 = t144 * t110;
t116 = -t109 * t128 - t112 * t129;
t94 = t109 * t113 + t110 * t112;
t115 = -t94 * qJ(4) + t116;
t78 = t123 * (-qJ(4) - t144);
t72 = t106 * t78 - t107 * t115;
t141 = t111 * t72;
t68 = t72 * t108;
t82 = t106 * t94 + t107 * t123;
t140 = MDP(26) * t82;
t139 = t108 * t111;
t83 = -t106 * t123 + t107 * t94;
t136 = t83 * MDP(19);
t134 = MDP(16) * t112;
t133 = t106 * MDP(19);
t104 = t108 ^ 2;
t126 = MDP(23) * t139;
t127 = t104 * MDP(22) + MDP(15) + 0.2e1 * t126;
t99 = pkin(2) * t112 + pkin(3);
t95 = t107 * t99;
t89 = -t106 * t143 + t95;
t87 = -pkin(4) - t89;
t90 = t106 * t99 + t107 * t143;
t88 = pkin(8) + t90;
t125 = -t82 * t88 + t83 * t87;
t97 = pkin(3) * t106 + pkin(8);
t98 = -pkin(3) * t107 - pkin(4);
t124 = -t82 * t97 + t83 * t98;
t122 = t107 * MDP(18) - t133;
t120 = -MDP(27) * t108 - MDP(28) * t111;
t119 = (MDP(24) * t111 - MDP(25) * t108) * t83;
t118 = 0.2e1 * t121;
t105 = t111 ^ 2;
t74 = t106 * t115 + t107 * t78;
t117 = -t72 * MDP(18) - t74 * MDP(19) + t116 * MDP(16) + (t109 * t129 - t112 * t128) * MDP(17) - t123 * MDP(14) + t94 * MDP(13) + ((-t104 + t105) * MDP(23) + MDP(22) * t139) * t83 + t135 * t82;
t86 = t123 * pkin(3) + t100;
t71 = t82 * pkin(4) - t83 * pkin(8) + t86;
t67 = t108 * t71 + t111 * t74;
t66 = -t108 * t74 + t111 * t71;
t1 = [MDP(1) + t123 * MDP(16) * t146 + 0.2e1 * t86 * t136 + (t72 ^ 2 + t74 ^ 2 + t86 ^ 2) * MDP(21) + (t105 * MDP(22) - 0.2e1 * t126) * t83 ^ 2 + (MDP(11) * t94 - 0.2e1 * t123 * MDP(12) + MDP(17) * t146) * t94 + (t86 * t145 + t140) * t82 + 0.2e1 * (t72 * t83 - t74 * t82) * MDP(20) + 0.2e1 * (t66 * t82 + t83 * t68) * MDP(27) + 0.2e1 * (t83 * t141 - t67 * t82) * MDP(28) + (MDP(4) * t110 + 0.2e1 * t113 * MDP(5)) * t110 + 0.2e1 * (-t110 * MDP(10) + t113 * MDP(9)) * pkin(1) + 0.2e1 * t119 * t82; t117 + (-t82 * t90 - t83 * t89) * MDP(20) + (-t72 * t89 + t74 * t90) * MDP(21) + t110 * MDP(6) + t113 * MDP(7) + (t125 * t111 + t68) * MDP(28) + (t125 * t108 - t141) * MDP(27) + (-t113 * MDP(10) - t110 * MDP(9)) * pkin(6); MDP(8) + (t89 ^ 2 + t90 ^ 2) * MDP(21) - t87 * t118 + 0.2e1 * (-MDP(17) * t109 + t134) * pkin(2) + t89 * t145 - 0.2e1 * t90 * MDP(19) + t127; (t124 * t108 - t141) * MDP(27) + (t124 * t111 + t68) * MDP(28) + ((-t106 * t82 - t107 * t83) * MDP(20) + (t106 * t74 - t107 * t72) * MDP(21)) * pkin(3) + t117; -t99 * t133 + t95 * MDP(18) + ((t106 * t90 + t107 * t89) * MDP(21) + t122) * pkin(3) + (t134 + (-MDP(18) * t106 - MDP(19) * t107 - MDP(17)) * t109) * pkin(2) + t127 - t121 * (t87 + t98); -t98 * t118 + t127 + (0.2e1 * t122 + (t106 ^ 2 + t107 ^ 2) * MDP(21) * pkin(3)) * pkin(3); t136 + MDP(21) * t86 + (MDP(18) + t121) * t82; 0; 0; MDP(21); t66 * MDP(27) - t67 * MDP(28) + t119 + t140; t120 * t88 + t135; t120 * t97 + t135; t121; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
