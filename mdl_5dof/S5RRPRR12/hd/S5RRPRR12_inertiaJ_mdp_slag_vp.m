% Calculate joint inertia matrix for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR12_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:38
% EndTime: 2019-12-31 20:30:39
% DurationCPUTime: 0.29s
% Computational Cost: add. (278->94), mult. (503->129), div. (0->0), fcn. (461->6), ass. (0->50)
t85 = sin(qJ(5));
t88 = cos(qJ(5));
t97 = MDP(27) * t85 + MDP(28) * t88;
t119 = pkin(6) - pkin(7);
t87 = sin(qJ(2));
t74 = t119 * t87;
t90 = cos(qJ(2));
t75 = t119 * t90;
t86 = sin(qJ(4));
t89 = cos(qJ(4));
t63 = -t89 * t74 + t86 * t75;
t64 = t86 * t74 + t89 * t75;
t68 = -t90 * t86 + t87 * t89;
t81 = t85 ^ 2;
t83 = t88 ^ 2;
t108 = t88 * MDP(27);
t98 = -t85 * MDP(28) + t108;
t96 = MDP(20) + t98;
t122 = -t96 * t63 - t64 * MDP(21) - (t81 - t83) * t68 * MDP(23);
t121 = -t86 * MDP(21) + t89 * t96;
t73 = -t90 * pkin(2) - t87 * qJ(3) - pkin(1);
t65 = t90 * pkin(3) - t73;
t120 = 0.2e1 * t65;
t91 = -pkin(2) - pkin(3);
t71 = t86 * qJ(3) - t89 * t91;
t69 = pkin(4) + t71;
t118 = pkin(4) + t69;
t117 = t63 * t68;
t82 = t87 ^ 2;
t116 = t90 ^ 2 + t82;
t67 = t87 * t86 + t90 * t89;
t113 = t67 * MDP(26);
t112 = t71 * MDP(20);
t72 = t89 * qJ(3) + t86 * t91;
t111 = t72 * MDP(21);
t109 = t88 * MDP(23);
t107 = t81 * MDP(22) + MDP(19);
t106 = t85 * t109;
t105 = 0.2e1 * t106 + t107;
t104 = -pkin(2) * MDP(14) - MDP(11);
t103 = -t88 * t85 * MDP(22) - MDP(17);
t100 = MDP(24) * t88 - MDP(25) * t85;
t99 = t85 * MDP(24) + t88 * MDP(25);
t95 = 0.2e1 * t98;
t94 = -pkin(8) * t97 + t99;
t93 = -(-pkin(8) + t72) * t97 - t99;
t61 = t67 * pkin(4) - t68 * pkin(8) + t65;
t60 = t85 * t61 + t88 * t64;
t59 = t88 * t61 - t85 * t64;
t1 = [MDP(1) + t82 * MDP(4) + (t116 * pkin(6) ^ 2 + t73 ^ 2) * MDP(14) + (MDP(20) * t120 + t113) * t67 + 0.2e1 * (t85 * t117 + t59 * t67) * MDP(27) + 0.2e1 * (t88 * t117 - t60 * t67) * MDP(28) + 0.2e1 * t116 * MDP(12) * pkin(6) + 0.2e1 * (-t73 * MDP(11) + pkin(1) * MDP(9)) * t90 + 0.2e1 * (-pkin(1) * MDP(10) - t73 * MDP(13) + t90 * MDP(5)) * t87 + (MDP(21) * t120 + 0.2e1 * (-MDP(16) + t100) * t67 + (t83 * MDP(22) + MDP(15) - 0.2e1 * t106) * t68) * t68; t87 * MDP(6) + t90 * MDP(7) + (-t87 * pkin(2) + t90 * qJ(3)) * MDP(12) + (t69 * t97 + t103) * t68 + (MDP(18) + t93) * t67 + ((MDP(14) * qJ(3) - MDP(10) + MDP(13)) * t90 + (-MDP(9) + t104) * t87) * pkin(6) - t122; MDP(8) + 0.2e1 * pkin(2) * MDP(11) + 0.2e1 * qJ(3) * MDP(13) + (pkin(2) ^ 2 + qJ(3) ^ 2) * MDP(14) + t69 * t95 + 0.2e1 * t112 + 0.2e1 * t111 + t105; (pkin(6) * MDP(14) + MDP(12)) * t87 + t97 * (-t67 * t86 - t68 * t89); t104 - t121; MDP(14); (-pkin(4) * t97 - t103) * t68 + (-MDP(18) + t94) * t67 + t122; -t112 - t111 - t118 * t108 + (t118 * MDP(28) - 0.2e1 * t109) * t85 - t107; t121; pkin(4) * t95 + t105; t59 * MDP(27) - t60 * MDP(28) + t100 * t68 + t113; t93; -t97 * t86; t94; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
