% Calculate joint inertia matrix for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP1_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:52:46
% EndTime: 2021-01-15 23:52:47
% DurationCPUTime: 0.26s
% Computational Cost: add. (571->106), mult. (1022->142), div. (0->0), fcn. (1082->6), ass. (0->50)
t87 = cos(qJ(2));
t80 = -t87 * pkin(2) - pkin(1);
t110 = 0.2e1 * t80;
t109 = 0.2e1 * t87;
t108 = 2 * MDP(25);
t107 = pkin(6) + pkin(7);
t83 = sin(qJ(3));
t106 = pkin(2) * t83;
t82 = sin(qJ(4));
t105 = t82 * pkin(3);
t85 = cos(qJ(4));
t81 = t85 * pkin(3);
t88 = 2 * pkin(4);
t104 = t88 + t81;
t103 = MDP(28) * pkin(4);
t78 = t81 + pkin(4);
t102 = MDP(28) * t78;
t86 = cos(qJ(3));
t79 = t86 * pkin(2) + pkin(3);
t77 = t85 * t79;
t67 = -t82 * t106 + t77;
t101 = t67 * MDP(23);
t100 = t85 * MDP(23);
t99 = t86 * MDP(16);
t98 = MDP(15) + MDP(22);
t97 = MDP(24) + MDP(26);
t96 = t85 * t106;
t84 = sin(qJ(2));
t73 = t83 * t87 + t86 * t84;
t75 = t107 * t84;
t76 = t107 * t87;
t94 = -t86 * t75 - t83 * t76;
t56 = -t73 * pkin(8) + t94;
t72 = t83 * t84 - t86 * t87;
t91 = t83 * t75 - t86 * t76;
t57 = -t72 * pkin(8) - t91;
t95 = t85 * t56 - t82 * t57;
t68 = t82 * t79 + t96;
t93 = t97 * t68;
t92 = -t82 * t56 - t85 * t57;
t65 = t72 * pkin(3) + t80;
t61 = -t82 * t72 + t85 * t73;
t50 = -t61 * qJ(5) + t95;
t60 = t85 * t72 + t82 * t73;
t51 = -t60 * qJ(5) - t92;
t90 = t61 * MDP(20) - t60 * MDP(21) + t95 * MDP(23) + t92 * MDP(24) + t50 * MDP(25) - t51 * MDP(26);
t89 = t73 * MDP(13) - t72 * MDP(14) + t94 * MDP(16) + t91 * MDP(17) + t90;
t66 = pkin(4) + t67;
t54 = t60 * pkin(4) + t65;
t1 = [MDP(1) + pkin(1) * MDP(9) * t109 + t72 * MDP(16) * t110 + (t50 ^ 2 + t51 ^ 2 + t54 ^ 2) * MDP(28) + 0.2e1 * (t65 * MDP(23) + t54 * MDP(25) - MDP(27) * t51) * t60 + (MDP(18) * t61 - 0.2e1 * t60 * MDP(19) + 0.2e1 * t65 * MDP(24) + 0.2e1 * t54 * MDP(26) - 0.2e1 * MDP(27) * t50) * t61 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t84 + MDP(5) * t109) * t84 + (MDP(11) * t73 - 0.2e1 * t72 * MDP(12) + MDP(17) * t110) * t73; t89 + (-t87 * MDP(10) - t84 * MDP(9)) * pkin(6) + t84 * MDP(6) + t87 * MDP(7) + (-t68 * t60 - t66 * t61) * MDP(27) + (t50 * t66 + t51 * t68) * MDP(28); MDP(8) + (t66 ^ 2 + t68 ^ 2) * MDP(28) + 0.2e1 * (-t83 * MDP(17) + t99) * pkin(2) + 0.2e1 * t101 + t66 * t108 - 0.2e1 * t93 + t98; (-t60 * t105 - t78 * t61) * MDP(27) + (t51 * t105 + t50 * t78) * MDP(28) + t89; (t77 + t81) * MDP(23) + (t77 + t104) * MDP(25) + (t68 * t105 + t66 * t78) * MDP(28) + t97 * (-t96 + (-pkin(3) - t79) * t82) + (t99 + (-MDP(17) + (-MDP(23) - MDP(25)) * t82) * t83) * pkin(2) + t98; (t108 + t102) * t78 + (0.2e1 * t100 + (MDP(28) * t105 - 0.2e1 * MDP(24) - 0.2e1 * MDP(26)) * t82) * pkin(3) + t98; (-t61 * MDP(27) + t50 * MDP(28)) * pkin(4) + t90; MDP(22) + t101 + (t67 + t88) * MDP(25) + t66 * t103 - t93; MDP(22) + t104 * MDP(25) + pkin(4) * t102 + (-t97 * t82 + t100) * pkin(3); MDP(22) + ((t108 + t103) * pkin(4)); t60 * MDP(25) + t61 * MDP(26) + t54 * MDP(28); 0; 0; 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
