% Calculate joint inertia matrix for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:03
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP10_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:02:46
% EndTime: 2021-01-15 21:02:48
% DurationCPUTime: 0.31s
% Computational Cost: add. (299->113), mult. (501->136), div. (0->0), fcn. (384->4), ass. (0->53)
t110 = pkin(3) + pkin(6);
t78 = -pkin(2) - pkin(7);
t75 = sin(qJ(2));
t108 = t75 * pkin(4);
t63 = t110 * t75;
t74 = sin(qJ(4));
t107 = t74 * t63;
t76 = cos(qJ(4));
t106 = t74 * t76;
t77 = cos(qJ(2));
t64 = t110 * t77;
t105 = MDP(24) * pkin(4);
t104 = MDP(25) * pkin(4);
t103 = pkin(2) * MDP(14);
t102 = qJ(5) * t77;
t101 = -qJ(5) + t78;
t67 = t74 * pkin(4) + qJ(3);
t100 = t67 * MDP(25);
t99 = t74 * MDP(20);
t98 = t74 * MDP(22);
t97 = t76 * MDP(18);
t96 = t76 * MDP(21);
t95 = t76 * MDP(23);
t94 = pkin(6) ^ 2 * MDP(14);
t93 = MDP(16) * t106;
t92 = -t75 * qJ(3) - pkin(1);
t58 = t78 * t77 + t92;
t54 = -t74 * t58 + t76 * t63;
t91 = pkin(6) * MDP(14) + MDP(11);
t90 = MDP(12) - t103;
t89 = MDP(22) + t104;
t88 = (-MDP(21) - MDP(23)) * t74;
t87 = MDP(20) * t78 + MDP(17);
t86 = t74 * t102 + t54;
t52 = t86 + t108;
t53 = t107 + (t58 - t102) * t76;
t85 = t52 * t76 + t53 * t74;
t84 = t76 * MDP(20) - t74 * MDP(21);
t83 = t76 * MDP(22) - t74 * MDP(23);
t82 = t54 * MDP(20) - (t76 * t58 + t107) * MDP(21) - t53 * MDP(23);
t81 = t95 + t98 + t100;
t60 = t101 * t74;
t80 = -t60 * MDP(23) + (-MDP(21) * t78 - MDP(18)) * t74;
t73 = t77 ^ 2;
t72 = t76 ^ 2;
t71 = t75 ^ 2;
t70 = t74 ^ 2;
t65 = t70 + t72;
t62 = -t77 * pkin(2) + t92;
t61 = t101 * t76;
t57 = t76 * t77 * pkin(4) + t64;
t56 = t60 * t74 + t61 * t76;
t1 = [MDP(1) + (t52 ^ 2 + t53 ^ 2 + t57 ^ 2) * MDP(25) + t62 ^ 2 * MDP(14) + (t70 * MDP(15) + 0.2e1 * t93 + t94) * t73 + (MDP(19) + MDP(4) + t94) * t71 + 0.2e1 * (t71 + t73) * MDP(11) * pkin(6) + 0.2e1 * (pkin(1) * MDP(9) + t62 * MDP(12) + (t52 * t74 - t53 * t76) * MDP(24) + t84 * t64 + t83 * t57) * t77 + 0.2e1 * (-pkin(1) * MDP(10) - t62 * MDP(13) + (-t74 * MDP(17) + MDP(5) - t97) * t77 + t52 * MDP(22) + t82) * t75; -t85 * MDP(24) + (t52 * t61 + t53 * t60) * MDP(25) + (t96 + t99) * t64 + t81 * t57 + (-pkin(2) * MDP(11) + t61 * MDP(22) + MDP(6) + t87 * t76 + (-MDP(9) + t90) * pkin(6) + t80) * t75 + (MDP(7) - MDP(15) * t106 + (t70 - t72) * MDP(16) + (-t60 * t76 + t61 * t74) * MDP(24) + t83 * t67 + (-MDP(10) + MDP(13)) * pkin(6) + (t84 + t91) * qJ(3)) * t77; MDP(8) + t72 * MDP(15) - 0.2e1 * t93 - 0.2e1 * t56 * MDP(24) + (t60 ^ 2 + t61 ^ 2) * MDP(25) + (0.2e1 * t95 + 0.2e1 * t98 + t100) * t67 + (-0.2e1 * MDP(12) + t103) * pkin(2) + (MDP(14) * qJ(3) + 0.2e1 * MDP(13) + 0.2e1 * t96 + 0.2e1 * t99) * qJ(3); t85 * MDP(25) + ((MDP(20) + MDP(22)) * t76 + t88 + t91) * t75; -t65 * MDP(24) + t56 * MDP(25) + t90; t65 * MDP(25) + MDP(14); t75 * MDP(19) + (t86 + 0.2e1 * t108) * MDP(22) + t52 * t104 + (-t97 + (-MDP(17) + t105) * t74) * t77 + t82; t89 * t61 + (t87 - t105) * t76 + t80; t88 + (MDP(20) + t89) * t76; MDP(19) + (0.2e1 * MDP(22) + t104) * pkin(4); t57 * MDP(25) + t77 * t83; t81; 0; 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
