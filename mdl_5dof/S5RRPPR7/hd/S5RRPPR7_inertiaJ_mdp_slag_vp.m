% Calculate joint inertia matrix for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR7_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:59:02
% EndTime: 2021-01-15 19:59:04
% DurationCPUTime: 0.26s
% Computational Cost: add. (298->90), mult. (548->124), div. (0->0), fcn. (530->6), ass. (0->45)
t73 = sin(pkin(8));
t74 = cos(pkin(8));
t76 = sin(qJ(2));
t78 = cos(qJ(2));
t61 = t73 * t76 - t74 * t78;
t62 = t73 * t78 + t74 * t76;
t70 = -t78 * pkin(2) - pkin(1);
t83 = -t62 * qJ(4) + t70;
t55 = t61 * pkin(3) + t83;
t102 = -0.2e1 * t55;
t101 = 2 * MDP(11);
t100 = 0.2e1 * t78;
t97 = -qJ(3) - pkin(6);
t64 = t97 * t78;
t88 = t97 * t76;
t58 = -t74 * t64 + t73 * t88;
t54 = -t61 * pkin(4) + t58;
t99 = t54 * t61;
t75 = sin(qJ(5));
t77 = cos(qJ(5));
t98 = t75 * t77;
t66 = t73 * pkin(2) + qJ(4);
t96 = MDP(18) * t66;
t95 = t62 * MDP(23);
t69 = -t74 * pkin(2) - pkin(3);
t94 = t69 * MDP(18);
t93 = t75 * MDP(24);
t92 = t77 * MDP(25);
t91 = MDP(17) - MDP(12);
t90 = MDP(20) * t98;
t56 = -t73 * t64 - t74 * t88;
t89 = t56 ^ 2 + t58 ^ 2;
t87 = MDP(16) + t94;
t85 = t77 * MDP(24) - t75 * MDP(25);
t84 = t92 + t93;
t82 = MDP(15) + t85;
t81 = (MDP(21) * t75 + MDP(22) * t77) * t61;
t80 = t77 * MDP(21) - t75 * MDP(22) + t85 * (-pkin(7) + t69);
t72 = t77 ^ 2;
t71 = t75 ^ 2;
t53 = t62 * pkin(4) + t56;
t52 = (pkin(3) + pkin(7)) * t61 + t83;
t51 = t77 * t52 + t75 * t53;
t50 = -t75 * t52 + t77 * t53;
t1 = [MDP(1) + pkin(1) * MDP(9) * t100 + (t70 ^ 2 + t89) * MDP(14) + (t55 ^ 2 + t89) * MDP(18) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t76 + MDP(5) * t100) * t76 + (0.2e1 * t70 * MDP(12) + MDP(17) * t102 + 0.2e1 * t81 + t95) * t62 + 0.2e1 * (t50 * t62 - t77 * t99) * MDP(24) + 0.2e1 * (-t51 * t62 + t75 * t99) * MDP(25) + 0.2e1 * (MDP(13) + MDP(15)) * (t56 * t62 - t58 * t61) + (t70 * t101 + MDP(16) * t102 + (t71 * MDP(19) + 0.2e1 * t90) * t61) * t61; t76 * MDP(6) + t78 * MDP(7) + (t91 + t96) * t58 + (-MDP(11) + t87) * t56 + t84 * t54 + (-t78 * MDP(10) - t76 * MDP(9)) * pkin(6) + (t69 * MDP(15) + t80) * t62 + (MDP(19) * t98 + (-t71 + t72) * MDP(20) - t82 * t66) * t61 + ((-t61 * t73 - t62 * t74) * MDP(13) + (-t56 * t74 + t58 * t73) * MDP(14)) * pkin(2); -0.2e1 * t90 + t72 * MDP(19) + MDP(8) + (0.2e1 * MDP(16) + t94) * t69 + (0.2e1 * MDP(17) + 0.2e1 * t92 + 0.2e1 * t93 + t96) * t66 + (t74 * t101 - 0.2e1 * t73 * MDP(12) + (t73 ^ 2 + t74 ^ 2) * MDP(14) * pkin(2)) * pkin(2); t70 * MDP(14) + t55 * MDP(18) + (MDP(11) - MDP(16)) * t61 + (-t84 - t91) * t62; 0; MDP(14) + MDP(18); t56 * MDP(18) + t62 * t82; t87; 0; MDP(18); t50 * MDP(24) - t51 * MDP(25) + t81 + t95; t80; -t84; t85; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
