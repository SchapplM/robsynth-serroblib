% Calculate joint inertia matrix for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRRP10_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:05
% EndTime: 2019-12-31 18:52:06
% DurationCPUTime: 0.21s
% Computational Cost: add. (341->88), mult. (648->125), div. (0->0), fcn. (655->6), ass. (0->46)
t73 = cos(pkin(8));
t66 = -t73 * pkin(2) - pkin(1);
t101 = 0.2e1 * t66;
t100 = cos(qJ(3));
t99 = pkin(1) * MDP(7);
t74 = sin(qJ(4));
t76 = cos(qJ(4));
t98 = t74 * t76;
t72 = sin(pkin(8));
t96 = pkin(6) + qJ(2);
t61 = t96 * t72;
t62 = t96 * t73;
t75 = sin(qJ(3));
t57 = t100 * t62 - t75 * t61;
t97 = t76 * t57;
t95 = -qJ(5) - pkin(7);
t70 = t74 ^ 2;
t71 = t76 ^ 2;
t93 = t70 + t71;
t92 = MDP(23) * pkin(4);
t60 = t100 * t72 + t75 * t73;
t91 = qJ(5) * t60;
t90 = t72 * MDP(5);
t89 = t73 * MDP(4);
t88 = MDP(18) * t74;
t59 = -t100 * t73 + t75 * t72;
t87 = t59 * MDP(19);
t86 = t74 * MDP(21);
t85 = MDP(16) * t98;
t55 = t59 * pkin(3) - t60 * pkin(7) + t66;
t51 = t76 * t55 - t74 * t57;
t84 = -MDP(22) * pkin(4) + MDP(17);
t49 = t59 * pkin(4) - t76 * t91 + t51;
t50 = t97 + (t55 - t91) * t74;
t83 = t49 * t76 + t50 * t74;
t63 = t95 * t74;
t64 = t95 * t76;
t82 = t76 * t63 - t74 * t64;
t81 = t51 * MDP(20) - (t74 * t55 + t97) * MDP(21);
t80 = t76 * MDP(20) - t86;
t79 = MDP(20) * t74 + MDP(21) * t76;
t56 = t100 * t61 + t75 * t62;
t78 = MDP(13) + t80;
t67 = -t76 * pkin(4) - pkin(3);
t53 = t74 * t60 * pkin(4) + t56;
t1 = [MDP(1) + (t49 ^ 2 + t50 ^ 2 + t53 ^ 2) * MDP(23) + (0.2e1 * t89 - 0.2e1 * t90 + t99) * pkin(1) + (MDP(14) * t101 - 0.2e1 * t83 * MDP(22) + 0.2e1 * t79 * t56 + (t71 * MDP(15) + MDP(8) - 0.2e1 * t85) * t60) * t60 + (MDP(13) * t101 + t87 + 0.2e1 * (MDP(17) * t76 - MDP(9) - t88) * t60 + 0.2e1 * t81) * t59 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t72 ^ 2 + t73 ^ 2) * qJ(2); -t89 + t90 - t99 + t83 * MDP(23) + (-t93 * MDP(22) + MDP(14)) * t60 + t78 * t59; t93 * MDP(23) + MDP(7); -t57 * MDP(14) + (-t49 * t74 + t50 * t76) * MDP(22) + (t49 * t63 - t50 * t64 + t53 * t67) * MDP(23) - t78 * t56 + (t74 * MDP(17) + t76 * MDP(18) - t79 * pkin(7) - MDP(11)) * t59 + (MDP(10) + MDP(15) * t98 + (-t70 + t71) * MDP(16) - t82 * MDP(22) - t79 * pkin(3)) * t60; t82 * MDP(23); MDP(12) + t70 * MDP(15) + 0.2e1 * t85 + 0.2e1 * (-t63 * t74 - t64 * t76) * MDP(22) + (t63 ^ 2 + t64 ^ 2 + t67 ^ 2) * MDP(23) + 0.2e1 * t80 * pkin(3); t49 * t92 + t87 + (t84 * t76 - t88) * t60 + t81; -t86 + (MDP(20) + t92) * t76; t63 * t92 + (-MDP(21) * pkin(7) + MDP(18)) * t76 + (-MDP(20) * pkin(7) + t84) * t74; MDP(23) * pkin(4) ^ 2 + MDP(19); t53 * MDP(23); 0; t67 * MDP(23); 0; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
