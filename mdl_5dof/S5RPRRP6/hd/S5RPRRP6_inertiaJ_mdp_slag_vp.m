% Calculate joint inertia matrix for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 18:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP6_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 18:08:28
% EndTime: 2021-01-15 18:08:29
% DurationCPUTime: 0.27s
% Computational Cost: add. (268->107), mult. (488->148), div. (0->0), fcn. (392->6), ass. (0->44)
t60 = sin(pkin(8));
t52 = t60 * pkin(1) + pkin(6);
t62 = sin(qJ(4));
t88 = t52 * t62;
t65 = cos(qJ(3));
t87 = t52 * t65;
t64 = cos(qJ(4));
t86 = t62 * t64;
t85 = -qJ(5) - pkin(7);
t84 = MDP(22) * pkin(4);
t63 = sin(qJ(3));
t83 = qJ(5) * t63;
t51 = t85 * t64;
t82 = t51 * MDP(20);
t55 = -t64 * pkin(4) - pkin(3);
t81 = t55 * MDP(22);
t80 = t62 * MDP(15);
t79 = t62 * MDP(20);
t78 = t64 * MDP(19);
t77 = -MDP(18) - MDP(20);
t76 = t64 * t87;
t75 = MDP(13) * t86;
t61 = cos(pkin(8));
t53 = -t61 * pkin(1) - pkin(2);
t74 = -MDP(21) * pkin(4) + MDP(14);
t73 = MDP(19) + t84;
t49 = -t65 * pkin(3) - t63 * pkin(7) + t53;
t47 = t64 * t49;
t43 = -t64 * t83 + t47 + (-pkin(4) - t88) * t65;
t44 = t76 + (t49 - t83) * t62;
t72 = -t43 * t62 + t44 * t64;
t50 = t85 * t62;
t71 = -t50 * t62 - t51 * t64;
t70 = MDP(17) * t62 + MDP(18) * t64;
t69 = t62 * MDP(19) + t64 * MDP(20);
t68 = (-t62 * t87 + t47) * MDP(17) - (t62 * t49 + t76) * MDP(18) - t44 * MDP(20);
t67 = -t78 + t79 + t81;
t59 = t65 ^ 2;
t58 = t64 ^ 2;
t57 = t63 ^ 2;
t56 = t62 ^ 2;
t54 = t58 * t63;
t48 = (pkin(4) * t62 + t52) * t63;
t1 = [MDP(1) - 0.2e1 * t53 * t65 * MDP(10) + t59 * MDP(16) + (t43 ^ 2 + t44 ^ 2 + t48 ^ 2) * MDP(22) + (t60 ^ 2 + t61 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * (t53 * MDP(11) + (-t64 * MDP(14) + MDP(6) + t80) * t65) * t63 + 0.2e1 * (-t43 * MDP(19) - t68) * t65 + 0.2e1 * ((-t43 * t64 - t44 * t62) * MDP(21) + t69 * t48) * t63 + (t58 * MDP(12) + 0.2e1 * t70 * t52 + MDP(5) - 0.2e1 * t75) * t57; (-t48 * t65 + t72 * t63) * MDP(22); MDP(4) + (t59 + (t56 + t58) * t57) * MDP(22); t54 * MDP(13) + t72 * MDP(21) + (t43 * t50 - t44 * t51) * MDP(22) + t67 * t48 + (-t52 * MDP(11) - t62 * MDP(14) - t64 * MDP(15) - t50 * MDP(19) + t70 * pkin(7) + MDP(8) - t82) * t65 + (MDP(7) - t52 * MDP(10) + MDP(12) * t86 - t56 * MDP(13) + (-pkin(3) * t62 - t52 * t64) * MDP(17) + (-pkin(3) * t64 + t88) * MDP(18) + (-t50 * t64 + t51 * t62) * MDP(21) + t69 * t55) * t63; t54 * MDP(21) + (t56 * MDP(21) + t71 * MDP(22) - MDP(11)) * t63 + (-t81 + MDP(10) + (MDP(17) + MDP(19)) * t64 + t77 * t62) * t65; MDP(9) + t56 * MDP(12) + 0.2e1 * t75 + 0.2e1 * t71 * MDP(21) + (t50 ^ 2 + t51 ^ 2) * MDP(22) + (-0.2e1 * t78 + 0.2e1 * t79 + t81) * t55 + 0.2e1 * (t64 * MDP(17) - t62 * MDP(18)) * pkin(3); t43 * t84 + t47 * MDP(19) + (-MDP(16) + (-0.2e1 * pkin(4) - t88) * MDP(19)) * t65 + (-t80 + (-MDP(19) * qJ(5) + t74) * t64) * t63 + t68; (t77 * t64 + (-MDP(17) - t73) * t62) * t63; t82 + (-MDP(18) * pkin(7) + MDP(15)) * t64 + t73 * t50 + (-MDP(17) * pkin(7) + t74) * t62; MDP(16) + (0.2e1 * MDP(19) + t84) * pkin(4); t48 * MDP(22) + t69 * t63; -t65 * MDP(22); t67; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
