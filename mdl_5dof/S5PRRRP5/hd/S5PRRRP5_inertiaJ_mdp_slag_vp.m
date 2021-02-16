% Calculate joint inertia matrix for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:34
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP5_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:33:44
% EndTime: 2021-01-15 16:33:45
% DurationCPUTime: 0.18s
% Computational Cost: add. (220->77), mult. (431->109), div. (0->0), fcn. (422->6), ass. (0->34)
t83 = -MDP(17) - MDP(19);
t71 = -MDP(18) - MDP(20);
t60 = sin(qJ(4));
t61 = sin(qJ(3));
t63 = cos(qJ(4));
t64 = cos(qJ(3));
t82 = -t60 * t61 + t63 * t64;
t81 = 2 * MDP(19);
t80 = pkin(6) + pkin(7);
t79 = pkin(3) * t60;
t76 = MDP(22) * pkin(4);
t58 = -t64 * pkin(3) - pkin(2);
t41 = -pkin(4) * t82 + t58;
t75 = t41 * MDP(22);
t59 = t63 * pkin(3);
t57 = t59 + pkin(4);
t74 = t57 * MDP(22);
t73 = t63 * MDP(17);
t72 = t64 * MDP(10);
t53 = t80 * t61;
t54 = t80 * t64;
t70 = -t63 * t53 - t60 * t54;
t51 = t60 * t64 + t63 * t61;
t62 = sin(qJ(2));
t46 = t51 * t62;
t47 = t82 * t62;
t69 = t83 * t46 + t71 * t47;
t68 = t60 * t53 - t63 * t54;
t37 = -t51 * qJ(5) + t70;
t38 = qJ(5) * t82 - t68;
t67 = t51 * MDP(14) + MDP(15) * t82 + t70 * MDP(17) + t68 * MDP(18) + t37 * MDP(19) - t38 * MDP(20);
t66 = -t61 * MDP(10) - t64 * MDP(11);
t65 = cos(qJ(2));
t1 = [MDP(1) + (t46 ^ 2 + t47 ^ 2 + t65 ^ 2) * MDP(22); -t62 * MDP(4) + (t46 * t51 + t47 * t82) * MDP(21) + (-t46 * t37 + t47 * t38) * MDP(22) + (-t61 * MDP(11) + t71 * t51 - t82 * t83 + MDP(3) + t72 - t75) * t65; MDP(2) + 0.2e1 * pkin(2) * t72 + (t37 ^ 2 + t38 ^ 2 + t41 ^ 2) * MDP(22) - 0.2e1 * (t58 * MDP(17) + t41 * MDP(19) - MDP(21) * t38) * t82 + (MDP(12) * t51 + 0.2e1 * MDP(13) * t82 + 0.2e1 * t58 * MDP(18) + 0.2e1 * t41 * MDP(20) - 0.2e1 * MDP(21) * t37) * t51 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t61 + 0.2e1 * t64 * MDP(6)) * t61; (-t46 * t57 + t47 * t79) * MDP(22) + t66 * t62 + t69; t61 * MDP(7) + t64 * MDP(8) + (-t57 * t51 + t79 * t82) * MDP(21) + (t37 * t57 + t38 * t79) * MDP(22) + t66 * pkin(6) + t67; MDP(16) + MDP(9) + (t81 + t74) * t57 + (0.2e1 * t73 + (MDP(22) * t79 - 0.2e1 * MDP(18) - 0.2e1 * MDP(20)) * t60) * pkin(3); -t46 * t76 + t69; (-t51 * MDP(21) + t37 * MDP(22)) * pkin(4) + t67; MDP(16) + (0.2e1 * pkin(4) + t59) * MDP(19) + pkin(4) * t74 + (t71 * t60 + t73) * pkin(3); MDP(16) + (t81 + t76) * pkin(4); -t65 * MDP(22); -MDP(19) * t82 + t51 * MDP(20) + t75; 0; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
