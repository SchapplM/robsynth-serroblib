% Calculate joint inertia matrix for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPPR1_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:11
% EndTime: 2019-12-05 15:22:11
% DurationCPUTime: 0.20s
% Computational Cost: add. (148->60), mult. (324->87), div. (0->0), fcn. (290->6), ass. (0->30)
t60 = sin(pkin(8));
t62 = cos(pkin(8));
t52 = -t62 * pkin(3) - t60 * qJ(4) - pkin(2);
t61 = cos(pkin(9));
t48 = t61 * t52;
t59 = sin(pkin(9));
t75 = pkin(6) * t60;
t41 = -t61 * t75 + t48 + (-qJ(3) * t59 - pkin(4)) * t62;
t72 = qJ(3) * t62;
t44 = t59 * t52 + t61 * t72;
t42 = -t59 * t75 + t44;
t63 = sin(qJ(5));
t64 = cos(qJ(5));
t50 = t64 * t59 + t63 * t61;
t45 = t50 * t60;
t77 = (t64 * t41 - t63 * t42) * MDP(18) - (t63 * t41 + t64 * t42) * MDP(19) - t45 * MDP(16);
t74 = pkin(2) * MDP(8);
t73 = t59 ^ 2 + t61 ^ 2;
t43 = -t59 * t72 + t48;
t70 = t43 * t61 + t44 * t59;
t49 = -t63 * t59 + t64 * t61;
t69 = t61 * MDP(10) + t59 * MDP(9);
t46 = t49 * t60;
t67 = t45 * MDP(18) + t46 * MDP(19);
t66 = t49 * MDP(18) - t50 * MDP(19);
t65 = qJ(3) ^ 2;
t58 = t62 ^ 2;
t56 = t60 ^ 2;
t54 = t56 * t65;
t1 = [MDP(1) + (t56 + t58) * MDP(8) + (t73 * t56 + t58) * MDP(12); (-t43 * t59 + t44 * t61 - t72) * t60 * MDP(12); MDP(2) + (t58 * t65 + t54) * MDP(8) + (t43 ^ 2 + t44 ^ 2 + t54) * MDP(12) + t58 * MDP(17) + t74 * pkin(2) + (MDP(13) * t46 - 0.2e1 * t45 * MDP(14) - 0.2e1 * t62 * MDP(15)) * t46 + 0.2e1 * (t58 * MDP(7) + (MDP(7) + t69) * t56) * qJ(3) + 0.2e1 * (t44 * MDP(10) + MDP(5) * pkin(2) - t43 * MDP(9) - t77) * t62 + 0.2e1 * (-pkin(2) * MDP(6) + t67 * (pkin(4) * t59 + qJ(3)) - t70 * MDP(11)) * t60; 0; -t74 + t70 * MDP(12) + (-t73 * MDP(11) + MDP(6)) * t60 + (t59 * MDP(10) - t61 * MDP(9) - MDP(5) - t66) * t62; t73 * MDP(12) + MDP(8); -t62 * MDP(12); (MDP(12) * qJ(3) + t69) * t60 + t67; 0; MDP(12); -t67; t46 * MDP(15) - t62 * MDP(17) + t77; t66; 0; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
