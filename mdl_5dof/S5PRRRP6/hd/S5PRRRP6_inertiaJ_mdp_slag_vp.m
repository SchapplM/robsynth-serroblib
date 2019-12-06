% Calculate joint inertia matrix for
% S5PRRRP6
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
%   see S5PRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP6_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:19
% EndTime: 2019-12-05 16:52:20
% DurationCPUTime: 0.24s
% Computational Cost: add. (221->81), mult. (431->119), div. (0->0), fcn. (406->6), ass. (0->31)
t72 = MDP(17) + MDP(19);
t78 = MDP(18) - MDP(21);
t77 = 2 * MDP(19);
t76 = 2 * MDP(21);
t75 = -pkin(7) - pkin(6);
t61 = sin(qJ(4));
t74 = t61 * MDP(18);
t65 = cos(qJ(3));
t73 = t65 * MDP(10);
t62 = sin(qJ(3));
t71 = t75 * t62;
t59 = -t65 * pkin(3) - pkin(2);
t64 = cos(qJ(4));
t52 = t61 * t65 + t64 * t62;
t63 = sin(qJ(2));
t47 = t52 * t63;
t51 = t61 * t62 - t64 * t65;
t48 = t51 * t63;
t70 = -t72 * t47 + t78 * t48;
t69 = pkin(4) * t77 + MDP(16);
t54 = t75 * t65;
t41 = -t61 * t54 - t64 * t71;
t42 = -t64 * t54 + t61 * t71;
t68 = t52 * MDP(14) - t51 * MDP(15) - t72 * t41 - t78 * t42;
t67 = -t62 * MDP(10) - t65 * MDP(11);
t66 = cos(qJ(2));
t60 = t61 * pkin(3);
t57 = t64 * pkin(3) + pkin(4);
t55 = t60 + qJ(5);
t36 = t51 * pkin(4) - t52 * qJ(5) + t59;
t1 = [MDP(1) + (t47 ^ 2 + t48 ^ 2 + t66 ^ 2) * MDP(22); -t63 * MDP(4) + (t47 * t52 + t48 * t51) * MDP(20) + (t47 * t41 - t48 * t42) * MDP(22) + (-t62 * MDP(11) - t36 * MDP(22) - t72 * t51 - t78 * t52 + MDP(3) + t73) * t66; MDP(2) + 0.2e1 * pkin(2) * t73 + (t36 ^ 2 + t41 ^ 2 + t42 ^ 2) * MDP(22) + 0.2e1 * (t59 * MDP(17) + t36 * MDP(19) - MDP(20) * t42) * t51 + (MDP(12) * t52 - 0.2e1 * t51 * MDP(13) + 0.2e1 * t59 * MDP(18) + 0.2e1 * MDP(20) * t41 - 0.2e1 * t36 * MDP(21)) * t52 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t62 + 0.2e1 * t65 * MDP(6)) * t62; (-t47 * t57 - t48 * t55) * MDP(22) + t67 * t63 + t70; t62 * MDP(7) + t65 * MDP(8) + (-t55 * t51 - t57 * t52) * MDP(20) + (-t41 * t57 + t42 * t55) * MDP(22) + t67 * pkin(6) + t68; MDP(9) + MDP(16) + (t55 ^ 2 + t57 ^ 2) * MDP(22) + 0.2e1 * (t64 * MDP(17) - t74) * pkin(3) + t57 * t77 + t55 * t76; (-t47 * pkin(4) - t48 * qJ(5)) * MDP(22) + t70; (-pkin(4) * t52 - t51 * qJ(5)) * MDP(20) + (-t41 * pkin(4) + t42 * qJ(5)) * MDP(22) + t68; (0.2e1 * qJ(5) + t60) * MDP(21) + (t57 * pkin(4) + t55 * qJ(5)) * MDP(22) + (t72 * t64 - t74) * pkin(3) + t69; qJ(5) * t76 + ((pkin(4) ^ 2) + qJ(5) ^ 2) * MDP(22) + t69; t47 * MDP(22); t52 * MDP(20) + t41 * MDP(22); -t57 * MDP(22) - MDP(19); -MDP(22) * pkin(4) - MDP(19); MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
