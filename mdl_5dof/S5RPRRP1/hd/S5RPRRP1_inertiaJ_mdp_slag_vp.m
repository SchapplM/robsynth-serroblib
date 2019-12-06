% Calculate joint inertia matrix for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP1_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:58
% EndTime: 2019-12-05 17:59:59
% DurationCPUTime: 0.14s
% Computational Cost: add. (183->60), mult. (284->89), div. (0->0), fcn. (266->4), ass. (0->28)
t50 = sin(qJ(4));
t51 = sin(qJ(3));
t52 = cos(qJ(4));
t53 = cos(qJ(3));
t43 = -t50 * t51 + t52 * t53;
t41 = t43 ^ 2;
t42 = t50 * t53 + t52 * t51;
t68 = t42 ^ 2 + t41;
t66 = pkin(3) * t50;
t54 = -pkin(1) - pkin(6);
t65 = -pkin(7) + t54;
t64 = (pkin(1) * MDP(6));
t63 = t43 * MDP(19) - t42 * MDP(20);
t62 = MDP(22) * pkin(4);
t61 = t51 * pkin(3) + qJ(2);
t60 = t52 * MDP(19);
t44 = t65 * t51;
t45 = t65 * t53;
t59 = -t50 * t44 + t52 * t45;
t56 = -t52 * t44 - t50 * t45;
t58 = t43 * MDP(16) - t42 * MDP(17) + t59 * MDP(19) + t56 * MDP(20);
t31 = -t43 * qJ(5) + t59;
t32 = -t42 * qJ(5) - t56;
t57 = t31 * t43 + t32 * t42;
t48 = t52 * pkin(3) + pkin(4);
t55 = t42 * t66 + t43 * t48;
t35 = t42 * pkin(4) + t61;
t1 = [MDP(1) + t41 * MDP(14) - 0.2e1 * t43 * t42 * MDP(15) - 0.2e1 * t57 * MDP(21) + (t31 ^ 2 + t32 ^ 2 + t35 ^ 2) * MDP(22) + (MDP(7) * t53 - 0.2e1 * t51 * MDP(8)) * t53 + 0.2e1 * (t42 * MDP(19) + t43 * MDP(20)) * t61 + ((-2 * MDP(4) + t64) * pkin(1)) + (0.2e1 * t51 * MDP(12) + 0.2e1 * t53 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2); -t68 * MDP(21) + t57 * MDP(22) + MDP(4) - t64; t68 * MDP(22) + MDP(6); -t55 * MDP(21) + (t31 * t48 + t32 * t66) * MDP(22) + (t54 * MDP(12) + MDP(9)) * t53 + (-t54 * MDP(13) - MDP(10)) * t51 + t58; t53 * MDP(12) - t51 * MDP(13) + MDP(22) * t55 + t63; t48 ^ 2 * MDP(22) + MDP(11) + MDP(18) + (0.2e1 * t60 + (MDP(22) * t66 - 0.2e1 * MDP(20)) * t50) * pkin(3); (-MDP(21) * t43 + t31 * MDP(22)) * pkin(4) + t58; t43 * t62 + t63; t48 * t62 + MDP(18) + (-t50 * MDP(20) + t60) * pkin(3); pkin(4) ^ 2 * MDP(22) + MDP(18); t35 * MDP(22); 0; 0; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
