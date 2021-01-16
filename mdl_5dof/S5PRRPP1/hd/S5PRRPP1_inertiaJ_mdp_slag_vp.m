% Calculate joint inertia matrix for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:23
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP1_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:22:41
% EndTime: 2021-01-15 15:22:42
% DurationCPUTime: 0.15s
% Computational Cost: add. (209->62), mult. (384->93), div. (0->0), fcn. (366->4), ass. (0->27)
t77 = 2 * MDP(12);
t76 = qJ(4) + pkin(6);
t75 = MDP(15) * pkin(3);
t74 = MDP(15) + MDP(19);
t73 = 2 * MDP(16);
t72 = cos(qJ(3));
t55 = sin(pkin(8));
t56 = cos(pkin(8));
t57 = sin(qJ(3));
t45 = t55 * t57 - t56 * t72;
t47 = t55 * t72 + t56 * t57;
t54 = -t72 * pkin(3) - pkin(2);
t38 = t45 * pkin(4) - t47 * qJ(5) + t54;
t71 = t38 * MDP(19);
t70 = t54 * MDP(15);
t69 = MDP(13) - MDP(18);
t66 = t76 * t57;
t62 = t72 * MDP(10);
t52 = t56 * pkin(3) + pkin(4);
t61 = -t52 * MDP(19) - MDP(16);
t60 = -MDP(12) + t61;
t50 = t55 * pkin(3) + qJ(5);
t59 = MDP(19) * t50 - t69;
t49 = t76 * t72;
t43 = t56 * t49 - t55 * t66;
t41 = t55 * t49 + t56 * t66;
t1 = [MDP(1) + t74 * (t45 ^ 2 + t47 ^ 2); t74 * (t45 * t41 + t47 * t43); MDP(2) + 0.2e1 * pkin(2) * t62 + (0.2e1 * t47 * MDP(13) + t45 * t77 + t70) * t54 + (-0.2e1 * t47 * MDP(18) + t45 * t73 + t71) * t38 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t57 + 0.2e1 * t72 * MDP(6)) * t57 + t74 * (t41 ^ 2 + t43 ^ 2) + 0.2e1 * (MDP(14) + MDP(17)) * (t41 * t47 - t43 * t45); t62 - t57 * MDP(11) + t59 * t47 + t60 * t45 + (-t45 * t56 + t47 * t55) * t75; t57 * MDP(7) + t72 * MDP(8) + (-t50 * t45 - t52 * t47) * MDP(17) + t59 * t43 + t60 * t41 + (-t57 * MDP(10) - t72 * MDP(11)) * pkin(6) + ((-t45 * t55 - t47 * t56) * MDP(14) + (-t41 * t56 + t43 * t55) * MDP(15)) * pkin(3); MDP(9) + (t50 ^ 2 + t52 ^ 2) * MDP(19) + t52 * t73 + 0.2e1 * t50 * MDP(18) + (t56 * t77 - 0.2e1 * t55 * MDP(13) + (t55 ^ 2 + t56 ^ 2) * t75) * pkin(3); 0; t70 + t71 + t69 * t47 + (MDP(12) + MDP(16)) * t45; 0; t74; t45 * MDP(19); t47 * MDP(17) + t41 * MDP(19); t61; 0; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
