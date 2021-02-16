% Calculate joint inertia matrix for
% S5PRRRP3
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
%   see S5PRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP3_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:23:29
% EndTime: 2021-01-15 16:23:29
% DurationCPUTime: 0.13s
% Computational Cost: add. (202->66), mult. (366->94), div. (0->0), fcn. (352->4), ass. (0->32)
t71 = pkin(6) + pkin(7);
t53 = sin(qJ(4));
t54 = sin(qJ(3));
t65 = cos(qJ(4));
t66 = cos(qJ(3));
t43 = t53 * t54 - t65 * t66;
t36 = t43 * MDP(19);
t45 = t53 * t66 + t65 * t54;
t40 = t45 * MDP(20);
t70 = t36 + t40;
t69 = t43 * MDP(17) + t45 * MDP(18);
t68 = 0.2e1 * MDP(19);
t67 = pkin(3) * t53;
t64 = MDP(22) * pkin(4);
t57 = -t66 * pkin(3) - pkin(2);
t35 = t43 * pkin(4) + t57;
t63 = t35 * MDP(22);
t52 = t65 * pkin(3);
t50 = t52 + pkin(4);
t62 = t50 * MDP(22);
t47 = t71 * t54;
t48 = t71 * t66;
t61 = -t65 * t47 - t53 * t48;
t60 = t66 * MDP(10);
t59 = t65 * MDP(17);
t58 = -t69 - t70;
t31 = -t45 * qJ(5) + t61;
t55 = t53 * t47 - t65 * t48;
t32 = -t43 * qJ(5) - t55;
t56 = t45 * MDP(14) - t43 * MDP(15) + t61 * MDP(17) + t55 * MDP(18) + t31 * MDP(19) - t32 * MDP(20);
t42 = t45 ^ 2;
t1 = [MDP(1) + (t43 ^ 2 + t42) * MDP(22); (-t43 * t31 + t45 * t32) * MDP(22); MDP(2) + 0.2e1 * pkin(2) * t60 + t42 * MDP(12) - 0.2e1 * t45 * t43 * MDP(13) + 0.2e1 * (-t31 * t45 - t32 * t43) * MDP(21) + (t31 ^ 2 + t32 ^ 2) * MDP(22) + 0.2e1 * t69 * t57 + (0.2e1 * t36 + 0.2e1 * t40 + t63) * t35 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t54 + 0.2e1 * t66 * MDP(6)) * t54; t60 - t54 * MDP(11) + (-t43 * t50 + t45 * t67) * MDP(22) + t58; t54 * MDP(7) + t66 * MDP(8) + (-t43 * t67 - t50 * t45) * MDP(21) + (t31 * t50 + t32 * t67) * MDP(22) + (-t54 * MDP(10) - t66 * MDP(11)) * pkin(6) + t56; MDP(16) + MDP(9) + (t68 + t62) * t50 + (0.2e1 * t59 + (MDP(22) * t67 - 0.2e1 * MDP(18) - 0.2e1 * MDP(20)) * t53) * pkin(3); -t43 * t64 + t58; (-t45 * MDP(21) + t31 * MDP(22)) * pkin(4) + t56; MDP(16) + (0.2e1 * pkin(4) + t52) * MDP(19) + pkin(4) * t62 + (t59 + (-MDP(18) - MDP(20)) * t53) * pkin(3); MDP(16) + (t68 + t64) * pkin(4); 0; t63 + t70; 0; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
