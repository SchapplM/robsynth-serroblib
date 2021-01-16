% Calculate joint inertia matrix for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRRP1_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:14:41
% EndTime: 2021-01-15 16:14:42
% DurationCPUTime: 0.17s
% Computational Cost: add. (138->65), mult. (229->81), div. (0->0), fcn. (154->4), ass. (0->35)
t44 = sin(qJ(4));
t39 = t44 * MDP(16);
t46 = cos(qJ(4));
t52 = -t46 * MDP(15) + t39;
t66 = 2 * MDP(17);
t47 = cos(qJ(3));
t65 = t47 * pkin(2);
t37 = -pkin(3) - t65;
t64 = pkin(3) - t37;
t63 = qJ(5) + pkin(7);
t38 = -t46 * pkin(4) - pkin(3);
t32 = t38 - t65;
t62 = t32 + t38;
t61 = t44 * MDP(10) + t46 * MDP(11);
t60 = MDP(18) * pkin(4);
t45 = sin(qJ(3));
t36 = t45 * pkin(2) + pkin(7);
t59 = qJ(5) + t36;
t58 = t32 * MDP(18);
t57 = t38 * MDP(18);
t56 = t44 * MDP(14);
t55 = t44 * MDP(17);
t43 = t44 ^ 2;
t53 = 0.2e1 * t44 * t46 * MDP(9) + t43 * MDP(8) + MDP(5);
t51 = t46 * MDP(13) - t56;
t50 = -MDP(13) * t44 - MDP(14) * t46;
t49 = 0.2e1 * t52;
t48 = (t47 * MDP(6) - t45 * MDP(7)) * pkin(2);
t34 = t63 * t46;
t33 = t63 * t44;
t31 = t34 * t46;
t30 = t59 * t46;
t29 = t59 * t44;
t28 = t30 * t46;
t1 = [MDP(1) + (t46 ^ 2 + t43) * MDP(18); (-t46 * t29 + t44 * t30) * MDP(18); MDP(2) + (t29 * t44 + t28) * t66 + (t29 ^ 2 + t30 ^ 2) * MDP(18) + (t49 + t58) * t32 + t53 - 0.2e1 * t51 * t37 + 0.2e1 * t48; (-t46 * t33 + t44 * t34) * MDP(18); (t28 + t31) * MDP(17) + (t29 * t33 + t30 * t34 + t32 * t38) * MDP(18) + t48 + (t64 * MDP(13) - t62 * MDP(15)) * t46 + (-t64 * MDP(14) + t62 * MDP(16) + (t29 + t33) * MDP(17)) * t44 + t53; (t33 * t44 + t31) * t66 + (t33 ^ 2 + t34 ^ 2) * MDP(18) + (t49 + t57) * t38 + 0.2e1 * t51 * pkin(3) + t53; -t56 - t39 + (MDP(13) + MDP(15) + t60) * t46; -t29 * MDP(15) - t30 * MDP(16) + t50 * t36 + (-t29 * MDP(18) - t55) * pkin(4) + t61; -t33 * MDP(15) - t34 * MDP(16) + t50 * pkin(7) + (-MDP(18) * t33 - t55) * pkin(4) + t61; MDP(12) + (0.2e1 * MDP(15) + t60) * pkin(4); 0; t52 + t58; t52 + t57; 0; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
