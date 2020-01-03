% Calculate joint inertia matrix for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPPR3_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:36
% EndTime: 2019-12-31 19:26:37
% DurationCPUTime: 0.12s
% Computational Cost: add. (120->51), mult. (194->60), div. (0->0), fcn. (128->6), ass. (0->25)
t46 = sin(qJ(5));
t48 = cos(qJ(5));
t64 = t46 * MDP(16) + t48 * MDP(17);
t63 = 2 * MDP(8);
t47 = sin(qJ(2));
t62 = pkin(1) * t47;
t49 = cos(qJ(2));
t41 = t49 * pkin(1) + pkin(2);
t44 = sin(pkin(8));
t61 = t44 * t41;
t59 = t49 * MDP(5);
t45 = cos(pkin(8));
t33 = t45 * t41 - t44 * t62;
t32 = -pkin(3) - t33;
t58 = t32 * MDP(10);
t40 = -t45 * pkin(2) - pkin(3);
t57 = t40 * MDP(10);
t54 = MDP(4) + (MDP(11) * t48 - 0.2e1 * MDP(12) * t46) * t48;
t53 = t48 * MDP(13) - t46 * MDP(14);
t52 = t48 * MDP(16) - t46 * MDP(17);
t34 = t45 * t62 + t61;
t51 = (2 * MDP(9)) + 0.2e1 * t64;
t38 = t44 * pkin(2) + qJ(4);
t30 = qJ(4) + t34;
t1 = [MDP(1) + (t33 ^ 2 + t34 ^ 2) * MDP(7) + (t63 + t58) * t32 + 0.2e1 * (-t47 * MDP(6) + t59) * pkin(1) + (t30 * MDP(10) + t51) * t30 + t54; (-0.2e1 * pkin(3) - t33) * MDP(8) + (0.2e1 * qJ(4) + t61) * MDP(9) + (t30 * t38 + t32 * t40) * MDP(10) + ((t33 * t45 + t34 * t44) * MDP(7) - t45 * MDP(8) + t44 * MDP(9)) * pkin(2) + (t59 + (MDP(9) * t45 - MDP(6)) * t47) * pkin(1) + t54 + t64 * (t30 + t38); (t44 ^ 2 + t45 ^ 2) * MDP(7) * pkin(2) ^ 2 + (t63 + t57) * t40 + (MDP(10) * t38 + t51) * t38 + t54; 0; 0; MDP(7) + MDP(10); MDP(8) + t58; MDP(8) + t57; 0; MDP(10); t52 * (-pkin(7) + t32) + t53; t52 * (-pkin(7) + t40) + t53; -t64; t52; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
