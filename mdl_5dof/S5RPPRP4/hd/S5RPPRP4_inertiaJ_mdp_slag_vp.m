% Calculate joint inertia matrix for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP4_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:13:16
% EndTime: 2021-01-15 17:13:17
% DurationCPUTime: 0.12s
% Computational Cost: add. (141->58), mult. (188->73), div. (0->0), fcn. (137->4), ass. (0->26)
t47 = sin(qJ(4));
t53 = MDP(14) + MDP(16);
t59 = t53 * t47;
t45 = sin(pkin(7));
t46 = cos(pkin(7));
t49 = -pkin(1) - pkin(2);
t36 = t46 * qJ(2) + t45 * t49;
t43 = t47 ^ 2;
t48 = cos(qJ(4));
t58 = t48 ^ 2 + t43;
t57 = MDP(18) * pkin(4);
t34 = -pkin(6) + t36;
t56 = qJ(5) - t34;
t35 = -t45 * qJ(2) + t46 * t49;
t51 = pkin(3) - t35;
t32 = t48 * pkin(4) + t51;
t55 = t32 * MDP(18);
t54 = t47 * MDP(16);
t40 = t48 * MDP(15);
t52 = MDP(15) + t57;
t30 = t56 * t47;
t31 = t56 * t48;
t50 = t30 * t47 + t31 * t48;
t42 = t46 ^ 2;
t41 = t45 ^ 2;
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + (t35 ^ 2 + t36 ^ 2) * MDP(7) + t43 * MDP(8) + 0.2e1 * t47 * t48 * MDP(9) + 0.2e1 * t50 * MDP(17) + (t30 ^ 2 + t31 ^ 2) * MDP(18) + 0.2e1 * (t48 * MDP(13) - t47 * MDP(14)) * t51 + (0.2e1 * t40 - 0.2e1 * t54 + t55) * t32; -pkin(1) * MDP(6) - MDP(4) + (-t58 * MDP(17) - t50 * MDP(18) + t36 * MDP(7)) * t45 + (-t55 + t35 * MDP(7) + (-MDP(13) - MDP(15)) * t48 + t59) * t46; MDP(6) + (t41 + t42) * MDP(7) + (t58 * t41 + t42) * MDP(18); (t30 * t48 - t31 * t47) * MDP(18); 0; t58 * MDP(18) + MDP(7); t31 * MDP(16) + (-MDP(14) * t34 - MDP(11)) * t48 + t52 * t30 + (-MDP(13) * t34 + MDP(17) * pkin(4) - MDP(10)) * t47; (-t53 * t48 + (-MDP(13) - t52) * t47) * t45; t40 + (MDP(13) + t57) * t48 - t59; MDP(12) + (0.2e1 * MDP(15) + t57) * pkin(4); t40 - t54 + t55; -t46 * MDP(18); 0; 0; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
