% Calculate joint inertia matrix for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP3_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:17
% EndTime: 2019-12-05 15:11:18
% DurationCPUTime: 0.15s
% Computational Cost: add. (91->47), mult. (212->75), div. (0->0), fcn. (179->6), ass. (0->26)
t55 = pkin(6) * MDP(16);
t59 = MDP(14) + t55;
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t33 = -pkin(4) * t42 - qJ(5) * t40 - pkin(3);
t52 = MDP(12) - MDP(15);
t58 = t33 * MDP(16) + t52 * t40 - (MDP(11) + MDP(13)) * t42 - MDP(4);
t38 = sin(pkin(8));
t43 = cos(qJ(3));
t57 = t38 * t43;
t35 = t40 ^ 2;
t56 = t42 ^ 2 + t35;
t41 = sin(qJ(3));
t54 = MDP(16) * t41;
t51 = t56 * MDP(14);
t50 = -MDP(16) * pkin(4) - MDP(13);
t39 = cos(pkin(8));
t31 = t39 * t42 + t40 * t57;
t32 = -t39 * t40 + t42 * t57;
t48 = t31 * t40 + t32 * t42;
t47 = -MDP(11) + t50;
t46 = MDP(16) * qJ(5) - t52;
t45 = t40 * t47 + t42 * t46;
t36 = t41 ^ 2;
t34 = t38 ^ 2;
t1 = [MDP(1) + (t39 ^ 2 + t34) * MDP(2) + (t31 ^ 2 + t32 ^ 2 + t34 * t36) * MDP(16); (t48 - t57) * t54; MDP(2) + (t36 * t56 + t43 ^ 2) * MDP(16); (-t43 * MDP(5) + t41 * t58) * t38 + t59 * t48; (t55 * t56 - MDP(5) + t51) * t41 - t58 * t43; MDP(3) + t35 * MDP(6) + (pkin(6) ^ 2 * t56 + t33 ^ 2) * MDP(16) + 0.2e1 * pkin(6) * t51 + 0.2e1 * (MDP(11) * pkin(3) - MDP(13) * t33) * t42 + 0.2e1 * (-MDP(12) * pkin(3) - MDP(15) * t33 + MDP(7) * t42) * t40; t31 * t47 + t32 * t46; t45 * t41; t40 * MDP(8) + t42 * MDP(9) + (-pkin(4) * t40 + qJ(5) * t42) * MDP(14) + t45 * pkin(6); MDP(10) + 0.2e1 * pkin(4) * MDP(13) + 0.2e1 * qJ(5) * MDP(15) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(16); t31 * MDP(16); t40 * t54; t59 * t40; t50; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
