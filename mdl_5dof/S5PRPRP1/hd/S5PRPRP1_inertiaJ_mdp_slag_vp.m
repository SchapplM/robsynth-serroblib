% Calculate joint inertia matrix for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRP1_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:40
% EndTime: 2019-12-05 15:28:41
% DurationCPUTime: 0.12s
% Computational Cost: add. (147->51), mult. (271->74), div. (0->0), fcn. (258->4), ass. (0->24)
t51 = cos(qJ(4));
t50 = pkin(2) * MDP(8);
t49 = pkin(6) + qJ(3);
t36 = sin(pkin(8));
t37 = cos(pkin(8));
t48 = t36 ^ 2 + t37 ^ 2;
t47 = t36 * MDP(6);
t46 = t37 * MDP(5);
t45 = MDP(15) - MDP(18);
t33 = -t37 * pkin(3) - pkin(2);
t44 = t49 * t36;
t43 = t48 * MDP(8);
t42 = -pkin(4) * MDP(19) - MDP(16);
t41 = -MDP(14) + t42;
t40 = MDP(19) * qJ(5) - t45;
t38 = sin(qJ(4));
t31 = t49 * t37;
t30 = t36 * t51 + t38 * t37;
t29 = t38 * t36 - t37 * t51;
t27 = t30 ^ 2;
t26 = t31 * t51 - t38 * t44;
t25 = t38 * t31 + t44 * t51;
t24 = t29 * pkin(4) - t30 * qJ(5) + t33;
t1 = [MDP(1) + t43 + (t29 ^ 2 + t27) * MDP(19); (t29 * t25 + t30 * t26) * MDP(19); MDP(2) + t27 * MDP(9) + (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) * MDP(19) + (0.2e1 * t46 - 0.2e1 * t47 + t50) * pkin(2) + 0.2e1 * (t25 * t30 - t26 * t29) * MDP(17) + (0.2e1 * MDP(7) * t48 + qJ(3) * t43) * qJ(3) + 0.2e1 * (t33 * MDP(15) - t24 * MDP(18)) * t30 + 0.2e1 * (-t30 * MDP(10) + t33 * MDP(14) + t24 * MDP(16)) * t29; 0; t24 * MDP(19) - t46 + t47 - t50 + t45 * t30 + (MDP(14) + MDP(16)) * t29; MDP(8) + MDP(19); t29 * t41 + t30 * t40; t30 * MDP(11) - t29 * MDP(12) + (-pkin(4) * t30 - t29 * qJ(5)) * MDP(17) + t40 * t26 + t41 * t25; 0; MDP(13) + 0.2e1 * pkin(4) * MDP(16) + 0.2e1 * qJ(5) * MDP(18) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(19); t29 * MDP(19); t30 * MDP(17) + t25 * MDP(19); 0; t42; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
