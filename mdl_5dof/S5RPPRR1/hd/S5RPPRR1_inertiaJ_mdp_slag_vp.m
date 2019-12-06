% Calculate joint inertia matrix for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR1_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:16
% EndTime: 2019-12-05 17:38:16
% DurationCPUTime: 0.10s
% Computational Cost: add. (90->46), mult. (130->57), div. (0->0), fcn. (103->4), ass. (0->22)
t37 = -pkin(6) + qJ(2);
t53 = -pkin(7) + t37;
t38 = pkin(1) + qJ(3);
t40 = sin(qJ(5));
t41 = sin(qJ(4));
t42 = cos(qJ(5));
t43 = cos(qJ(4));
t30 = t40 * t41 - t42 * t43;
t31 = t40 * t43 + t42 * t41;
t52 = -t30 * MDP(22) - t31 * MDP(23);
t51 = t38 * MDP(9);
t50 = t31 * MDP(22);
t49 = t41 * MDP(15);
t48 = t43 * MDP(16);
t33 = t53 * t41;
t34 = t53 * t43;
t47 = -t30 * MDP(19) - t31 * MDP(20) + (-t40 * t33 + t42 * t34) * MDP(22) + (-t42 * t33 - t40 * t34) * MDP(23);
t46 = t43 * MDP(15) - t41 * MDP(16);
t45 = (MDP(22) * t42 - MDP(23) * t40) * pkin(4);
t44 = (qJ(2) ^ 2);
t35 = t41 * pkin(4) + t38;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2 + t44) * MDP(6)) + (t44 * MDP(9)) + 0.2e1 * t35 * t50 + (MDP(10) * t43 - 0.2e1 * t41 * MDP(11)) * t43 + (2 * (MDP(5) + MDP(7)) * qJ(2)) + ((2 * MDP(8)) + 0.2e1 * t48 + 0.2e1 * t49 + t51) * t38 + (MDP(17) * t30 + 0.2e1 * t31 * MDP(18) - 0.2e1 * t35 * MDP(23)) * t30; t30 * MDP(23) - (pkin(1) * MDP(6)) + MDP(4) - MDP(8) - t48 - t49 - t50 - t51; MDP(6) + MDP(9); qJ(2) * MDP(9) + MDP(7); 0; MDP(9); t43 * MDP(12) - t41 * MDP(13) + t46 * t37 + t47; 0; t46 + t52; MDP(14) + MDP(21) + 0.2e1 * t45; t47; 0; t52; MDP(21) + t45; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
