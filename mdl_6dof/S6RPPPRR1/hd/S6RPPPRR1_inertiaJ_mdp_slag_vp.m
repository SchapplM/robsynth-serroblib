% Calculate joint inertia matrix for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPPPRR1_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:30:15
% EndTime: 2019-03-09 01:30:15
% DurationCPUTime: 0.18s
% Computational Cost: add. (130->67), mult. (214->85), div. (0->0), fcn. (149->6), ass. (0->31)
t51 = sin(qJ(6));
t53 = cos(qJ(6));
t59 = t53 * MDP(23) - t51 * MDP(24);
t75 = MDP(16) + t59;
t74 = t51 * MDP(23) + t53 * MDP(24);
t49 = sin(pkin(9));
t42 = t49 * pkin(1) + qJ(3);
t73 = t42 ^ 2;
t38 = -pkin(7) + t42;
t71 = t38 * t51;
t70 = t38 * t53;
t69 = t51 * t53;
t54 = cos(qJ(5));
t63 = t54 * MDP(17);
t62 = MDP(7) + MDP(10);
t61 = MDP(19) * t69;
t50 = cos(pkin(9));
t44 = -t50 * pkin(1) - pkin(2);
t39 = qJ(4) - t44;
t60 = MDP(20) * t53 - MDP(21) * t51;
t52 = sin(qJ(5));
t57 = -t75 * t52 - t63;
t56 = t51 * MDP(20) + t53 * MDP(21) - pkin(8) * t74;
t48 = t54 ^ 2;
t47 = t53 ^ 2;
t46 = t52 ^ 2;
t45 = t51 ^ 2;
t37 = t52 * pkin(5) - t54 * pkin(8) + t39;
t36 = t51 * t37 + t52 * t70;
t35 = t53 * t37 - t52 * t71;
t1 = [MDP(1) + (t44 ^ 2 + t73) * MDP(7) + (t39 ^ 2 + t73) * MDP(10) + t46 * MDP(22) + (t49 ^ 2 + t50 ^ 2) * MDP(4) * pkin(1) ^ 2 + (t47 * MDP(18) + MDP(11) - 0.2e1 * t61) * t48 + 0.2e1 * (-MDP(12) + t60) * t54 * t52 + 0.2e1 * t44 * MDP(5) + 0.2e1 * (t35 * t52 - t48 * t71) * MDP(23) + 0.2e1 * (-t36 * t52 - t48 * t70) * MDP(24) + 0.2e1 * (MDP(6) + MDP(8)) * t42 + 0.2e1 * (t52 * MDP(16) + MDP(9) + t63) * t39; 0; MDP(4) + t62; -t39 * MDP(10) + t44 * MDP(7) + MDP(5) - MDP(9) + t57; 0; t62; t42 * MDP(10) + MDP(8) + t74 * (-t46 - t48); 0; 0; MDP(10); (-t38 * MDP(17) - MDP(14) + t56) * t52 + (MDP(13) + t38 * MDP(16) + MDP(18) * t69 + (-t45 + t47) * MDP(19) + (-pkin(5) * t51 + t70) * MDP(23) + (-pkin(5) * t53 - t71) * MDP(24)) * t54; t57; 0; -t52 * MDP(17) + t75 * t54; t45 * MDP(18) + 0.2e1 * pkin(5) * t59 + MDP(15) + 0.2e1 * t61; t52 * MDP(22) + t35 * MDP(23) - t36 * MDP(24) + t60 * t54; -t74 * t54; -t59; -t74 * t52; t56; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
