% Calculate joint inertia matrix for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RPRR6_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:42
% EndTime: 2019-12-31 16:52:42
% DurationCPUTime: 0.10s
% Computational Cost: add. (136->44), mult. (258->65), div. (0->0), fcn. (270->6), ass. (0->28)
t52 = sin(pkin(7));
t53 = cos(pkin(7));
t55 = sin(qJ(3));
t69 = cos(qJ(3));
t44 = t55 * t52 - t69 * t53;
t49 = -t53 * pkin(2) - pkin(1);
t71 = 0.2e1 * t44 * pkin(3) + 0.2e1 * t49;
t70 = 0.2e1 * t49;
t68 = pkin(1) * MDP(7);
t67 = pkin(5) + qJ(2);
t65 = t52 * MDP(5);
t64 = t53 * MDP(4);
t45 = t69 * t52 + t55 * t53;
t54 = sin(qJ(4));
t56 = cos(qJ(4));
t39 = t56 * t44 + t54 * t45;
t63 = t39 * MDP(20);
t62 = t44 * MDP(13);
t46 = t67 * t52;
t47 = t67 * t53;
t61 = -t69 * t46 - t55 * t47;
t35 = -t45 * pkin(6) + t61;
t59 = t55 * t46 - t69 * t47;
t36 = -t44 * pkin(6) - t59;
t40 = -t54 * t44 + t56 * t45;
t60 = t40 * MDP(17) - t39 * MDP(18) + (t56 * t35 - t54 * t36) * MDP(20) + (-t54 * t35 - t56 * t36) * MDP(21);
t58 = (MDP(20) * t56 - MDP(21) * t54) * pkin(3);
t1 = [t62 * t70 + t63 * t71 + MDP(1) + (0.2e1 * t64 - 0.2e1 * t65 + t68) * pkin(1) + (MDP(14) * t70 + MDP(8) * t45 - 0.2e1 * t44 * MDP(9)) * t45 + (MDP(15) * t40 - 0.2e1 * t39 * MDP(16) + MDP(21) * t71) * t40 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t52 ^ 2 + t53 ^ 2) * qJ(2); t45 * MDP(14) + t40 * MDP(21) + t62 + t63 - t64 + t65 - t68; MDP(7); t45 * MDP(10) - t44 * MDP(11) + t61 * MDP(13) + t59 * MDP(14) + t60; 0; MDP(12) + MDP(19) + 0.2e1 * t58; t60; 0; MDP(19) + t58; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
