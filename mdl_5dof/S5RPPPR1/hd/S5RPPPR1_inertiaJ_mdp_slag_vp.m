% Calculate joint inertia matrix for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPPPR1_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:45
% EndTime: 2022-01-20 09:12:46
% DurationCPUTime: 0.17s
% Computational Cost: add. (184->60), mult. (351->90), div. (0->0), fcn. (318->8), ass. (0->32)
t70 = cos(pkin(7));
t60 = -t70 * pkin(1) - pkin(2);
t66 = sin(pkin(8));
t69 = cos(pkin(8));
t53 = -t69 * pkin(3) - t66 * qJ(4) + t60;
t68 = cos(pkin(9));
t51 = t68 * t53;
t67 = sin(pkin(7));
t59 = t67 * pkin(1) + qJ(3);
t65 = sin(pkin(9));
t81 = pkin(6) * t66;
t44 = -t68 * t81 + t51 + (-t59 * t65 - pkin(4)) * t69;
t80 = t59 * t69;
t47 = t65 * t53 + t68 * t80;
t45 = -t65 * t81 + t47;
t71 = sin(qJ(5));
t72 = cos(qJ(5));
t55 = t72 * t65 + t71 * t68;
t48 = t55 * t66;
t83 = (t72 * t44 - t71 * t45) * MDP(16) - (t71 * t44 + t72 * t45) * MDP(17) - t48 * MDP(14);
t79 = t65 ^ 2 + t68 ^ 2;
t54 = -t71 * t65 + t72 * t68;
t77 = t65 * MDP(8) + t68 * MDP(9);
t49 = t54 * t66;
t75 = t48 * MDP(16) + t49 * MDP(17);
t74 = t54 * MDP(16) - t55 * MDP(17);
t64 = t69 ^ 2;
t62 = t66 ^ 2;
t58 = t59 ^ 2;
t57 = t62 * t58;
t46 = -t65 * t80 + t51;
t1 = [MDP(1) + (t64 * t58 + t60 ^ 2 + t57) * MDP(7) + (t46 ^ 2 + t47 ^ 2 + t57) * MDP(10) + t64 * MDP(15) + (t67 ^ 2 + t70 ^ 2) * MDP(4) * pkin(1) ^ 2 + (MDP(11) * t49 - 0.2e1 * t48 * MDP(12) - 0.2e1 * t69 * MDP(13)) * t49 + 0.2e1 * t75 * (pkin(4) * t65 + t59) * t66 + 0.2e1 * (t64 * MDP(6) + (MDP(6) + t77) * t62) * t59 + 0.2e1 * (-t60 * MDP(5) - t46 * MDP(8) + t47 * MDP(9) - t83) * t69; (-t46 * t65 + t47 * t68 - t80) * t66 * MDP(10); MDP(4) + (t62 + t64) * MDP(7) + (t79 * t62 + t64) * MDP(10); t60 * MDP(7) + (t46 * t68 + t47 * t65) * MDP(10) + (-t68 * MDP(8) + t65 * MDP(9) - MDP(5) - t74) * t69; 0; t79 * MDP(10) + MDP(7); (MDP(10) * t59 + t77) * t66 + t75; -t69 * MDP(10); 0; MDP(10); t49 * MDP(13) - t69 * MDP(15) + t83; -t75; t74; 0; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
