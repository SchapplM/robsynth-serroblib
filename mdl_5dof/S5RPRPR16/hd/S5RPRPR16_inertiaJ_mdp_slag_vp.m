% Calculate joint inertia matrix for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR16_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR16_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR16_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:33
% EndTime: 2019-12-31 18:39:34
% DurationCPUTime: 0.14s
% Computational Cost: add. (139->75), mult. (233->94), div. (0->0), fcn. (153->4), ass. (0->38)
t54 = cos(qJ(3));
t76 = 0.2e1 * t54;
t51 = sin(qJ(5));
t53 = cos(qJ(5));
t60 = t53 * MDP(23) - t51 * MDP(24);
t75 = -2 * MDP(15);
t56 = -pkin(1) - pkin(6);
t74 = pkin(4) - t56;
t73 = (pkin(1) * MDP(6));
t52 = sin(qJ(3));
t41 = t74 * t52;
t72 = t41 * t52;
t71 = t51 * t53;
t47 = t52 ^ 2;
t49 = t54 ^ 2;
t44 = t47 + t49;
t70 = (MDP(17) * pkin(3));
t69 = t56 ^ 2 * MDP(17);
t68 = t51 * MDP(23);
t65 = t53 * MDP(24);
t64 = -MDP(13) + MDP(16);
t63 = MDP(19) * t71;
t62 = MDP(15) - t70;
t40 = t52 * pkin(3) - t54 * qJ(4) + qJ(2);
t61 = MDP(20) * t51 + MDP(21) * t53;
t59 = t65 + t68;
t58 = t56 * MDP(17) - t60;
t55 = -pkin(3) - pkin(7);
t57 = (MDP(23) * t55 + MDP(20)) * t53 + (-MDP(24) * t55 - MDP(21)) * t51;
t48 = t53 ^ 2;
t46 = t51 ^ 2;
t43 = t54 * pkin(3) + t52 * qJ(4);
t42 = t74 * t54;
t39 = t44 * t56;
t38 = t52 * pkin(7) + t40;
t37 = t53 * t38 + t51 * t42;
t36 = -t51 * t38 + t53 * t42;
t1 = [MDP(1) + (-0.2e1 * t54 * MDP(16) + MDP(17) * t40) * t40 + (MDP(22) + MDP(7) + t69) * t49 + (t46 * MDP(18) + 0.2e1 * t63 + t69) * t47 + ((-2 * MDP(4) + t73) * pkin(1)) + (MDP(13) * t76 + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + (0.2e1 * qJ(2) * MDP(12) + t40 * t75 + (-MDP(8) + t61) * t76) * t52 - 0.2e1 * t39 * MDP(14) + 0.2e1 * (t36 * t54 + t53 * t72) * MDP(23) + 0.2e1 * (-t37 * t54 - t51 * t72) * MDP(24); t39 * MDP(17) + MDP(4) - t73 + (-MDP(14) - t60) * t44; t44 * MDP(17) + MDP(6); -t43 * MDP(14) - t59 * t41 + (MDP(9) + (MDP(12) - t62) * t56 + t57) * t54 + (-MDP(10) + MDP(18) * t71 + (-t46 + t48) * MDP(19) + t64 * t56 + t58 * qJ(4)) * t52; t43 * MDP(17) + (MDP(12) - MDP(15)) * t54 + (t59 + t64) * t52; -0.2e1 * t63 + t48 * MDP(18) + MDP(11) + (t75 + t70) * pkin(3) + (MDP(17) * qJ(4) + 0.2e1 * MDP(16) + 0.2e1 * t65 + 0.2e1 * t68) * qJ(4); (MDP(14) - t58) * t54; -t54 * MDP(17); t62; MDP(17); t54 * MDP(22) + t36 * MDP(23) - t37 * MDP(24) + t61 * t52; -t60 * t54; t57; t60; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
