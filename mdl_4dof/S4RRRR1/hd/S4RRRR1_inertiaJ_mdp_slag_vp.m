% Calculate joint inertia matrix for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RRRR1_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:14
% EndTime: 2019-12-31 17:22:15
% DurationCPUTime: 0.09s
% Computational Cost: add. (103->49), mult. (180->56), div. (0->0), fcn. (122->6), ass. (0->30)
t44 = sin(qJ(4));
t47 = cos(qJ(4));
t56 = MDP(15) * t47;
t50 = -0.2e1 * MDP(16) * t44 + 0.2e1 * t56;
t46 = sin(qJ(2));
t64 = pkin(1) * t46;
t63 = pkin(3) * t44;
t48 = cos(qJ(3));
t62 = t48 * pkin(2);
t61 = t44 * MDP(12) + t47 * MDP(13);
t49 = cos(qJ(2));
t39 = t49 * pkin(1) + pkin(2);
t34 = t48 * t39;
t45 = sin(qJ(3));
t32 = -t45 * t64 + t34;
t60 = t32 * MDP(8);
t33 = -t45 * t39 - t48 * t64;
t59 = t33 * MDP(9);
t58 = t45 * MDP(9);
t57 = t49 * MDP(5);
t55 = MDP(7) + (MDP(10) * t44 + 0.2e1 * MDP(11) * t47) * t44;
t54 = MDP(4) + t55;
t52 = -MDP(15) * t44 - MDP(16) * t47;
t51 = (t48 * MDP(8) - t58) * pkin(2);
t43 = pkin(3) * t47;
t38 = -pkin(3) - t62;
t35 = t38 * t44;
t30 = -pkin(3) - t32;
t29 = t30 * t44;
t1 = [MDP(1) - t30 * t50 + 0.2e1 * (-t46 * MDP(6) + t57) * pkin(1) + 0.2e1 * t60 + 0.2e1 * t59 + t54; (t34 + t62) * MDP(8) + (t35 + t29) * MDP(16) + (-t30 - t38) * t56 + (-pkin(2) - t39) * t58 + (t57 + (-MDP(8) * t45 - MDP(9) * t48 - MDP(6)) * t46) * pkin(1) + t54; -t38 * t50 + 0.2e1 * t51 + t54; t60 + t59 + (-t30 * t47 + t43) * MDP(15) + (t29 - t63) * MDP(16) + t55; (-t38 * t47 + t43) * MDP(15) + (t35 - t63) * MDP(16) + t51 + t55; pkin(3) * t50 + t55; t52 * (pkin(7) - t33) + t61; t52 * (t45 * pkin(2) + pkin(7)) + t61; t52 * pkin(7) + t61; MDP(14);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
