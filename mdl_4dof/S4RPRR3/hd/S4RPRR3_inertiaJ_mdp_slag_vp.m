% Calculate joint inertia matrix for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPRR3_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:18
% EndTime: 2019-12-31 16:49:19
% DurationCPUTime: 0.05s
% Computational Cost: add. (70->32), mult. (125->51), div. (0->0), fcn. (109->6), ass. (0->20)
t38 = cos(pkin(7));
t35 = -t38 * pkin(1) - pkin(2);
t42 = cos(qJ(3));
t49 = -0.2e1 * t42 * pkin(3) + 0.2e1 * t35;
t37 = sin(pkin(7));
t34 = t37 * pkin(1) + pkin(5);
t48 = pkin(6) + t34;
t39 = sin(qJ(4));
t40 = sin(qJ(3));
t41 = cos(qJ(4));
t31 = t39 * t40 - t41 * t42;
t25 = t31 * MDP(17);
t32 = t39 * t42 + t41 * t40;
t47 = -t32 * MDP(18) - t25;
t46 = t42 * MDP(10);
t29 = t48 * t40;
t30 = t48 * t42;
t45 = t32 * MDP(14) - t31 * MDP(15) + (-t41 * t29 - t39 * t30) * MDP(17) + (t39 * t29 - t41 * t30) * MDP(18);
t44 = (MDP(17) * t41 - MDP(18) * t39) * pkin(3);
t1 = [-0.2e1 * t35 * t46 + t25 * t49 + MDP(1) + (t37 ^ 2 + t38 ^ 2) * MDP(4) * pkin(1) ^ 2 + (0.2e1 * t35 * MDP(11) + MDP(5) * t40 + 0.2e1 * t42 * MDP(6)) * t40 + (MDP(12) * t32 - 0.2e1 * t31 * MDP(13) + MDP(18) * t49) * t32; 0; MDP(4); t40 * MDP(7) + t42 * MDP(8) + (-MDP(10) * t40 - MDP(11) * t42) * t34 + t45; -t40 * MDP(11) + t46 + t47; MDP(16) + MDP(9) + 0.2e1 * t44; t45; t47; MDP(16) + t44; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
