% Calculate joint inertia matrix for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4PRRR7_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:48
% EndTime: 2019-12-31 16:36:49
% DurationCPUTime: 0.14s
% Computational Cost: add. (75->52), mult. (188->85), div. (0->0), fcn. (168->8), ass. (0->31)
t43 = sin(qJ(4));
t46 = cos(qJ(4));
t62 = -(MDP(17) * t43 + MDP(18) * t46) * pkin(7) + t43 * MDP(14) + t46 * MDP(15);
t61 = pkin(6) * t43;
t60 = pkin(6) * t46;
t47 = cos(qJ(3));
t59 = pkin(6) * t47;
t42 = cos(pkin(4));
t44 = sin(qJ(3));
t41 = sin(pkin(4));
t45 = sin(qJ(2));
t57 = t41 * t45;
t35 = -t42 * t47 + t44 * t57;
t58 = t35 * t44;
t48 = cos(qJ(2));
t56 = t41 * t48;
t55 = t43 * t46;
t54 = t47 * MDP(16);
t53 = MDP(13) * t55;
t52 = MDP(14) * t46 - MDP(15) * t43;
t50 = t46 * MDP(17) - t43 * MDP(18);
t40 = t46 ^ 2;
t39 = t44 ^ 2;
t38 = t43 ^ 2;
t37 = -t47 * pkin(3) - t44 * pkin(7) - pkin(2);
t36 = t42 * t44 + t47 * t57;
t34 = t43 * t37 + t46 * t59;
t33 = t46 * t37 - t43 * t59;
t32 = t36 * t46 - t43 * t56;
t31 = -t36 * t43 - t46 * t56;
t1 = [MDP(1); (-t31 * t47 + t43 * t58) * MDP(17) + (t32 * t47 + t46 * t58) * MDP(18) + (-t45 * MDP(4) + (MDP(10) * t47 - MDP(11) * t44 + MDP(3)) * t48) * t41; MDP(2) + (0.2e1 * pkin(2) * MDP(10) + t54) * t47 + (t40 * MDP(12) + MDP(5) - 0.2e1 * t53) * t39 + 0.2e1 * (-t33 * t47 + t39 * t61) * MDP(17) + 0.2e1 * (t34 * t47 + t39 * t60) * MDP(18) + 0.2e1 * (-pkin(2) * MDP(11) + (MDP(6) - t52) * t47) * t44; -t36 * MDP(11) + (-MDP(10) - t50) * t35; (-pkin(6) * MDP(11) + MDP(8) - t62) * t47 + (MDP(7) - pkin(6) * MDP(10) + MDP(12) * t55 + (-t38 + t40) * MDP(13) + (-pkin(3) * t43 - t60) * MDP(17) + (-pkin(3) * t46 + t61) * MDP(18)) * t44; t38 * MDP(12) + 0.2e1 * pkin(3) * t50 + MDP(9) + 0.2e1 * t53; t31 * MDP(17) - t32 * MDP(18); t33 * MDP(17) - t34 * MDP(18) + t52 * t44 - t54; t62; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
