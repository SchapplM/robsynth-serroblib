% Calculate joint inertia matrix for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP5_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:06
% EndTime: 2019-12-31 17:17:06
% DurationCPUTime: 0.16s
% Computational Cost: add. (156->62), mult. (278->90), div. (0->0), fcn. (229->4), ass. (0->24)
t58 = MDP(16) + MDP(18);
t49 = cos(qJ(2));
t57 = 0.2e1 * t49;
t56 = 2 * MDP(18);
t55 = 2 * MDP(20);
t54 = -pkin(6) - pkin(5);
t46 = sin(qJ(3));
t53 = t46 * MDP(17);
t47 = sin(qJ(2));
t52 = t54 * t47;
t44 = -t49 * pkin(2) - pkin(1);
t51 = pkin(3) * t56 + MDP(15);
t39 = t54 * t49;
t48 = cos(qJ(3));
t32 = -t46 * t39 - t48 * t52;
t33 = -t48 * t39 + t46 * t52;
t36 = t46 * t47 - t48 * t49;
t37 = t46 * t49 + t48 * t47;
t50 = t37 * MDP(13) - t36 * MDP(14) + (-MDP(17) + MDP(20)) * t33 - t58 * t32;
t45 = t46 * pkin(2);
t42 = t48 * pkin(2) + pkin(3);
t40 = t45 + qJ(4);
t27 = t36 * pkin(3) - t37 * qJ(4) + t44;
t1 = [MDP(1) + pkin(1) * MDP(9) * t57 + (t27 ^ 2 + t32 ^ 2 + t33 ^ 2) * MDP(21) + 0.2e1 * (t44 * MDP(16) + t27 * MDP(18) - MDP(19) * t33) * t36 + (MDP(11) * t37 - 0.2e1 * t36 * MDP(12) + 0.2e1 * t44 * MDP(17) + 0.2e1 * MDP(19) * t32 - 0.2e1 * t27 * MDP(20)) * t37 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t47 + MDP(5) * t57) * t47; t47 * MDP(6) + t49 * MDP(7) + (-t40 * t36 - t42 * t37) * MDP(19) + (-t32 * t42 + t33 * t40) * MDP(21) + (-t49 * MDP(10) - t47 * MDP(9)) * pkin(5) + t50; MDP(8) + MDP(15) + (t40 ^ 2 + t42 ^ 2) * MDP(21) + 0.2e1 * (t48 * MDP(16) - t53) * pkin(2) + t42 * t56 + t40 * t55; (-pkin(3) * t37 - t36 * qJ(4)) * MDP(19) + (-t32 * pkin(3) + t33 * qJ(4)) * MDP(21) + t50; (0.2e1 * qJ(4) + t45) * MDP(20) + (t42 * pkin(3) + t40 * qJ(4)) * MDP(21) + (t58 * t48 - t53) * pkin(2) + t51; qJ(4) * t55 + ((pkin(3) ^ 2) + qJ(4) ^ 2) * MDP(21) + t51; t37 * MDP(19) + t32 * MDP(21); -t42 * MDP(21) - MDP(18); -pkin(3) * MDP(21) - MDP(18); MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
