% Calculate joint inertia matrix for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRRP6_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:11
% EndTime: 2019-12-31 17:19:11
% DurationCPUTime: 0.16s
% Computational Cost: add. (128->71), mult. (262->107), div. (0->0), fcn. (197->4), ass. (0->31)
t60 = 2 * MDP(18);
t44 = sin(qJ(3));
t59 = pkin(5) * t44;
t46 = cos(qJ(3));
t58 = pkin(5) * t46;
t47 = cos(qJ(2));
t57 = pkin(5) * t47;
t56 = t44 * t46;
t55 = -qJ(4) - pkin(6);
t54 = MDP(19) * pkin(3);
t45 = sin(qJ(2));
t53 = qJ(4) * t45;
t52 = t44 * MDP(14);
t51 = t47 * MDP(15);
t50 = t46 * t57;
t49 = MDP(12) * t56;
t48 = -MDP(18) * pkin(3) + MDP(13);
t43 = t46 ^ 2;
t42 = t45 ^ 2;
t41 = t44 ^ 2;
t40 = -t46 * pkin(3) - pkin(2);
t39 = t55 * t46;
t38 = t55 * t44;
t37 = -t47 * pkin(2) - t45 * pkin(6) - pkin(1);
t36 = (pkin(3) * t44 + pkin(5)) * t45;
t35 = t46 * t37;
t34 = t44 * t37 + t50;
t33 = -t44 * t57 + t35;
t32 = t50 + (t37 - t53) * t44;
t31 = -t46 * t53 + t35 + (-pkin(3) - t59) * t47;
t1 = [MDP(1) + (t31 ^ 2 + t32 ^ 2 + t36 ^ 2) * MDP(19) + (0.2e1 * pkin(1) * MDP(9) + t51) * t47 + (t43 * MDP(11) + MDP(4) - 0.2e1 * t49) * t42 + 0.2e1 * (-t33 * t47 + t42 * t59) * MDP(16) + 0.2e1 * (t34 * t47 + t42 * t58) * MDP(17) + (-0.2e1 * pkin(1) * MDP(10) + (-t31 * t46 - t32 * t44) * t60 + 0.2e1 * (-t46 * MDP(13) + MDP(5) + t52) * t47) * t45; (-t31 * t44 + t32 * t46) * MDP(18) + (t31 * t38 - t32 * t39 + t36 * t40) * MDP(19) + (-pkin(5) * MDP(10) - t44 * MDP(13) - t46 * MDP(14) + MDP(7) + (MDP(16) * t44 + MDP(17) * t46) * pkin(6)) * t47 + (MDP(6) - pkin(5) * MDP(9) + MDP(11) * t56 + (-t41 + t43) * MDP(12) + (-pkin(2) * t44 - t58) * MDP(16) + (-pkin(2) * t46 + t59) * MDP(17) + (-t38 * t46 + t39 * t44) * MDP(18)) * t45; MDP(8) + t41 * MDP(11) + 0.2e1 * t49 + (-t38 * t44 - t39 * t46) * t60 + (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) * MDP(19) + 0.2e1 * (t46 * MDP(16) - t44 * MDP(17)) * pkin(2); t31 * t54 - t51 + t33 * MDP(16) - t34 * MDP(17) + (t46 * t48 - t52) * t45; t38 * t54 + (-MDP(17) * pkin(6) + MDP(14)) * t46 + (-MDP(16) * pkin(6) + t48) * t44; (pkin(3) ^ 2) * MDP(19) + MDP(15); t36 * MDP(19); t40 * MDP(19); 0; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
