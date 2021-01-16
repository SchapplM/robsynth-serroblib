% Calculate joint inertia matrix for
% S4RRRP4
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
%   see S4RRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP4_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:30:23
% EndTime: 2021-01-15 14:30:24
% DurationCPUTime: 0.11s
% Computational Cost: add. (156->57), mult. (284->83), div. (0->0), fcn. (252->4), ass. (0->24)
t47 = cos(qJ(2));
t56 = 0.2e1 * t47;
t55 = 2 * MDP(18);
t54 = pkin(5) + pkin(6);
t44 = sin(qJ(3));
t53 = pkin(2) * t44;
t46 = cos(qJ(3));
t43 = t46 * pkin(2);
t41 = t43 + pkin(3);
t52 = t41 * MDP(21);
t51 = t46 * MDP(16);
t42 = -t47 * pkin(2) - pkin(1);
t45 = sin(qJ(2));
t38 = t54 * t45;
t39 = t54 * t47;
t50 = -t46 * t38 - t44 * t39;
t49 = t44 * t38 - t46 * t39;
t36 = t44 * t47 + t46 * t45;
t28 = -t36 * qJ(4) + t50;
t35 = t44 * t45 - t46 * t47;
t29 = -t35 * qJ(4) - t49;
t48 = t36 * MDP(13) - t35 * MDP(14) + t50 * MDP(16) + t49 * MDP(17) + t28 * MDP(18) - t29 * MDP(19);
t32 = t35 * pkin(3) + t42;
t1 = [MDP(1) + pkin(1) * MDP(9) * t56 + (t28 ^ 2 + t29 ^ 2 + t32 ^ 2) * MDP(21) + 0.2e1 * (t42 * MDP(16) + t32 * MDP(18) - MDP(20) * t29) * t35 + (MDP(11) * t36 - 0.2e1 * t35 * MDP(12) + 0.2e1 * t42 * MDP(17) + 0.2e1 * t32 * MDP(19) - 0.2e1 * MDP(20) * t28) * t36 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t45 + MDP(5) * t56) * t45; t45 * MDP(6) + t47 * MDP(7) + (-t35 * t53 - t41 * t36) * MDP(20) + (t28 * t41 + t29 * t53) * MDP(21) + (-t47 * MDP(10) - t45 * MDP(9)) * pkin(5) + t48; MDP(15) + MDP(8) + (t55 + t52) * t41 + (0.2e1 * t51 + (MDP(21) * t53 - 0.2e1 * MDP(17) - 0.2e1 * MDP(19)) * t44) * pkin(2); (-t36 * MDP(20) + t28 * MDP(21)) * pkin(3) + t48; MDP(15) + (0.2e1 * pkin(3) + t43) * MDP(18) + pkin(3) * t52 + (t51 + (-MDP(17) - MDP(19)) * t44) * pkin(2); MDP(15) + (MDP(21) * pkin(3) + t55) * pkin(3); t35 * MDP(18) + t36 * MDP(19) + t32 * MDP(21); 0; 0; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
