% Calculate joint inertia matrix for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RRPP3_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:36:06
% EndTime: 2021-01-15 10:36:07
% DurationCPUTime: 0.15s
% Computational Cost: add. (150->54), mult. (275->79), div. (0->0), fcn. (236->4), ass. (0->24)
t60 = 2 * MDP(11);
t59 = MDP(14) + MDP(18);
t47 = cos(qJ(2));
t58 = 0.2e1 * t47;
t57 = 2 * MDP(15);
t56 = -qJ(3) - pkin(5);
t44 = sin(pkin(6));
t45 = cos(pkin(6));
t46 = sin(qJ(2));
t35 = t44 * t46 - t45 * t47;
t36 = t44 * t47 + t45 * t46;
t43 = -t47 * pkin(2) - pkin(1);
t30 = t35 * pkin(3) - t36 * qJ(4) + t43;
t55 = t30 * MDP(18);
t54 = t43 * MDP(14);
t53 = MDP(12) - MDP(17);
t51 = t56 * t46;
t41 = t45 * pkin(2) + pkin(3);
t49 = -t41 * MDP(18) - MDP(15);
t39 = t44 * pkin(2) + qJ(4);
t38 = t56 * t47;
t34 = -t45 * t38 + t44 * t51;
t32 = -t44 * t38 - t45 * t51;
t1 = [MDP(1) + pkin(1) * MDP(9) * t58 + (0.2e1 * t36 * MDP(12) + t35 * t60 + t54) * t43 + (-0.2e1 * t36 * MDP(17) + t35 * t57 + t55) * t30 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t46 + MDP(5) * t58) * t46 + t59 * (t32 ^ 2 + t34 ^ 2) + 0.2e1 * (MDP(13) + MDP(16)) * (t32 * t36 - t34 * t35); t46 * MDP(6) + t47 * MDP(7) + (-t39 * t35 - t41 * t36) * MDP(16) + (MDP(18) * t39 - t53) * t34 + (-MDP(11) + t49) * t32 + (-t47 * MDP(10) - t46 * MDP(9)) * pkin(5) + ((-t35 * t44 - t36 * t45) * MDP(13) + (-t32 * t45 + t34 * t44) * MDP(14)) * pkin(2); MDP(8) + (t39 ^ 2 + t41 ^ 2) * MDP(18) + t41 * t57 + 0.2e1 * t39 * MDP(17) + (t45 * t60 - 0.2e1 * t44 * MDP(12) + (t44 ^ 2 + t45 ^ 2) * MDP(14) * pkin(2)) * pkin(2); t54 + t55 + t53 * t36 + (MDP(11) + MDP(15)) * t35; 0; t59; t36 * MDP(16) + t32 * MDP(18); t49; 0; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
