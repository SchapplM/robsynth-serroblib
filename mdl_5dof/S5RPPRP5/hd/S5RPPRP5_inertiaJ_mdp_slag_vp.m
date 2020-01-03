% Calculate joint inertia matrix for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPRP5_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:44
% EndTime: 2019-12-31 17:53:44
% DurationCPUTime: 0.18s
% Computational Cost: add. (166->66), mult. (289->90), div. (0->0), fcn. (247->4), ass. (0->24)
t60 = (MDP(7) + MDP(11));
t44 = sin(pkin(7));
t59 = -0.2e1 * t44;
t58 = pkin(1) * MDP(7);
t57 = -pkin(6) + qJ(2);
t45 = cos(pkin(7));
t52 = t44 * qJ(3) + pkin(1);
t36 = -t45 * pkin(2) - t52;
t55 = t36 * MDP(11);
t54 = -MDP(18) + MDP(21);
t53 = t57 * t44;
t51 = -pkin(4) * MDP(22) - MDP(19);
t50 = MDP(17) - t51;
t49 = MDP(22) * qJ(5) + t54;
t32 = (pkin(2) + pkin(3)) * t45 + t52;
t47 = cos(qJ(4));
t46 = sin(qJ(4));
t38 = t57 * t45;
t35 = t44 * t47 - t45 * t46;
t34 = t44 * t46 + t45 * t47;
t31 = t47 * t38 + t46 * t53;
t30 = t46 * t38 - t47 * t53;
t29 = t34 * pkin(4) - t35 * qJ(5) + t32;
t1 = [MDP(1) + (t29 ^ 2 + t30 ^ 2 + t31 ^ 2) * MDP(22) + (MDP(10) * t59 - 0.2e1 * t45 * MDP(8) + t55) * t36 + 0.2e1 * (t32 * MDP(17) + t29 * MDP(19)) * t34 + (0.2e1 * t45 * MDP(4) + MDP(5) * t59 + t58) * pkin(1) + (MDP(12) * t35 - 0.2e1 * t34 * MDP(13) + 0.2e1 * t32 * MDP(18) - 0.2e1 * t29 * MDP(21)) * t35 + 0.2e1 * (t30 * t35 - t31 * t34) * MDP(20) + (t60 * qJ(2) + 2 * MDP(6) + 2 * MDP(9)) * (t44 ^ 2 + t45 ^ 2) * qJ(2); t55 - t29 * MDP(22) - t58 + (-MDP(4) - MDP(8)) * t45 + (-MDP(10) + MDP(5)) * t44 + t54 * t35 + (-MDP(17) - MDP(19)) * t34; MDP(22) + t60; (-t46 * t34 - t47 * t35) * MDP(20) + (-t30 * t47 + t31 * t46) * MDP(22) + (qJ(2) * MDP(11) + MDP(9)) * t44; 0; MDP(11) + (t46 ^ 2 + t47 ^ 2) * MDP(22); t35 * MDP(14) - t34 * MDP(15) + (-t35 * pkin(4) - t34 * qJ(5)) * MDP(20) + t49 * t31 - t50 * t30; 0; t49 * t46 + t47 * t50; MDP(16) + 0.2e1 * pkin(4) * MDP(19) + 0.2e1 * qJ(5) * MDP(21) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(22); t35 * MDP(20) + t30 * MDP(22); 0; -t47 * MDP(22); t51; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
