% Calculate joint inertia matrix for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP4_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:20
% EndTime: 2019-12-31 17:52:20
% DurationCPUTime: 0.10s
% Computational Cost: add. (120->50), mult. (164->70), div. (0->0), fcn. (115->4), ass. (0->24)
t42 = sin(pkin(7));
t43 = cos(pkin(7));
t46 = -pkin(1) - pkin(2);
t34 = t43 * qJ(2) + t42 * t46;
t44 = sin(qJ(4));
t40 = t44 ^ 2;
t45 = cos(qJ(4));
t55 = t45 ^ 2 + t40;
t54 = MDP(18) * pkin(4);
t31 = -pkin(6) + t34;
t53 = qJ(5) - t31;
t32 = t42 * qJ(2) - t43 * t46;
t49 = pkin(3) + t32;
t29 = t45 * pkin(4) + t49;
t52 = t29 * MDP(18);
t51 = t44 * MDP(16);
t50 = MDP(15) + t54;
t27 = t53 * t44;
t28 = t53 * t45;
t48 = t27 * t44 + t28 * t45;
t47 = t45 * MDP(15) - t51;
t39 = t43 ^ 2;
t38 = t42 ^ 2;
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + (t32 ^ 2 + t34 ^ 2) * MDP(9) + t40 * MDP(10) + 0.2e1 * t44 * t45 * MDP(11) + (t27 ^ 2 + t28 ^ 2 + t29 ^ 2) * MDP(18) + 0.2e1 * t47 * t49 + 0.2e1 * t32 * MDP(7) + 0.2e1 * t34 * MDP(8) + 0.2e1 * t48 * MDP(17); -pkin(1) * MDP(6) - MDP(4) + (-t32 * MDP(9) - MDP(7) - t47 - t52) * t43 + (-t55 * MDP(17) - t48 * MDP(18) + t34 * MDP(9) + MDP(8)) * t42; MDP(6) + (t38 + t39) * MDP(9) + (t55 * t38 + t39) * MDP(18); (t27 * t45 - t28 * t44) * MDP(18); 0; t55 * MDP(18) + MDP(9); t27 * t54 + (-MDP(16) * t31 - MDP(13)) * t45 + (-MDP(15) * t31 + MDP(17) * pkin(4) - MDP(12)) * t44; (-MDP(16) * t45 - t50 * t44) * t42; t50 * t45 - t51; pkin(4) ^ 2 * MDP(18) + MDP(14); t52; -t43 * MDP(18); 0; 0; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
