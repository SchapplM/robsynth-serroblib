% Calculate joint inertia matrix for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPPR4_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:55
% EndTime: 2019-12-31 17:36:56
% DurationCPUTime: 0.08s
% Computational Cost: add. (68->38), mult. (127->52), div. (0->0), fcn. (95->4), ass. (0->21)
t44 = MDP(8) + MDP(12);
t38 = sin(pkin(8));
t39 = cos(pkin(8));
t46 = t38 ^ 2 + t39 ^ 2;
t53 = t44 * t46;
t43 = t38 * qJ(4) + pkin(2);
t52 = 0.2e1 * (pkin(3) + pkin(4)) * t39 + 0.2e1 * t43;
t51 = -0.2e1 * t38;
t50 = pkin(2) * MDP(8);
t49 = -pkin(6) + qJ(3);
t40 = sin(qJ(5));
t41 = cos(qJ(5));
t27 = t38 * t40 + t39 * t41;
t25 = t27 * MDP(18);
t28 = t38 * t41 - t39 * t40;
t48 = -t28 * MDP(19) - t25;
t29 = -t39 * pkin(3) - t43;
t45 = t29 * MDP(12);
t32 = t49 * t39;
t31 = t49 * t38;
t1 = [MDP(1) + t53; 0; MDP(2) + t25 * t52 + (MDP(11) * t51 - 0.2e1 * t39 * MDP(9) + t45) * t29 + (0.2e1 * t39 * MDP(5) + MDP(6) * t51 + t50) * pkin(2) + (MDP(13) * t28 - 0.2e1 * t27 * MDP(14) + MDP(19) * t52) * t28 + (0.2e1 * (MDP(10) + MDP(7)) * t46 + t53 * qJ(3)) * qJ(3); 0; t45 - t50 + (-MDP(5) - MDP(9)) * t39 + (-MDP(11) + MDP(6)) * t38 + t48; t44; -t39 * MDP(12); (MDP(12) * qJ(3) + MDP(10)) * t38; 0; MDP(12); t48; t28 * MDP(15) - t27 * MDP(16) + (t41 * t31 - t40 * t32) * MDP(18) + (-t40 * t31 - t41 * t32) * MDP(19); 0; t41 * MDP(18) - t40 * MDP(19); MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
