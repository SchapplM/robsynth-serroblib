% Calculate joint inertia matrix for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPRR8_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:14
% EndTime: 2019-12-31 18:01:14
% DurationCPUTime: 0.11s
% Computational Cost: add. (137->46), mult. (188->62), div. (0->0), fcn. (151->6), ass. (0->27)
t41 = sin(pkin(8));
t42 = cos(pkin(8));
t44 = sin(qJ(4));
t46 = cos(qJ(4));
t32 = t41 * t46 + t42 * t44;
t43 = sin(qJ(5));
t45 = cos(qJ(5));
t54 = t45 * MDP(18);
t50 = -MDP(19) * t43 + t54;
t61 = -(MDP(11) + t50) * (t41 * t44 - t42 * t46) - t32 * MDP(12);
t60 = -pkin(1) - pkin(2);
t34 = qJ(2) * t41 - t42 * t60;
t33 = -pkin(3) - t34;
t36 = t42 * qJ(2) + t41 * t60;
t29 = -t33 * t46 + t36 * t44;
t27 = pkin(4) + t29;
t59 = pkin(4) + t27;
t58 = t29 * MDP(11);
t30 = t33 * t44 + t36 * t46;
t57 = t30 * MDP(12);
t55 = t45 * MDP(14);
t53 = t43 ^ 2 * MDP(13) + MDP(10);
t52 = 0.2e1 * t43 * t55 + t53;
t51 = -t43 * MDP(15) - t45 * MDP(16);
t49 = -MDP(18) * t43 - MDP(19) * t45;
t47 = 0.2e1 * t50;
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + (t34 ^ 2 + t36 ^ 2) * MDP(9) + t27 * t47 + 0.2e1 * t58 + 0.2e1 * t57 + 0.2e1 * t34 * MDP(7) + 0.2e1 * t36 * MDP(8) + t52; -MDP(4) - pkin(1) * MDP(6) - t42 * MDP(7) + t41 * MDP(8) + (-t34 * t42 + t36 * t41) * MDP(9) - t61; MDP(6) + (t41 ^ 2 + t42 ^ 2) * MDP(9); 0; 0; MDP(9); -t58 - t57 - t59 * t54 + (MDP(19) * t59 - 0.2e1 * t55) * t43 - t53; t61; 0; pkin(4) * t47 + t52; t49 * (-pkin(7) + t30) + t51; t49 * t32; t50; pkin(7) * t49 - t51; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
