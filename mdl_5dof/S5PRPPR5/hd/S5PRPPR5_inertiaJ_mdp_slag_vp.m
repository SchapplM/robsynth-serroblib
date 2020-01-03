% Calculate joint inertia matrix for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRPPR5_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:38
% EndTime: 2019-12-31 17:38:38
% DurationCPUTime: 0.06s
% Computational Cost: add. (84->43), mult. (130->62), div. (0->0), fcn. (112->6), ass. (0->19)
t45 = -pkin(2) - pkin(3);
t33 = sin(pkin(8));
t34 = cos(pkin(8));
t27 = t34 * qJ(3) + t33 * t45;
t37 = cos(qJ(5));
t44 = t37 * MDP(16);
t43 = -pkin(2) * MDP(7) - MDP(5);
t25 = t33 * qJ(3) - t34 * t45;
t42 = t27 * MDP(10) + MDP(9);
t35 = sin(qJ(5));
t41 = -t35 * MDP(17) + t44;
t40 = -MDP(16) * t35 - MDP(17) * t37;
t39 = -t25 * MDP(10) - MDP(8) - t41;
t38 = cos(qJ(2));
t36 = sin(qJ(2));
t23 = pkin(4) + t25;
t22 = -t38 * t33 + t36 * t34;
t21 = -t36 * t33 - t38 * t34;
t1 = [MDP(1) + (t36 ^ 2 + t38 ^ 2) * MDP(7) + (t21 ^ 2 + t22 ^ 2) * MDP(10); (MDP(3) - t43) * t38 + (qJ(3) * MDP(7) - MDP(4) + MDP(6)) * t36 + t42 * t22 + t39 * t21; MDP(2) + (2 * pkin(2) * MDP(5)) + 0.2e1 * qJ(3) * MDP(6) + ((pkin(2) ^ 2) + qJ(3) ^ 2) * MDP(7) + (t25 ^ 2 + t27 ^ 2) * MDP(10) + 0.2e1 * t23 * t44 + 0.2e1 * t25 * MDP(8) + 0.2e1 * t27 * MDP(9) + (MDP(11) * t35 + 0.2e1 * t37 * MDP(12) - 0.2e1 * t23 * MDP(17)) * t35; -t38 * MDP(7) + (t21 * t34 + t22 * t33) * MDP(10); t42 * t33 + t39 * t34 + t43; MDP(7) + (t33 ^ 2 + t34 ^ 2) * MDP(10); 0; 0; 0; MDP(10); t40 * t22; -t35 * MDP(13) - t37 * MDP(14) + t40 * (-pkin(6) + t27); t40 * t33; t41; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
