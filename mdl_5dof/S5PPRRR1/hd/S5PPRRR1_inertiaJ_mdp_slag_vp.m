% Calculate joint inertia matrix for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRRR1_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:51
% EndTime: 2019-12-05 15:12:51
% DurationCPUTime: 0.06s
% Computational Cost: add. (75->27), mult. (148->39), div. (0->0), fcn. (153->8), ass. (0->19)
t39 = sin(qJ(5));
t42 = cos(qJ(5));
t48 = t42 * MDP(14) - t39 * MDP(15);
t52 = t39 * MDP(11) + t42 * MDP(12);
t49 = MDP(6) + (0.2e1 * MDP(10) * t42 + MDP(9) * t39) * t39;
t47 = -MDP(14) * t39 - MDP(15) * t42;
t37 = sin(pkin(9));
t38 = cos(pkin(9));
t41 = sin(qJ(3));
t44 = cos(qJ(3));
t29 = -t41 * t37 + t44 * t38;
t30 = t44 * t37 + t41 * t38;
t40 = sin(qJ(4));
t43 = cos(qJ(4));
t28 = t40 * t29 + t43 * t30;
t46 = -t28 * MDP(8) + (-MDP(7) - t48) * (-t43 * t29 + t40 * t30);
t45 = (t43 * MDP(7) - t40 * MDP(8)) * pkin(3);
t33 = -t43 * pkin(3) - pkin(4);
t1 = [MDP(1) + (t37 ^ 2 + t38 ^ 2) * MDP(2); 0; MDP(2); t29 * MDP(4) - t30 * MDP(5) + t46; 0; -0.2e1 * t33 * t48 + MDP(3) + 0.2e1 * t45 + t49; t46; 0; t45 + t49 + t48 * (pkin(4) - t33); 0.2e1 * pkin(4) * t48 + t49; t47 * t28; t48; t47 * (t40 * pkin(3) + pkin(7)) + t52; t47 * pkin(7) + t52; MDP(13);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
