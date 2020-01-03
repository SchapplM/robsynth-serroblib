% Calculate joint inertia matrix for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPPR3_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:57
% EndTime: 2019-12-31 16:37:57
% DurationCPUTime: 0.05s
% Computational Cost: add. (56->28), mult. (103->45), div. (0->0), fcn. (84->6), ass. (0->22)
t36 = cos(pkin(6));
t30 = -t36 * pkin(1) - pkin(2);
t35 = cos(pkin(7));
t48 = -0.2e1 * t35 * pkin(3) + 0.2e1 * t30;
t34 = sin(pkin(6));
t28 = t34 * pkin(1) + qJ(3);
t47 = pkin(5) + t28;
t33 = sin(pkin(7));
t46 = t33 ^ 2 + t35 ^ 2;
t45 = t30 * MDP(8);
t44 = t33 * MDP(6);
t43 = t35 * MDP(5);
t37 = sin(qJ(4));
t38 = cos(qJ(4));
t24 = t37 * t33 - t38 * t35;
t42 = t24 * MDP(14);
t41 = t46 * MDP(8);
t25 = t38 * t33 + t37 * t35;
t40 = -t25 * MDP(15) - t42;
t23 = t47 * t35;
t22 = t47 * t33;
t1 = [t42 * t48 + MDP(1) + (t34 ^ 2 + t36 ^ 2) * MDP(4) * pkin(1) ^ 2 + (-0.2e1 * t43 + 0.2e1 * t44 + t45) * t30 + (-0.2e1 * t24 * MDP(10) + MDP(15) * t48 + MDP(9) * t25) * t25 + (0.2e1 * t46 * MDP(7) + t41 * t28) * t28; 0; MDP(4) + t41; -t40 - t43 + t44 + t45; 0; MDP(8); t25 * MDP(11) - t24 * MDP(12) + (-t38 * t22 - t37 * t23) * MDP(14) + (t37 * t22 - t38 * t23) * MDP(15); t40; 0; MDP(13);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
