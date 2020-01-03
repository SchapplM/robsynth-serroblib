% Calculate joint inertia matrix for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPRR4_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:35
% EndTime: 2019-12-31 16:50:36
% DurationCPUTime: 0.09s
% Computational Cost: add. (72->45), mult. (153->69), div. (0->0), fcn. (112->6), ass. (0->25)
t35 = sin(qJ(4));
t37 = cos(qJ(4));
t40 = MDP(17) * t35 + MDP(18) * t37;
t50 = t35 * MDP(14) + t37 * MDP(15) - t40 * pkin(6);
t33 = sin(pkin(7));
t28 = t33 * pkin(1) + pkin(5);
t49 = t28 * t35;
t48 = t28 * t37;
t38 = cos(qJ(3));
t47 = t28 * t38;
t46 = t35 * t37;
t45 = t38 * MDP(16);
t44 = MDP(13) * t46;
t34 = cos(pkin(7));
t29 = -t34 * pkin(1) - pkin(2);
t43 = MDP(14) * t37 - MDP(15) * t35;
t41 = t37 * MDP(17) - t35 * MDP(18);
t36 = sin(qJ(3));
t32 = t37 ^ 2;
t31 = t36 ^ 2;
t30 = t35 ^ 2;
t27 = -t38 * pkin(3) - t36 * pkin(6) + t29;
t26 = t35 * t27 + t37 * t47;
t25 = t37 * t27 - t35 * t47;
t1 = [MDP(1) + (t33 ^ 2 + t34 ^ 2) * MDP(4) * pkin(1) ^ 2 + (-0.2e1 * t29 * MDP(10) + t45) * t38 + (t32 * MDP(12) + MDP(5) - 0.2e1 * t44) * t31 + 0.2e1 * (-t25 * t38 + t31 * t49) * MDP(17) + 0.2e1 * (t26 * t38 + t31 * t48) * MDP(18) + 0.2e1 * (t29 * MDP(11) + (MDP(6) - t43) * t38) * t36; 0; MDP(4); (-t28 * MDP(11) + MDP(8) - t50) * t38 + (MDP(7) - t28 * MDP(10) + MDP(12) * t46 + (-t30 + t32) * MDP(13) + (-pkin(3) * t35 - t48) * MDP(17) + (-pkin(3) * t37 + t49) * MDP(18)) * t36; -t36 * MDP(11) + (MDP(10) + t41) * t38; t30 * MDP(12) + 0.2e1 * pkin(3) * t41 + MDP(9) + 0.2e1 * t44; t25 * MDP(17) - t26 * MDP(18) + t43 * t36 - t45; -t40 * t36; t50; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
