% Calculate joint inertia matrix for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4PRRR6_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:07
% EndTime: 2019-12-31 16:35:07
% DurationCPUTime: 0.06s
% Computational Cost: add. (61->32), mult. (128->50), div. (0->0), fcn. (119->6), ass. (0->18)
t38 = cos(qJ(3));
t47 = -0.2e1 * t38 * pkin(3) - (2 * pkin(2));
t46 = pkin(5) + pkin(6);
t34 = sin(qJ(4));
t35 = sin(qJ(3));
t37 = cos(qJ(4));
t28 = t34 * t35 - t37 * t38;
t29 = t34 * t38 + t37 * t35;
t36 = sin(qJ(2));
t45 = (-MDP(17) * t29 + MDP(18) * t28) * t36;
t44 = MDP(10) * t38;
t43 = MDP(17) * t28;
t30 = t46 * t35;
t31 = t46 * t38;
t42 = t29 * MDP(14) - t28 * MDP(15) + (-t37 * t30 - t34 * t31) * MDP(17) + (t34 * t30 - t37 * t31) * MDP(18);
t41 = -MDP(10) * t35 - MDP(11) * t38;
t40 = (MDP(17) * t37 - MDP(18) * t34) * pkin(3);
t1 = [MDP(1); -t36 * MDP(4) + (-MDP(11) * t35 - MDP(18) * t29 + MDP(3) - t43 + t44) * cos(qJ(2)); 0.2e1 * pkin(2) * t44 + t43 * t47 + MDP(2) + (-0.2e1 * MDP(11) * pkin(2) + MDP(5) * t35 + 0.2e1 * MDP(6) * t38) * t35 + (MDP(12) * t29 - 0.2e1 * MDP(13) * t28 + MDP(18) * t47) * t29; t41 * t36 + t45; t35 * MDP(7) + t38 * MDP(8) + t41 * pkin(5) + t42; MDP(16) + MDP(9) + 0.2e1 * t40; t45; t42; MDP(16) + t40; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
