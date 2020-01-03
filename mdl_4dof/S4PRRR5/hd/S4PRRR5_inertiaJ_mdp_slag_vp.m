% Calculate joint inertia matrix for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRRR5_inertiaJ_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:46
% EndTime: 2019-12-31 16:33:46
% DurationCPUTime: 0.05s
% Computational Cost: add. (50->23), mult. (99->32), div. (0->0), fcn. (83->6), ass. (0->15)
t31 = sin(qJ(4));
t34 = cos(qJ(4));
t40 = t34 * MDP(13) - t31 * MDP(14);
t44 = t31 * MDP(10) + t34 * MDP(11);
t41 = MDP(5) + (MDP(8) * t31 + 0.2e1 * MDP(9) * t34) * t31;
t39 = -MDP(13) * t31 - MDP(14) * t34;
t32 = sin(qJ(3));
t33 = sin(qJ(2));
t35 = cos(qJ(3));
t36 = cos(qJ(2));
t24 = t32 * t36 + t35 * t33;
t38 = -t24 * MDP(7) + (-MDP(6) - t40) * (t32 * t33 - t35 * t36);
t37 = (t35 * MDP(6) - t32 * MDP(7)) * pkin(2);
t27 = -t35 * pkin(2) - pkin(3);
t1 = [MDP(1); t36 * MDP(3) - t33 * MDP(4) + t38; -0.2e1 * t27 * t40 + MDP(2) + 0.2e1 * t37 + t41; t38; t37 + t41 + t40 * (pkin(3) - t27); 0.2e1 * pkin(3) * t40 + t41; t39 * t24; t39 * (t32 * pkin(2) + pkin(6)) + t44; t39 * pkin(6) + t44; MDP(12);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
