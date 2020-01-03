% Calculate joint inertia matrix for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4RPRR2_inertiaJ_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:13
% EndTime: 2019-12-31 16:48:13
% DurationCPUTime: 0.05s
% Computational Cost: add. (54->23), mult. (102->32), div. (0->0), fcn. (69->6), ass. (0->18)
t29 = sin(qJ(4));
t31 = cos(qJ(4));
t35 = t31 * MDP(13) - t29 * MDP(14);
t27 = sin(pkin(7));
t43 = pkin(1) * t27;
t41 = t29 * MDP(10) + t31 * MDP(11);
t28 = cos(pkin(7));
t23 = t28 * pkin(1) + pkin(2);
t30 = sin(qJ(3));
t32 = cos(qJ(3));
t20 = t32 * t23 - t30 * t43;
t40 = t20 * MDP(6);
t21 = -t30 * t23 - t32 * t43;
t39 = t21 * MDP(7);
t36 = MDP(5) + (MDP(8) * t29 + 0.2e1 * MDP(9) * t31) * t29;
t34 = -MDP(13) * t29 - MDP(14) * t31;
t18 = -pkin(3) - t20;
t1 = [MDP(1) + (t27 ^ 2 + t28 ^ 2) * MDP(4) * pkin(1) ^ 2 - 0.2e1 * t35 * t18 + 0.2e1 * t40 + 0.2e1 * t39 + t36; 0; MDP(4); t36 + t39 + t40 + t35 * (pkin(3) - t18); 0; 0.2e1 * pkin(3) * t35 + t36; t34 * (pkin(6) - t21) + t41; t35; t34 * pkin(6) + t41; MDP(12);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
