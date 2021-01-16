% Calculate joint inertia matrix for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP3_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:04:42
% EndTime: 2021-01-15 17:04:42
% DurationCPUTime: 0.09s
% Computational Cost: add. (99->44), mult. (139->61), div. (0->0), fcn. (88->4), ass. (0->20)
t42 = MDP(18) * pkin(4);
t35 = cos(qJ(4));
t31 = t35 ^ 2;
t34 = sin(qJ(4));
t26 = -t34 ^ 2 - t31;
t41 = -t26 * MDP(18) + MDP(7);
t33 = cos(pkin(7));
t30 = -pkin(1) * t33 - pkin(2);
t27 = -pkin(6) + t30;
t40 = -qJ(5) + t27;
t39 = -MDP(14) - MDP(16);
t32 = sin(pkin(7));
t28 = pkin(1) * t32 + qJ(3);
t38 = MDP(15) + t42;
t37 = MDP(13) + t38;
t24 = pkin(4) * t34 + t28;
t23 = t40 * t35;
t22 = t40 * t34;
t21 = t22 * t34 + t23 * t35;
t1 = [MDP(1) + (t28 ^ 2 + t30 ^ 2) * MDP(7) + t31 * MDP(8) + (t22 ^ 2 + t23 ^ 2 + t24 ^ 2) * MDP(18) + (t32 ^ 2 + t33 ^ 2) * MDP(4) * pkin(1) ^ 2 - 0.2e1 * t21 * MDP(17) + 0.2e1 * t30 * MDP(5) + 0.2e1 * t28 * MDP(6) + 0.2e1 * (MDP(14) * t28 + MDP(16) * t24) * t35 + 0.2e1 * (MDP(13) * t28 + MDP(15) * t24 - MDP(9) * t35) * t34; (t22 * t35 - t23 * t34) * MDP(18); MDP(4) + t41; MDP(17) * t26 + MDP(18) * t21 + MDP(7) * t30 + MDP(5); 0; t41; -t22 * MDP(16) + (-MDP(14) * t27 - MDP(11)) * t34 + t38 * t23 + (MDP(13) * t27 - MDP(17) * pkin(4) + MDP(10)) * t35; -t34 * t37 + t35 * t39; t34 * t39 + t35 * t37; MDP(12) + (0.2e1 * MDP(15) + t42) * pkin(4); t34 * MDP(15) + t35 * MDP(16) + MDP(18) * t24; 0; 0; 0; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
