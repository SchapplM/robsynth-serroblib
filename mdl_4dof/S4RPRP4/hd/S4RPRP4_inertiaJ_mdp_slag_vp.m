% Calculate joint inertia matrix for
% S4RPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPRP4_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:54
% EndTime: 2019-12-31 16:43:54
% DurationCPUTime: 0.06s
% Computational Cost: add. (58->33), mult. (101->51), div. (0->0), fcn. (56->4), ass. (0->13)
t27 = sin(qJ(3));
t23 = t27 ^ 2;
t28 = cos(qJ(3));
t33 = t28 ^ 2 + t23;
t26 = cos(pkin(6));
t22 = -t26 * pkin(1) - pkin(2);
t32 = -pkin(3) * MDP(15) - MDP(12);
t31 = MDP(10) - t32;
t30 = MDP(15) * qJ(4) - MDP(11) + MDP(14);
t25 = sin(pkin(6));
t21 = t25 * pkin(1) + pkin(5);
t19 = -t28 * pkin(3) - t27 * qJ(4) + t22;
t1 = [MDP(1) + t23 * MDP(5) + (t33 * t21 ^ 2 + t19 ^ 2) * MDP(15) + (t25 ^ 2 + t26 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * t33 * MDP(13) * t21 + 0.2e1 * (-t22 * MDP(10) - t19 * MDP(12)) * t28 + 0.2e1 * (t22 * MDP(11) - t19 * MDP(14) + t28 * MDP(6)) * t27; 0; t33 * MDP(15) + MDP(4); t27 * MDP(7) + t28 * MDP(8) + (-t27 * pkin(3) + t28 * qJ(4)) * MDP(13) + (-t31 * t27 + t30 * t28) * t21; t30 * t27 + t31 * t28; MDP(9) + 0.2e1 * pkin(3) * MDP(12) + 0.2e1 * qJ(4) * MDP(14) + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(15); (MDP(15) * t21 + MDP(13)) * t27; -t28 * MDP(15); t32; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
