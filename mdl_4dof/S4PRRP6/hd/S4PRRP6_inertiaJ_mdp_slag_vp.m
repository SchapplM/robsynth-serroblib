% Calculate joint inertia matrix for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP6_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:43
% EndTime: 2019-12-31 16:30:43
% DurationCPUTime: 0.07s
% Computational Cost: add. (52->35), mult. (115->53), div. (0->0), fcn. (69->4), ass. (0->13)
t25 = sin(qJ(3));
t22 = t25 ^ 2;
t27 = cos(qJ(3));
t35 = t27 ^ 2 + t22;
t34 = pkin(5) * MDP(15);
t33 = -MDP(11) + MDP(14);
t32 = t35 * MDP(13);
t31 = -pkin(3) * MDP(15) - MDP(12);
t30 = (MDP(15) * qJ(4) + t33) * t27 + (-MDP(10) + t31) * t25;
t28 = cos(qJ(2));
t26 = sin(qJ(2));
t21 = -t27 * pkin(3) - t25 * qJ(4) - pkin(2);
t1 = [MDP(1) + (t35 * t26 ^ 2 + t28 ^ 2) * MDP(15); (t35 * t34 - MDP(4) + t32) * t26 + (-t21 * MDP(15) + MDP(3) + (MDP(10) + MDP(12)) * t27 + t33 * t25) * t28; MDP(2) + t22 * MDP(5) + (t35 * pkin(5) ^ 2 + t21 ^ 2) * MDP(15) + 0.2e1 * pkin(5) * t32 + 0.2e1 * (pkin(2) * MDP(10) - t21 * MDP(12)) * t27 + 0.2e1 * (-pkin(2) * MDP(11) - t21 * MDP(14) + t27 * MDP(6)) * t25; t30 * t26; t25 * MDP(7) + t27 * MDP(8) + (-t25 * pkin(3) + t27 * qJ(4)) * MDP(13) + t30 * pkin(5); MDP(9) + 0.2e1 * pkin(3) * MDP(12) + 0.2e1 * qJ(4) * MDP(14) + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(15); t25 * t26 * MDP(15); (MDP(13) + t34) * t25; t31; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
