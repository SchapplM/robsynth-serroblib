% Calculate joint inertia matrix for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRPRR9_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:45
% EndTime: 2019-12-31 17:39:46
% DurationCPUTime: 0.06s
% Computational Cost: add. (66->33), mult. (99->39), div. (0->0), fcn. (63->4), ass. (0->21)
t26 = sin(qJ(4));
t28 = cos(qJ(4));
t25 = sin(qJ(5));
t27 = cos(qJ(5));
t36 = t27 * MDP(16);
t33 = -t25 * MDP(17) + t36;
t43 = (MDP(9) + t33) * t28 - t26 * MDP(10);
t29 = -pkin(2) - pkin(3);
t20 = t26 * qJ(3) - t28 * t29;
t18 = pkin(4) + t20;
t42 = pkin(4) + t18;
t41 = t20 * MDP(9);
t40 = t25 ^ 2 * MDP(11) + MDP(8);
t21 = t28 * qJ(3) + t26 * t29;
t39 = t21 * MDP(10);
t37 = t27 * MDP(12);
t35 = 0.2e1 * t25 * t37 + t40;
t34 = -t25 * MDP(13) - t27 * MDP(14);
t32 = -MDP(16) * t25 - MDP(17) * t27;
t30 = 0.2e1 * t33;
t1 = [MDP(1) + MDP(7); 0; MDP(2) + (2 * pkin(2) * MDP(5)) + 0.2e1 * qJ(3) * MDP(6) + ((pkin(2) ^ 2) + qJ(3) ^ 2) * MDP(7) + t18 * t30 + 0.2e1 * t39 + 0.2e1 * t41 + t35; 0; -pkin(2) * MDP(7) - MDP(5) - t43; MDP(7); 0; -t39 - t41 - t42 * t36 + (t42 * MDP(17) - 0.2e1 * t37) * t25 - t40; t43; pkin(4) * t30 + t35; -t33; t32 * (-pkin(7) + t21) + t34; t32 * t26; t32 * pkin(7) - t34; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
