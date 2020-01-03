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
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5RPPRP3_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:02
% EndTime: 2019-12-31 17:51:03
% DurationCPUTime: 0.09s
% Computational Cost: add. (78->36), mult. (113->53), div. (0->0), fcn. (70->4), ass. (0->20)
t32 = cos(pkin(7));
t29 = -t32 * pkin(1) - pkin(2);
t41 = t29 * MDP(7);
t40 = MDP(16) * pkin(4);
t34 = cos(qJ(4));
t30 = t34 ^ 2;
t33 = sin(qJ(4));
t25 = -t33 ^ 2 - t30;
t39 = -t25 * MDP(16) + MDP(7);
t26 = -pkin(6) + t29;
t38 = -qJ(5) + t26;
t37 = t34 * MDP(14);
t31 = sin(pkin(7));
t27 = t31 * pkin(1) + qJ(3);
t36 = MDP(13) + t40;
t23 = t33 * pkin(4) + t27;
t22 = t38 * t34;
t21 = t38 * t33;
t20 = t21 * t33 + t22 * t34;
t1 = [MDP(1) + t30 * MDP(8) - 0.2e1 * t34 * t33 * MDP(9) + (t21 ^ 2 + t22 ^ 2 + t23 ^ 2) * MDP(16) + (t31 ^ 2 + t32 ^ 2) * MDP(4) * pkin(1) ^ 2 - 0.2e1 * t20 * MDP(15) + ((2 * MDP(5)) + t41) * t29 + (0.2e1 * t33 * MDP(13) + MDP(7) * t27 + (2 * MDP(6)) + 0.2e1 * t37) * t27; (t21 * t34 - t22 * t33) * MDP(16); MDP(4) + t39; t25 * MDP(15) + t20 * MDP(16) + MDP(5) + t41; 0; t39; t22 * t40 + (-MDP(14) * t26 - MDP(11)) * t33 + (MDP(13) * t26 - MDP(15) * pkin(4) + MDP(10)) * t34; -t36 * t33 - t37; -t33 * MDP(14) + t36 * t34; pkin(4) ^ 2 * MDP(16) + MDP(12); t23 * MDP(16); 0; 0; 0; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
