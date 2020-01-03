% Calculate joint inertia matrix for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPRP5_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:02
% EndTime: 2019-12-31 16:45:02
% DurationCPUTime: 0.08s
% Computational Cost: add. (114->46), mult. (209->65), div. (0->0), fcn. (184->4), ass. (0->19)
t45 = cos(qJ(3));
t44 = pkin(1) * MDP(7);
t43 = pkin(5) + qJ(2);
t33 = sin(pkin(6));
t41 = t33 * MDP(5);
t34 = cos(pkin(6));
t40 = t34 * MDP(4);
t39 = MDP(14) - MDP(17);
t30 = -t34 * pkin(2) - pkin(1);
t38 = t43 * t33;
t37 = -pkin(3) * MDP(18) - MDP(15);
t35 = sin(qJ(3));
t28 = t43 * t34;
t27 = t45 * t33 + t35 * t34;
t26 = t35 * t33 - t45 * t34;
t24 = t45 * t28 - t35 * t38;
t23 = t35 * t28 + t45 * t38;
t22 = t26 * pkin(3) - t27 * qJ(4) + t30;
t1 = [MDP(1) + (t22 ^ 2 + t23 ^ 2 + t24 ^ 2) * MDP(18) + (0.2e1 * t40 - 0.2e1 * t41 + t44) * pkin(1) + 0.2e1 * (t30 * MDP(13) + t22 * MDP(15) - MDP(16) * t24) * t26 + (0.2e1 * t30 * MDP(14) + 0.2e1 * MDP(16) * t23 - 0.2e1 * t22 * MDP(17) + MDP(8) * t27 - 0.2e1 * t26 * MDP(9)) * t27 + (MDP(7) * qJ(2) + 2 * MDP(6)) * (t33 ^ 2 + t34 ^ 2) * qJ(2); t22 * MDP(18) - t40 + t41 - t44 + t39 * t27 + (MDP(13) + MDP(15)) * t26; MDP(7) + MDP(18); t27 * MDP(10) - t26 * MDP(11) + (-pkin(3) * t27 - t26 * qJ(4)) * MDP(16) + (MDP(18) * qJ(4) - t39) * t24 + (-MDP(13) + t37) * t23; 0; MDP(12) + 0.2e1 * pkin(3) * MDP(15) + 0.2e1 * qJ(4) * MDP(17) + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(18); t27 * MDP(16) + t23 * MDP(18); 0; t37; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
