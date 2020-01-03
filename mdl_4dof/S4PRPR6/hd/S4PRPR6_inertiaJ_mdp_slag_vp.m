% Calculate joint inertia matrix for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRPR6_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:41
% EndTime: 2019-12-31 16:24:41
% DurationCPUTime: 0.06s
% Computational Cost: add. (54->31), mult. (119->51), div. (0->0), fcn. (101->6), ass. (0->20)
t31 = cos(pkin(7));
t45 = -0.2e1 * t31 * pkin(3) - (2 * pkin(2));
t44 = pkin(2) * MDP(8);
t43 = pkin(5) + qJ(3);
t30 = sin(pkin(7));
t42 = t30 ^ 2 + t31 ^ 2;
t41 = t30 * MDP(6);
t40 = t31 * MDP(5);
t32 = sin(qJ(4));
t34 = cos(qJ(4));
t21 = t32 * t30 - t34 * t31;
t39 = t21 * MDP(14);
t38 = t42 * qJ(3);
t22 = t34 * t30 + t32 * t31;
t37 = t22 * MDP(15) + t39 - t40 + t41 - t44;
t35 = cos(qJ(2));
t33 = sin(qJ(2));
t24 = t43 * t31;
t23 = t43 * t30;
t1 = [MDP(1) + (t42 * t33 ^ 2 + t35 ^ 2) * MDP(8); (t42 * MDP(7) + MDP(8) * t38 - MDP(4)) * t33 + (MDP(3) - t37) * t35; t39 * t45 + MDP(2) + t42 * MDP(8) * qJ(3) ^ 2 + 0.2e1 * MDP(7) * t38 + (0.2e1 * t40 - 0.2e1 * t41 + t44) * pkin(2) + (-0.2e1 * t21 * MDP(10) + MDP(15) * t45 + MDP(9) * t22) * t22; -t35 * MDP(8); t37; MDP(8); (-t22 * MDP(14) + t21 * MDP(15)) * t33; t22 * MDP(11) - t21 * MDP(12) + (-t34 * t23 - t32 * t24) * MDP(14) + (t32 * t23 - t34 * t24) * MDP(15); 0; MDP(13);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
