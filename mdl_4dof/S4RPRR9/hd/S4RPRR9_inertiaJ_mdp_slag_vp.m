% Calculate joint inertia matrix for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S4RPRR9_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:29
% EndTime: 2019-12-31 16:56:29
% DurationCPUTime: 0.13s
% Computational Cost: add. (73->49), mult. (148->69), div. (0->0), fcn. (102->4), ass. (0->22)
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t51 = t34 * MDP(19) + t36 * MDP(20);
t50 = (pkin(1) * MDP(6));
t49 = t34 * t36;
t38 = -pkin(1) - pkin(5);
t48 = t34 * t38;
t47 = t36 * t38;
t43 = MDP(15) * t49;
t42 = MDP(16) * t36 - MDP(17) * t34;
t41 = t36 * MDP(19) - t34 * MDP(20);
t39 = t34 * MDP(16) + t36 * MDP(17) - pkin(6) * t51;
t37 = cos(qJ(3));
t35 = sin(qJ(3));
t33 = t37 ^ 2;
t32 = t36 ^ 2;
t31 = t35 ^ 2;
t30 = t34 ^ 2;
t29 = t35 * pkin(3) - t37 * pkin(6) + qJ(2);
t28 = t34 * t29 + t35 * t47;
t27 = t29 * t36 - t35 * t48;
t1 = [t31 * MDP(18) + MDP(1) + ((-2 * MDP(4) + t50) * pkin(1)) + (0.2e1 * MDP(12) * t35 + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + (MDP(14) * t32 + MDP(7) - 0.2e1 * t43) * t33 + 0.2e1 * (t27 * t35 - t33 * t48) * MDP(19) + 0.2e1 * (-t28 * t35 - t33 * t47) * MDP(20) + 0.2e1 * (qJ(2) * MDP(13) + (-MDP(8) + t42) * t35) * t37; MDP(4) - t50 + t51 * (-t31 - t33); MDP(6); (-t38 * MDP(13) - MDP(10) + t39) * t35 + (MDP(9) + t38 * MDP(12) + MDP(14) * t49 + (-t30 + t32) * MDP(15) + (-pkin(3) * t34 + t47) * MDP(19) + (-pkin(3) * t36 - t48) * MDP(20)) * t37; -t35 * MDP(13) + (MDP(12) + t41) * t37; t30 * MDP(14) + 0.2e1 * pkin(3) * t41 + MDP(11) + 0.2e1 * t43; t35 * MDP(18) + t27 * MDP(19) - t28 * MDP(20) + t37 * t42; -t51 * t35; t39; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
