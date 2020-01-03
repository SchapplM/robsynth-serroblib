% Calculate joint inertia matrix for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RRRP2_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:10
% EndTime: 2019-12-31 17:13:10
% DurationCPUTime: 0.08s
% Computational Cost: add. (89->46), mult. (154->61), div. (0->0), fcn. (99->4), ass. (0->26)
t53 = 2 * MDP(14);
t43 = cos(qJ(2));
t52 = t43 * pkin(1);
t34 = -pkin(2) - t52;
t51 = pkin(2) - t34;
t40 = sin(qJ(3));
t42 = cos(qJ(3));
t50 = t42 * MDP(10) + t40 * MDP(9);
t49 = t40 * MDP(14);
t48 = t42 * MDP(12);
t47 = MDP(4) + (MDP(7) * t40 + 0.2e1 * MDP(8) * t42) * t40;
t35 = -t42 * pkin(3) - pkin(2);
t46 = -t40 * MDP(13) + t48;
t45 = -MDP(12) * t40 - MDP(13) * t42;
t41 = sin(qJ(2));
t44 = (t43 * MDP(5) - t41 * MDP(6)) * pkin(1);
t39 = t42 * qJ(4);
t33 = t41 * pkin(1) + pkin(6);
t31 = t42 * pkin(6) + t39;
t30 = (-qJ(4) - pkin(6)) * t40;
t29 = t35 - t52;
t28 = t31 * t42;
t27 = t42 * t33 + t39;
t26 = (-qJ(4) - t33) * t40;
t25 = t27 * t42;
t1 = [MDP(1) + (-t26 * t40 + t25) * t53 + (t26 ^ 2 + t27 ^ 2 + t29 ^ 2) * MDP(15) + t47 - 0.2e1 * t46 * t34 + 0.2e1 * t44; (t25 + t28) * MDP(14) + (t26 * t30 + t27 * t31 + t29 * t35) * MDP(15) + t51 * t48 + t44 + (-t51 * MDP(13) + (-t26 - t30) * MDP(14)) * t40 + t47; (-t30 * t40 + t28) * t53 + (t30 ^ 2 + t31 ^ 2 + t35 ^ 2) * MDP(15) + 0.2e1 * t46 * pkin(2) + t47; t45 * t33 + (t26 * MDP(15) - t49) * pkin(3) + t50; t45 * pkin(6) + (MDP(15) * t30 - t49) * pkin(3) + t50; pkin(3) ^ 2 * MDP(15) + MDP(11); t29 * MDP(15); t35 * MDP(15); 0; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
