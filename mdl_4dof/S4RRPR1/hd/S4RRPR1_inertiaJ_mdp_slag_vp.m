% Calculate joint inertia matrix for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RRPR1_inertiaJ_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:35:01
% EndTime: 2019-03-08 18:35:01
% DurationCPUTime: 0.05s
% Computational Cost: add. (82->34), mult. (143->45), div. (0->0), fcn. (108->6), ass. (0->24)
t39 = sin(qJ(2));
t51 = pkin(1) * t39;
t36 = sin(pkin(7));
t50 = pkin(2) * t36;
t37 = cos(pkin(7));
t34 = t37 * pkin(2) + pkin(3);
t38 = sin(qJ(4));
t49 = t38 * t34;
t41 = cos(qJ(2));
t35 = t41 * pkin(1) + pkin(2);
t28 = t37 * t35 - t36 * t51;
t27 = pkin(3) + t28;
t30 = t36 * t35 + t37 * t51;
t40 = cos(qJ(4));
t24 = t40 * t27 - t38 * t30;
t48 = t24 * MDP(9);
t33 = t40 * t34;
t47 = (-t38 * t50 + t33) * MDP(9);
t25 = -t38 * t27 - t40 * t30;
t46 = t25 * MDP(10);
t45 = (-t40 * t50 - t49) * MDP(10);
t44 = MDP(4) + MDP(8);
t43 = (t41 * MDP(5) - t39 * MDP(6)) * pkin(1);
t1 = [MDP(1) + (t28 ^ 2 + t30 ^ 2) * MDP(7) + 0.2e1 * t43 + 0.2e1 * t46 + 0.2e1 * t48 + t44; (t24 + t33) * MDP(9) + (t25 - t49) * MDP(10) + (t28 * t37 * MDP(7) + (-t40 * MDP(10) + t30 * MDP(7) - t38 * MDP(9)) * t36) * pkin(2) + t43 + t44; (t36 ^ 2 + t37 ^ 2) * MDP(7) * pkin(2) ^ 2 + 0.2e1 * t45 + 0.2e1 * t47 + t44; 0; 0; MDP(7); MDP(8) + t46 + t48; MDP(8) + t45 + t47; 0; MDP(8);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
