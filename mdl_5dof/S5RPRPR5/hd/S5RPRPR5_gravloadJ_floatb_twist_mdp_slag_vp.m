% Calculate Gravitation load on the joints for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:57:12
% EndTime: 2019-12-05 17:57:13
% DurationCPUTime: 0.19s
% Computational Cost: add. (139->52), mult. (182->77), div. (0->0), fcn. (173->8), ass. (0->33)
t62 = sin(qJ(3));
t80 = pkin(3) * t62;
t59 = sin(pkin(8));
t79 = g(1) * t59;
t57 = qJ(3) + pkin(9) + qJ(5);
t54 = sin(t57);
t63 = sin(qJ(1));
t76 = t63 * t54;
t55 = cos(t57);
t75 = t63 * t55;
t74 = t63 * t62;
t64 = cos(qJ(3));
t73 = t63 * t64;
t65 = cos(qJ(1));
t72 = t65 * t54;
t71 = t65 * t55;
t70 = t65 * t62;
t69 = t65 * t64;
t60 = cos(pkin(8));
t43 = t60 * t76 + t71;
t44 = t60 * t75 - t72;
t45 = t60 * t72 - t75;
t46 = -t60 * t71 - t76;
t68 = (-g(2) * t43 + g(3) * t45 + t54 * t79) * MDP(22) + (-g(2) * t44 - g(3) * t46 + t55 * t79) * MDP(23);
t53 = g(2) * t65 + g(3) * t63;
t52 = g(2) * t63 - g(3) * t65;
t66 = (t64 * pkin(3) + pkin(2)) * t60 - t59 * (-qJ(4) - pkin(6)) + pkin(1);
t49 = t60 * t70 - t73;
t47 = t60 * t74 + t69;
t58 = t65 * qJ(2);
t50 = -t60 * t69 - t74;
t48 = t60 * t73 - t70;
t1 = [(-g(2) * (-t65 * pkin(1) - t63 * qJ(2)) - g(3) * (-t63 * pkin(1) + t58)) * MDP(7) + (-g(2) * t50 + g(3) * t48) * MDP(13) + (-g(2) * t49 - g(3) * t47) * MDP(14) + (-g(3) * t58 + (g(2) * t66 - g(3) * t80) * t65 + (-g(2) * (-qJ(2) - t80) + g(3) * t66) * t63) * MDP(16) + (-g(2) * t46 + g(3) * t44) * MDP(22) + (-g(2) * t45 - g(3) * t43) * MDP(23) + (-MDP(3) + MDP(6)) * t52 + (t60 * MDP(4) + MDP(2) + (-MDP(5) + MDP(15)) * t59) * t53; (-MDP(16) - MDP(7)) * t53; (-g(2) * t48 - g(3) * t50 + t64 * t79) * MDP(14) + t68 + (pkin(3) * MDP(16) + MDP(13)) * (-g(2) * t47 + g(3) * t49 + t62 * t79); (g(1) * t60 + t52 * t59) * MDP(16); t68;];
taug = t1;
