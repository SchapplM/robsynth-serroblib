% Calculate Gravitation load on the joints for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:36
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:35:12
% EndTime: 2021-01-15 21:35:13
% DurationCPUTime: 0.16s
% Computational Cost: add. (166->46), mult. (184->64), div. (0->0), fcn. (159->10), ass. (0->27)
t48 = qJ(2) + pkin(9);
t47 = qJ(4) + t48;
t42 = sin(t47);
t65 = g(3) * t42;
t50 = sin(qJ(5));
t52 = sin(qJ(1));
t63 = t52 * t50;
t53 = cos(qJ(5));
t62 = t52 * t53;
t55 = cos(qJ(1));
t61 = t55 * t50;
t60 = t55 * t53;
t43 = cos(t47);
t58 = g(1) * t55 + g(2) * t52;
t59 = (t58 * t43 + t65) * MDP(21) + (t53 * MDP(27) - t50 * MDP(28) + MDP(20)) * (-g(3) * t43 + t58 * t42);
t40 = g(1) * t52 - g(2) * t55;
t54 = cos(qJ(2));
t51 = sin(qJ(2));
t49 = -qJ(3) - pkin(6);
t46 = cos(t48);
t45 = sin(t48);
t44 = t54 * pkin(2) + pkin(1);
t39 = t43 * t60 + t63;
t38 = -t43 * t61 + t62;
t37 = -t43 * t62 + t61;
t36 = t43 * t63 + t60;
t1 = [(-g(1) * (-t52 * t44 - t49 * t55) - g(2) * (t55 * t44 - t52 * t49)) * MDP(14) + (-g(1) * t37 - g(2) * t39) * MDP(27) + (-g(1) * t36 - g(2) * t38) * MDP(28) + (MDP(3) - MDP(13)) * t58 + (-t51 * MDP(10) + t46 * MDP(11) - t45 * MDP(12) + MDP(20) * t43 - MDP(21) * t42 + MDP(9) * t54 + MDP(2)) * t40; (g(3) * t51 + t58 * t54) * MDP(10) + (-g(3) * t46 + t58 * t45) * MDP(11) + (g(3) * t45 + t58 * t46) * MDP(12) + t59 + (MDP(14) * pkin(2) + MDP(9)) * (-g(3) * t54 + t58 * t51); -t40 * MDP(14); t59; (-g(1) * t38 + g(2) * t36 + t50 * t65) * MDP(27) + (g(1) * t39 - g(2) * t37 + t53 * t65) * MDP(28);];
taug = t1;
