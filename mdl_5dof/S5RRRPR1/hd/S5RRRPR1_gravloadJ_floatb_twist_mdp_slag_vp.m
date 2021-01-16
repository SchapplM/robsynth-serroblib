% Calculate Gravitation load on the joints for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:50
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:48:58
% EndTime: 2021-01-15 22:48:59
% DurationCPUTime: 0.13s
% Computational Cost: add. (210->42), mult. (180->55), div. (0->0), fcn. (141->10), ass. (0->22)
t49 = qJ(2) + qJ(3);
t44 = pkin(9) + t49;
t43 = qJ(5) + t44;
t38 = sin(t43);
t39 = cos(t43);
t51 = sin(qJ(1));
t53 = cos(qJ(1));
t56 = g(1) * t53 + g(2) * t51;
t58 = (-g(3) * t39 + t56 * t38) * MDP(27) + (g(3) * t38 + t56 * t39) * MDP(28);
t46 = cos(t49);
t52 = cos(qJ(2));
t57 = t52 * pkin(2) + pkin(3) * t46;
t36 = g(1) * t51 - g(2) * t53;
t40 = sin(t44);
t41 = cos(t44);
t45 = sin(t49);
t54 = -g(3) * t46 + t56 * t45;
t55 = t54 * MDP(16) + (g(3) * t45 + t56 * t46) * MDP(17) + (-g(3) * t41 + t56 * t40) * MDP(18) + (g(3) * t40 + t56 * t41) * MDP(19) + t58;
t50 = sin(qJ(2));
t48 = -qJ(4) - pkin(7) - pkin(6);
t34 = pkin(1) + t57;
t1 = [(-g(1) * (-t51 * t34 - t53 * t48) - g(2) * (t53 * t34 - t51 * t48)) * MDP(21) + (MDP(3) - MDP(20)) * t56 + (-t50 * MDP(10) + MDP(16) * t46 - MDP(17) * t45 + MDP(18) * t41 - MDP(19) * t40 + MDP(27) * t39 - MDP(28) * t38 + t52 * MDP(9) + MDP(2)) * t36; (-g(3) * t52 + t56 * t50) * MDP(9) + (g(3) * t50 + t56 * t52) * MDP(10) + (-g(3) * t57 - t56 * (-t50 * pkin(2) - pkin(3) * t45)) * MDP(21) + t55; t54 * MDP(21) * pkin(3) + t55; -t36 * MDP(21); t58;];
taug = t1;
