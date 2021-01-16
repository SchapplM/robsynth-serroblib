% Calculate Gravitation load on the joints for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:52:43
% EndTime: 2021-01-15 23:52:44
% DurationCPUTime: 0.14s
% Computational Cost: add. (241->41), mult. (209->52), div. (0->0), fcn. (164->8), ass. (0->23)
t65 = MDP(23) + MDP(25);
t64 = MDP(24) + MDP(26);
t53 = qJ(2) + qJ(3);
t51 = qJ(4) + t53;
t47 = cos(t51);
t49 = cos(t53);
t63 = pkin(3) * t49 + pkin(4) * t47;
t56 = cos(qJ(2));
t62 = t56 * pkin(2) + t63;
t46 = sin(t51);
t55 = sin(qJ(1));
t57 = cos(qJ(1));
t59 = g(1) * t57 + g(2) * t55;
t33 = -g(3) * t47 + t59 * t46;
t61 = t64 * (g(3) * t46 + t59 * t47) + t65 * t33;
t48 = sin(t53);
t60 = -pkin(3) * t48 - pkin(4) * t46;
t42 = g(1) * t55 - g(2) * t57;
t58 = (-g(3) * t49 + t59 * t48) * MDP(16) + (g(3) * t48 + t59 * t49) * MDP(17) + t61;
t54 = sin(qJ(2));
t50 = -qJ(5) - pkin(8) - pkin(7) - pkin(6);
t39 = pkin(1) + t62;
t1 = [(-g(1) * (-t55 * t39 - t57 * t50) - g(2) * (t57 * t39 - t55 * t50)) * MDP(28) + (MDP(3) - MDP(27)) * t59 + (-t54 * MDP(10) + MDP(16) * t49 - MDP(17) * t48 + t56 * MDP(9) - t64 * t46 + t65 * t47 + MDP(2)) * t42; (-g(3) * t56 + t59 * t54) * MDP(9) + (g(3) * t54 + t59 * t56) * MDP(10) + (-g(3) * t62 - t59 * (-t54 * pkin(2) + t60)) * MDP(28) + t58; (-g(3) * t63 - t59 * t60) * MDP(28) + t58; t33 * MDP(28) * pkin(4) + t61; -t42 * MDP(28);];
taug = t1;
