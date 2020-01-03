% Calculate Gravitation load on the joints for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:06
% EndTime: 2019-12-31 17:17:06
% DurationCPUTime: 0.12s
% Computational Cost: add. (125->38), mult. (155->50), div. (0->0), fcn. (124->6), ass. (0->22)
t47 = cos(qJ(2));
t44 = qJ(2) + qJ(3);
t41 = sin(t44);
t42 = cos(t44);
t55 = t42 * pkin(3) + t41 * qJ(4);
t59 = t47 * pkin(2) + t55;
t58 = MDP(16) + MDP(18);
t57 = MDP(17) - MDP(20);
t56 = pkin(3) * t41;
t54 = qJ(4) * t42;
t46 = sin(qJ(1));
t48 = cos(qJ(1));
t35 = g(1) * t48 + g(2) * t46;
t31 = -g(3) * t42 + t35 * t41;
t53 = t57 * (g(3) * t41 + t35 * t42) + t58 * t31;
t45 = sin(qJ(2));
t52 = -pkin(2) * t45 - t56;
t50 = pkin(1) + t59;
t49 = -pkin(6) - pkin(5);
t37 = t48 * t54;
t36 = t46 * t54;
t1 = [((g(1) * t49 - g(2) * t50) * t48 + (g(1) * t50 + g(2) * t49) * t46) * MDP(21) + (MDP(3) - MDP(19)) * t35 + (-t45 * MDP(10) + t47 * MDP(9) - t57 * t41 + t58 * t42 + MDP(2)) * (g(1) * t46 - g(2) * t48); (-g(3) * t47 + t35 * t45) * MDP(9) + (g(3) * t45 + t35 * t47) * MDP(10) + (-g(1) * (t52 * t48 + t37) - g(2) * (t52 * t46 + t36) - g(3) * t59) * MDP(21) + t53; (-g(1) * (-t48 * t56 + t37) - g(2) * (-t46 * t56 + t36) - g(3) * t55) * MDP(21) + t53; -t31 * MDP(21);];
taug = t1;
