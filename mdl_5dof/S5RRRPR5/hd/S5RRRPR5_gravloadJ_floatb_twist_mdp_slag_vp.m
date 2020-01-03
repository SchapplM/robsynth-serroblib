% Calculate Gravitation load on the joints for
% S5RRRPR5
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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:14:17
% EndTime: 2019-12-31 21:14:18
% DurationCPUTime: 0.14s
% Computational Cost: add. (148->45), mult. (176->67), div. (0->0), fcn. (151->10), ass. (0->29)
t56 = qJ(2) + qJ(3);
t51 = pkin(9) + t56;
t48 = sin(t51);
t72 = g(3) * t48;
t57 = sin(qJ(5));
t59 = sin(qJ(1));
t70 = t59 * t57;
t60 = cos(qJ(5));
t69 = t59 * t60;
t62 = cos(qJ(1));
t68 = t62 * t57;
t67 = t62 * t60;
t53 = cos(t56);
t61 = cos(qJ(2));
t66 = t61 * pkin(2) + pkin(3) * t53;
t49 = cos(t51);
t52 = sin(t56);
t64 = g(1) * t62 + g(2) * t59;
t63 = -g(3) * t53 + t52 * t64;
t65 = t63 * MDP(16) + (g(3) * t52 + t53 * t64) * MDP(17) + (MDP(25) * t60 - t57 * MDP(26)) * (-g(3) * t49 + t48 * t64);
t46 = g(1) * t59 - g(2) * t62;
t58 = sin(qJ(2));
t55 = -qJ(4) - pkin(7) - pkin(6);
t44 = pkin(1) + t66;
t43 = t49 * t67 + t70;
t42 = -t49 * t68 + t69;
t41 = -t49 * t69 + t68;
t40 = t49 * t70 + t67;
t1 = [(-g(1) * (-t44 * t59 - t55 * t62) - g(2) * (t62 * t44 - t59 * t55)) * MDP(19) + (-g(1) * t41 - g(2) * t43) * MDP(25) + (-g(1) * t40 - g(2) * t42) * MDP(26) + (MDP(3) - MDP(18)) * t64 + (-t58 * MDP(10) + MDP(16) * t53 - MDP(17) * t52 + t61 * MDP(9) + MDP(2)) * t46; (-g(3) * t61 + t58 * t64) * MDP(9) + (g(3) * t58 + t61 * t64) * MDP(10) + (-g(3) * t66 - t64 * (-pkin(2) * t58 - pkin(3) * t52)) * MDP(19) + t65; MDP(19) * pkin(3) * t63 + t65; -t46 * MDP(19); (-g(1) * t42 + g(2) * t40 + t57 * t72) * MDP(25) + (g(1) * t43 - g(2) * t41 + t60 * t72) * MDP(26);];
taug = t1;
