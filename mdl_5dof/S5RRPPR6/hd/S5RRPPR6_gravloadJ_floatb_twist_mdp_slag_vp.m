% Calculate Gravitation load on the joints for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:06
% EndTime: 2019-12-31 19:33:07
% DurationCPUTime: 0.29s
% Computational Cost: add. (161->61), mult. (202->88), div. (0->0), fcn. (178->10), ass. (0->34)
t61 = sin(qJ(1));
t63 = cos(qJ(1));
t46 = g(1) * t63 + g(2) * t61;
t56 = qJ(2) + pkin(8);
t50 = sin(t56);
t52 = cos(t56);
t66 = -g(3) * t52 + t46 * t50;
t78 = g(3) * t50;
t55 = pkin(9) + qJ(5);
t49 = sin(t55);
t76 = t61 * t49;
t51 = cos(t55);
t75 = t61 * t51;
t57 = sin(pkin(9));
t74 = t61 * t57;
t58 = cos(pkin(9));
t73 = t61 * t58;
t72 = t63 * t49;
t71 = t63 * t51;
t70 = t63 * t57;
t69 = t63 * t58;
t45 = g(1) * t61 - g(2) * t63;
t68 = t52 * pkin(3) + t50 * qJ(4);
t62 = cos(qJ(2));
t60 = sin(qJ(2));
t59 = -qJ(3) - pkin(6);
t53 = t62 * pkin(2);
t48 = t53 + pkin(1);
t47 = t63 * t48;
t44 = t52 * t71 + t76;
t43 = -t52 * t72 + t75;
t42 = -t52 * t75 + t72;
t41 = t52 * t76 + t71;
t1 = [(-g(1) * (-t61 * t48 - t63 * t59) - g(2) * (-t61 * t59 + t47)) * MDP(12) + (-g(1) * (-t52 * t73 + t70) - g(2) * (t52 * t69 + t74)) * MDP(13) + (-g(1) * (t52 * t74 + t69) - g(2) * (-t52 * t70 + t73)) * MDP(14) + (-g(2) * t47 + (g(1) * t59 - g(2) * t68) * t63 + (-g(1) * (-t48 - t68) + g(2) * t59) * t61) * MDP(16) + (-g(1) * t42 - g(2) * t44) * MDP(22) + (-g(1) * t41 - g(2) * t43) * MDP(23) + (MDP(3) - MDP(11)) * t46 + (-t60 * MDP(10) + t50 * MDP(15) + t62 * MDP(9) + MDP(2)) * t45; (g(3) * t60 + t46 * t62) * MDP(10) + (-t46 * t52 - t78) * MDP(15) + (-g(3) * (t53 + t68) + t46 * (pkin(2) * t60 + pkin(3) * t50 - qJ(4) * t52)) * MDP(16) + (pkin(2) * MDP(12) + MDP(9)) * (-g(3) * t62 + t46 * t60) + (MDP(13) * t58 - MDP(14) * t57 + MDP(22) * t51 - MDP(23) * t49) * t66; (-MDP(12) - MDP(16)) * t45; -t66 * MDP(16); (-g(1) * t43 + g(2) * t41 + t49 * t78) * MDP(22) + (g(1) * t44 - g(2) * t42 + t51 * t78) * MDP(23);];
taug = t1;
