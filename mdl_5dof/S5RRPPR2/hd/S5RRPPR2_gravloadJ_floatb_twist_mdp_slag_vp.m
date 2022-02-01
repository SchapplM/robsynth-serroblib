% Calculate Gravitation load on the joints for
% S5RRPPR2
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
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:52
% EndTime: 2022-01-20 10:05:52
% DurationCPUTime: 0.12s
% Computational Cost: add. (165->43), mult. (129->65), div. (0->0), fcn. (110->10), ass. (0->28)
t56 = qJ(1) + qJ(2);
t53 = sin(t56);
t72 = pkin(2) * t53;
t71 = g(3) * sin(pkin(9));
t60 = sin(qJ(1));
t70 = t60 * pkin(1);
t58 = cos(pkin(9));
t59 = sin(qJ(5));
t69 = t58 * t59;
t61 = cos(qJ(5));
t68 = t58 * t61;
t52 = pkin(8) + t56;
t49 = sin(t52);
t50 = cos(t52);
t54 = cos(t56);
t51 = pkin(2) * t54;
t67 = t50 * pkin(3) + t49 * qJ(4) + t51;
t66 = g(1) * t49 - g(2) * t50;
t65 = g(1) * t53 - g(2) * t54;
t39 = t49 * t69 + t50 * t61;
t40 = -t49 * t68 + t50 * t59;
t41 = t49 * t61 - t50 * t69;
t42 = t49 * t59 + t50 * t68;
t64 = (-g(1) * t39 - g(2) * t41) * MDP(17) + (-g(1) * t40 - g(2) * t42) * MDP(16) + t66 * t58 * MDP(8) + (-g(1) * t50 - g(2) * t49) * MDP(9) + t65 * MDP(5) + (g(1) * t54 + g(2) * t53) * MDP(6);
t63 = -t49 * pkin(3) + t50 * qJ(4) - t72;
t62 = cos(qJ(1));
t55 = t62 * pkin(1);
t1 = [(g(1) * t60 - g(2) * t62) * MDP(2) + (g(1) * t62 + g(2) * t60) * MDP(3) + (-g(1) * (-t70 - t72) - g(2) * (t51 + t55)) * MDP(7) + (-g(1) * (t63 - t70) - g(2) * (t55 + t67)) * MDP(10) + t64; t65 * pkin(2) * MDP(7) + (-g(1) * t63 - g(2) * t67) * MDP(10) + t64; (-MDP(10) - MDP(7)) * g(3); -t66 * MDP(10); (-g(1) * t41 + g(2) * t39 + t59 * t71) * MDP(16) + (g(1) * t42 - g(2) * t40 + t61 * t71) * MDP(17);];
taug = t1;
