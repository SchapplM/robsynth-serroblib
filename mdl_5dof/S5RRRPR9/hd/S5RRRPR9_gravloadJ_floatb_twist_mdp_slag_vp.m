% Calculate Gravitation load on the joints for
% S5RRRPR9
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
%   see S5RRRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:24:39
% EndTime: 2019-12-31 21:24:39
% DurationCPUTime: 0.20s
% Computational Cost: add. (174->55), mult. (237->81), div. (0->0), fcn. (223->8), ass. (0->33)
t65 = sin(qJ(1));
t68 = cos(qJ(1));
t74 = g(1) * t68 + g(2) * t65;
t90 = MDP(10) - MDP(18);
t64 = sin(qJ(2));
t67 = cos(qJ(2));
t51 = -g(3) * t67 + t74 * t64;
t84 = g(3) * t64;
t82 = t65 * t67;
t61 = qJ(3) + pkin(9) + qJ(5);
t58 = sin(t61);
t81 = t68 * t58;
t59 = cos(t61);
t80 = t68 * t59;
t63 = sin(qJ(3));
t79 = t68 * t63;
t66 = cos(qJ(3));
t78 = t68 * t66;
t47 = t58 * t82 + t80;
t48 = -t59 * t82 + t81;
t49 = t65 * t59 - t67 * t81;
t50 = t65 * t58 + t67 * t80;
t77 = (-g(1) * t49 + g(2) * t47 + t58 * t84) * MDP(25) + (g(1) * t50 - g(2) * t48 + t59 * t84) * MDP(26);
t75 = pkin(3) * t63 + pkin(6);
t60 = t66 * pkin(3) + pkin(2);
t62 = -qJ(4) - pkin(7);
t72 = t67 * t60 - t64 * t62;
t70 = pkin(1) + t72;
t55 = t65 * t66 - t67 * t79;
t53 = t63 * t82 + t78;
t56 = t65 * t63 + t67 * t78;
t54 = -t66 * t82 + t79;
t1 = [t74 * MDP(3) + (-g(1) * t54 - g(2) * t56) * MDP(16) + (-g(1) * t53 - g(2) * t55) * MDP(17) + ((-g(1) * t75 - g(2) * t70) * t68 + (g(1) * t70 - g(2) * t75) * t65) * MDP(19) + (-g(1) * t48 - g(2) * t50) * MDP(25) + (-g(1) * t47 - g(2) * t49) * MDP(26) + (t67 * MDP(9) - t90 * t64 + MDP(2)) * (g(1) * t65 - g(2) * t68); (-g(3) * t72 + t74 * (t60 * t64 + t62 * t67)) * MDP(19) + t90 * (t74 * t67 + t84) + (MDP(16) * t66 - MDP(17) * t63 + MDP(25) * t59 - MDP(26) * t58 + MDP(9)) * t51; (g(1) * t56 - g(2) * t54 + t66 * t84) * MDP(17) + t77 + (pkin(3) * MDP(19) + MDP(16)) * (-g(1) * t55 + g(2) * t53 + t63 * t84); -t51 * MDP(19); t77;];
taug = t1;
