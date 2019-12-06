% Calculate Gravitation load on the joints for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:31:33
% EndTime: 2019-12-05 17:31:35
% DurationCPUTime: 0.36s
% Computational Cost: add. (122->61), mult. (280->98), div. (0->0), fcn. (306->10), ass. (0->36)
t63 = cos(pkin(7));
t62 = cos(pkin(8));
t67 = cos(qJ(1));
t72 = t67 * t62;
t59 = sin(pkin(8));
t65 = sin(qJ(1));
t75 = t65 * t59;
t52 = t63 * t72 + t75;
t58 = sin(pkin(9));
t61 = cos(pkin(9));
t60 = sin(pkin(7));
t76 = t60 * t67;
t45 = t52 * t61 + t58 * t76;
t73 = t67 * t59;
t74 = t65 * t62;
t51 = t63 * t73 - t74;
t64 = sin(qJ(5));
t66 = cos(qJ(5));
t83 = t45 * t64 - t51 * t66;
t82 = t45 * t66 + t51 * t64;
t81 = g(2) * t67;
t78 = t59 * t60;
t77 = t60 * t65;
t71 = MDP(11) + MDP(15);
t55 = g(3) * t65 + t81;
t54 = g(2) * t65 - g(3) * t67;
t69 = pkin(2) * t63 + qJ(3) * t60 + pkin(1);
t68 = (g(2) * qJ(2) + g(3) * t69) * t65 + t69 * t81;
t57 = t67 * qJ(2);
t50 = -t63 * t74 + t73;
t49 = t63 * t75 + t72;
t48 = t60 * t62 * t61 - t63 * t58;
t44 = t50 * t61 - t58 * t77;
t43 = t44 * t66 - t49 * t64;
t42 = -t44 * t64 - t49 * t66;
t1 = [(-g(2) * (-t67 * pkin(1) - t65 * qJ(2)) - g(3) * (-t65 * pkin(1) + t57)) * MDP(7) + (g(2) * t52 - g(3) * t50) * MDP(8) + (-g(3) * t57 + t68) * MDP(11) + (g(2) * t45 - g(3) * t44) * MDP(12) + (-g(2) * (t52 * t58 - t61 * t76) - g(3) * (-t50 * t58 - t61 * t77)) * MDP(13) + (-g(2) * (-t52 * pkin(3) - t51 * qJ(4)) - g(3) * (t50 * pkin(3) - t49 * qJ(4) + t57) + t68) * MDP(15) + (g(2) * t82 - g(3) * t43) * MDP(21) + (-g(2) * t83 - g(3) * t42) * MDP(22) + (-MDP(9) + MDP(14)) * (g(2) * t51 + g(3) * t49) + (-MDP(3) + MDP(6)) * t54 + (t63 * MDP(4) + MDP(2) + (-MDP(5) + MDP(10)) * t60) * t55; (-MDP(7) - t71) * t55; t71 * (g(1) * t63 + t54 * t60); (-g(1) * t78 + g(2) * t49 - g(3) * t51) * MDP(15); (-g(1) * (-t48 * t64 + t66 * t78) - g(2) * t42 + g(3) * t83) * MDP(21) + (-g(1) * (-t48 * t66 - t64 * t78) + g(2) * t43 + g(3) * t82) * MDP(22);];
taug = t1;
