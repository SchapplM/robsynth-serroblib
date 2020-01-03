% Calculate Gravitation load on the joints for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPPRR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:18
% EndTime: 2019-12-31 18:04:18
% DurationCPUTime: 0.16s
% Computational Cost: add. (114->45), mult. (194->70), div. (0->0), fcn. (193->8), ass. (0->29)
t59 = sin(qJ(1));
t69 = g(1) * t59;
t55 = qJ(4) + qJ(5);
t49 = sin(t55);
t50 = cos(t55);
t56 = sin(pkin(8));
t57 = cos(pkin(8));
t65 = t49 * t57 - t50 * t56;
t37 = t65 * t59;
t64 = t49 * t56 + t50 * t57;
t38 = t64 * t59;
t61 = cos(qJ(1));
t39 = t65 * t61;
t40 = t64 * t61;
t68 = (g(1) * t39 + g(2) * t37 + g(3) * t64) * MDP(24) + (g(1) * t40 + g(2) * t38 - g(3) * t65) * MDP(25);
t67 = t61 * pkin(1) + t59 * qJ(2);
t48 = g(1) * t61 + g(2) * t59;
t47 = -g(2) * t61 + t69;
t66 = pkin(2) * t57 + qJ(3) * t56;
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t63 = t56 * t60 - t57 * t58;
t62 = t56 * t58 + t57 * t60;
t52 = t61 * qJ(2);
t44 = t62 * t61;
t43 = t63 * t61;
t42 = t62 * t59;
t41 = t63 * t59;
t1 = [(-g(1) * (-pkin(1) * t59 + t52) - g(2) * t67) * MDP(7) + (-g(1) * t52 - g(2) * (t61 * t66 + t67) - (-pkin(1) - t66) * t69) * MDP(11) + (g(1) * t42 - g(2) * t44) * MDP(17) + (g(1) * t41 - g(2) * t43) * MDP(18) + (g(1) * t38 - g(2) * t40) * MDP(24) + (-g(1) * t37 + g(2) * t39) * MDP(25) + (MDP(3) - MDP(6) - MDP(9)) * t48 + (MDP(2) + (MDP(4) + MDP(8)) * t57 + (-MDP(5) + MDP(10)) * t56) * t47; (-MDP(11) - MDP(7)) * t47; (g(3) * t57 - t48 * t56) * MDP(11); (-g(1) * t43 - g(2) * t41 + g(3) * t62) * MDP(17) + (g(1) * t44 + g(2) * t42 + g(3) * t63) * MDP(18) + t68; t68;];
taug = t1;
