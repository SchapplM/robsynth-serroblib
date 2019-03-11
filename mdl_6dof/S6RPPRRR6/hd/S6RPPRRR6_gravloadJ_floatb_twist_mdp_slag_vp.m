% Calculate Gravitation load on the joints for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPPRRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:31:29
% EndTime: 2019-03-09 02:31:29
% DurationCPUTime: 0.20s
% Computational Cost: add. (124->55), mult. (203->77), div. (0->0), fcn. (190->8), ass. (0->32)
t62 = cos(qJ(4));
t78 = MDP(16) * t62;
t57 = qJ(5) + qJ(6);
t51 = sin(t57);
t52 = cos(t57);
t58 = sin(qJ(5));
t61 = cos(qJ(5));
t77 = MDP(22) * t61 - MDP(23) * t58 + MDP(29) * t52 - MDP(30) * t51 + MDP(15);
t76 = g(3) * t62;
t59 = sin(qJ(4));
t63 = cos(qJ(1));
t75 = t59 * t63;
t60 = sin(qJ(1));
t74 = t60 * t51;
t73 = t60 * t52;
t72 = t60 * t58;
t71 = t60 * t61;
t70 = t61 * t63;
t41 = -t52 * t63 + t59 * t74;
t42 = -t51 * t63 - t59 * t73;
t43 = -t51 * t75 - t73;
t44 = t52 * t75 - t74;
t69 = (-g(1) * t43 + g(2) * t41 + t51 * t76) * MDP(29) + (g(1) * t44 - g(2) * t42 + t52 * t76) * MDP(30);
t68 = t63 * pkin(1) + t60 * qJ(2);
t50 = g(1) * t63 + g(2) * t60;
t49 = g(1) * t60 - g(2) * t63;
t54 = t63 * qJ(2);
t48 = t59 * t70 - t72;
t47 = -t58 * t75 - t71;
t46 = -t58 * t63 - t59 * t71;
t45 = t59 * t72 - t70;
t1 = [(-g(1) * (-pkin(1) * t60 + t54) - g(2) * t68) * MDP(6) + (-g(1) * (t54 + (-pkin(1) - qJ(3)) * t60) - g(2) * (qJ(3) * t63 + t68)) * MDP(9) + (-g(1) * t46 - g(2) * t48) * MDP(22) + (-g(1) * t45 - g(2) * t47) * MDP(23) + (-g(1) * t42 - g(2) * t44) * MDP(29) + (-g(1) * t41 - g(2) * t43) * MDP(30) + (MDP(3) - MDP(5) - MDP(7)) * t50 + (MDP(15) * t59 + MDP(2) - MDP(4) + MDP(8) + t78) * t49; (-MDP(6) - MDP(9)) * t49; -t50 * MDP(9); (t77 * t59 + t78) * g(3) + (MDP(16) * t59 - t77 * t62) * t50; (-g(1) * t47 + g(2) * t45 + t58 * t76) * MDP(22) + (g(1) * t48 - g(2) * t46 + t61 * t76) * MDP(23) + t69; t69;];
taug  = t1;
