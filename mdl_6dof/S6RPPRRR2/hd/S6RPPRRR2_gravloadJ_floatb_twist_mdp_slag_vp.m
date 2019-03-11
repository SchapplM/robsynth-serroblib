% Calculate Gravitation load on the joints for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPPRRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:21:22
% EndTime: 2019-03-09 02:21:23
% DurationCPUTime: 0.18s
% Computational Cost: add. (234->56), mult. (198->82), div. (0->0), fcn. (186->12), ass. (0->36)
t57 = pkin(11) + qJ(4);
t51 = sin(t57);
t84 = t51 * MDP(15);
t59 = qJ(5) + qJ(6);
t55 = sin(t59);
t56 = cos(t59);
t62 = sin(qJ(5));
t64 = cos(qJ(5));
t83 = t64 * MDP(21) - t62 * MDP(22) + t56 * MDP(28) - t55 * MDP(29) + MDP(14);
t82 = g(3) * t51;
t58 = qJ(1) + pkin(10);
t52 = sin(t58);
t81 = t52 * t55;
t80 = t52 * t56;
t79 = t52 * t62;
t78 = t52 * t64;
t54 = cos(t58);
t77 = t54 * t55;
t76 = t54 * t56;
t75 = t54 * t62;
t74 = t54 * t64;
t53 = cos(t57);
t43 = t53 * t81 + t76;
t44 = -t53 * t80 + t77;
t45 = -t53 * t77 + t80;
t46 = t53 * t76 + t81;
t73 = (-g(1) * t45 + g(2) * t43 + t55 * t82) * MDP(28) + (g(1) * t46 - g(2) * t44 + t56 * t82) * MDP(29);
t68 = g(1) * t54 + g(2) * t52;
t67 = g(1) * t52 - g(2) * t54;
t65 = cos(qJ(1));
t63 = sin(qJ(1));
t50 = t53 * t74 + t79;
t49 = -t53 * t75 + t78;
t48 = -t53 * t78 + t75;
t47 = t53 * t79 + t74;
t1 = [(g(1) * t65 + g(2) * t63) * MDP(3) - t68 * MDP(7) + (-g(1) * (-t63 * pkin(1) - t52 * pkin(2) + t54 * qJ(3)) - g(2) * (t65 * pkin(1) + t54 * pkin(2) + t52 * qJ(3))) * MDP(8) + (-g(1) * t48 - g(2) * t50) * MDP(21) + (-g(1) * t47 - g(2) * t49) * MDP(22) + (-g(1) * t44 - g(2) * t46) * MDP(28) + (-g(1) * t43 - g(2) * t45) * MDP(29) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t63 - g(2) * t65) + (t53 * MDP(14) - t84 + MDP(5) * cos(pkin(11)) - MDP(6) * sin(pkin(11))) * t67; (-MDP(4) - MDP(8)) * g(3); -t67 * MDP(8); (-t83 * t53 + t84) * g(3) + (MDP(15) * t53 + t83 * t51) * t68; (-g(1) * t49 + g(2) * t47 + t62 * t82) * MDP(21) + (g(1) * t50 - g(2) * t48 + t64 * t82) * MDP(22) + t73; t73;];
taug  = t1;
