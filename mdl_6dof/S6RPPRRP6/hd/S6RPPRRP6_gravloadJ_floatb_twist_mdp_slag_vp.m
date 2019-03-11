% Calculate Gravitation load on the joints for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:11:18
% EndTime: 2019-03-09 02:11:19
% DurationCPUTime: 0.39s
% Computational Cost: add. (136->67), mult. (298->85), div. (0->0), fcn. (277->6), ass. (0->30)
t96 = MDP(16) - MDP(25);
t93 = MDP(22) + MDP(24);
t92 = MDP(23) - MDP(26);
t67 = sin(qJ(4));
t70 = cos(qJ(4));
t85 = t70 * pkin(8);
t95 = t67 * pkin(4) - t85;
t68 = sin(qJ(1));
t71 = cos(qJ(1));
t94 = -g(1) * t71 - g(2) * t68;
t87 = g(3) * t70;
t66 = sin(qJ(5));
t84 = t68 * t66;
t69 = cos(qJ(5));
t83 = t68 * t69;
t82 = t71 * t66;
t81 = t71 * t69;
t80 = -pkin(1) - qJ(3);
t79 = t71 * pkin(1) + t68 * qJ(2);
t78 = -MDP(27) - MDP(9);
t77 = t71 * qJ(3) + t79;
t57 = g(1) * t68 - g(2) * t71;
t75 = pkin(5) * t69 + qJ(6) * t66 + pkin(4);
t52 = t67 * t84 - t81;
t54 = t67 * t82 + t83;
t47 = g(1) * t54 + g(2) * t52 + t66 * t87;
t63 = t71 * qJ(2);
t55 = t67 * t81 - t84;
t53 = t67 * t83 + t82;
t1 = [(-g(1) * (-t68 * pkin(1) + t63) - g(2) * t79) * MDP(6) + (-g(1) * (t68 * t80 + t63) - g(2) * t77) * MDP(9) + (-g(1) * (-t53 * pkin(5) - t71 * pkin(7) - t52 * qJ(6) + t63) - g(2) * (t55 * pkin(5) + t54 * qJ(6) + t95 * t71 + t77) + (-g(1) * (t80 - t95) + g(2) * pkin(7)) * t68) * MDP(27) - t92 * (g(1) * t52 - g(2) * t54) + t93 * (g(1) * t53 - g(2) * t55) - (MDP(3) - MDP(5) - MDP(7)) * t94 + (t67 * MDP(15) + t96 * t70 + MDP(2) - MDP(4) + MDP(8)) * t57; (-MDP(6) + t78) * t57; -t78 * t94; (-g(3) * (-t67 * t75 + t85) + t94 * (pkin(8) * t67 + t70 * t75)) * MDP(27) + t96 * (-t67 * t94 + t87) + (t92 * t66 - t93 * t69 - MDP(15)) * (-g(3) * t67 - t94 * t70); (-g(1) * (-t54 * pkin(5) + t55 * qJ(6)) - g(2) * (-t52 * pkin(5) + t53 * qJ(6)) - (-pkin(5) * t66 + qJ(6) * t69) * t87) * MDP(27) + t92 * (g(1) * t55 + g(2) * t53 + t69 * t87) + t93 * t47; -t47 * MDP(27);];
taug  = t1;
