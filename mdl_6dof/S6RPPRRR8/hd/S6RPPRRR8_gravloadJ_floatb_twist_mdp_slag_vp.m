% Calculate Gravitation load on the joints for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPPRRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:36:20
% EndTime: 2019-03-09 02:36:21
% DurationCPUTime: 0.23s
% Computational Cost: add. (174->57), mult. (210->78), div. (0->0), fcn. (196->10), ass. (0->35)
t62 = pkin(10) + qJ(4);
t55 = cos(t62);
t86 = t55 * MDP(17);
t63 = qJ(5) + qJ(6);
t56 = sin(t63);
t57 = cos(t63);
t66 = sin(qJ(5));
t68 = cos(qJ(5));
t85 = t68 * MDP(23) - t66 * MDP(24) + t57 * MDP(30) - t56 * MDP(31) + MDP(16);
t84 = g(3) * t55;
t67 = sin(qJ(1));
t83 = t67 * t56;
t82 = t67 * t57;
t81 = t67 * t66;
t80 = t67 * t68;
t69 = cos(qJ(1));
t79 = t69 * t56;
t78 = t69 * t57;
t77 = t69 * t66;
t76 = t69 * t68;
t54 = sin(t62);
t44 = -t54 * t83 + t78;
t45 = t54 * t82 + t79;
t46 = t54 * t79 + t82;
t47 = t54 * t78 - t83;
t75 = (-g(1) * t44 - g(2) * t46 + t56 * t84) * MDP(30) + (g(1) * t45 - g(2) * t47 + t57 * t84) * MDP(31);
t74 = t69 * pkin(1) + t67 * qJ(2);
t53 = g(1) * t69 + g(2) * t67;
t52 = g(1) * t67 - g(2) * t69;
t59 = t69 * qJ(2);
t51 = t54 * t76 - t81;
t50 = t54 * t77 + t80;
t49 = t54 * t80 + t77;
t48 = -t54 * t81 + t76;
t1 = [(-g(1) * (-t67 * pkin(1) + t59) - g(2) * t74) * MDP(6) + (-g(1) * (t59 + (-pkin(1) - qJ(3)) * t67) - g(2) * (t69 * qJ(3) + t74)) * MDP(10) + (-g(1) * t51 - g(2) * t49) * MDP(23) + (g(1) * t50 - g(2) * t48) * MDP(24) + (-g(1) * t47 - g(2) * t45) * MDP(30) + (g(1) * t46 - g(2) * t44) * MDP(31) + (MDP(2) - MDP(4) + MDP(9)) * t52 + (-t54 * MDP(16) - t86 - MDP(7) * sin(pkin(10)) - MDP(8) * cos(pkin(10)) + MDP(3) - MDP(5)) * t53; (-MDP(10) - MDP(6)) * t52; -t53 * MDP(10); (t85 * t54 + t86) * g(3) + (MDP(17) * t54 - t85 * t55) * t52; (-g(1) * t48 - g(2) * t50 + t66 * t84) * MDP(23) + (g(1) * t49 - g(2) * t51 + t68 * t84) * MDP(24) + t75; t75;];
taug  = t1;
