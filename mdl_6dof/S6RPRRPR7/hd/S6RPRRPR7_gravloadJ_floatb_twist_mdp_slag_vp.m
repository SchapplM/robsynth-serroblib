% Calculate Gravitation load on the joints for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:21:01
% EndTime: 2019-03-09 05:21:02
% DurationCPUTime: 0.20s
% Computational Cost: add. (161->54), mult. (199->75), div. (0->0), fcn. (167->10), ass. (0->30)
t60 = qJ(3) + qJ(4);
t52 = pkin(10) + t60;
t51 = cos(t52);
t74 = g(3) * t51;
t61 = sin(qJ(6));
t66 = cos(qJ(1));
t73 = t61 * t66;
t63 = sin(qJ(1));
t72 = t63 * t61;
t64 = cos(qJ(6));
t71 = t63 * t64;
t70 = t64 * t66;
t69 = t66 * pkin(1) + t63 * qJ(2);
t48 = g(1) * t63 - g(2) * t66;
t50 = sin(t52);
t53 = sin(t60);
t54 = cos(t60);
t67 = g(3) * t53 - t48 * t54;
t68 = (g(3) * t54 + t48 * t53) * MDP(20) + t67 * MDP(19) + (-t64 * MDP(28) + t61 * MDP(29)) * (-g(3) * t50 + t48 * t51);
t49 = g(1) * t66 + g(2) * t63;
t65 = cos(qJ(3));
t62 = sin(qJ(3));
t59 = -qJ(5) - pkin(8) - pkin(7);
t56 = t66 * qJ(2);
t46 = pkin(3) * t62 + pkin(4) * t53;
t45 = t50 * t70 - t72;
t44 = t50 * t73 + t71;
t43 = t50 * t71 + t73;
t42 = -t50 * t72 + t70;
t1 = [(-g(1) * (-t63 * pkin(1) + t56) - g(2) * t69) * MDP(6) + (-g(1) * (t46 * t66 + t56 + (-pkin(1) + t59) * t63) - g(2) * (t63 * t46 - t59 * t66 + t69)) * MDP(22) + (-g(1) * t45 - g(2) * t43) * MDP(28) + (g(1) * t44 - g(2) * t42) * MDP(29) + (MDP(2) - MDP(4) + MDP(21)) * t48 + (-t62 * MDP(12) - t65 * MDP(13) - MDP(19) * t53 - MDP(20) * t54 + MDP(3) - MDP(5)) * t49; (-MDP(22) - MDP(6)) * t48; (g(3) * t62 - t48 * t65) * MDP(12) + (g(3) * t65 + t48 * t62) * MDP(13) + (g(3) * t46 - t48 * (pkin(3) * t65 + pkin(4) * t54)) * MDP(22) + t68; t67 * MDP(22) * pkin(4) + t68; -t49 * MDP(22); (-g(1) * t42 - g(2) * t44 + t61 * t74) * MDP(28) + (g(1) * t43 - g(2) * t45 + t64 * t74) * MDP(29);];
taug  = t1;
