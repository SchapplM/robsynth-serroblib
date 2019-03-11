% Calculate Gravitation load on the joints for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% MDP [34x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(34,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [34 1]), ...
  'S6RPRRRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [34x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:21:13
% EndTime: 2019-03-09 07:21:14
% DurationCPUTime: 0.20s
% Computational Cost: add. (210->54), mult. (252->78), div. (0->0), fcn. (236->10), ass. (0->33)
t67 = qJ(3) + qJ(4);
t64 = cos(t67);
t87 = g(3) * t64;
t66 = qJ(5) + qJ(6);
t61 = sin(t66);
t70 = sin(qJ(1));
t86 = t70 * t61;
t63 = cos(t66);
t85 = t70 * t63;
t68 = sin(qJ(5));
t84 = t70 * t68;
t71 = cos(qJ(5));
t83 = t70 * t71;
t73 = cos(qJ(1));
t82 = t73 * t61;
t81 = t73 * t63;
t80 = t73 * t68;
t79 = t73 * t71;
t62 = sin(t67);
t51 = -t62 * t86 + t81;
t52 = t62 * t85 + t82;
t53 = t62 * t82 + t85;
t54 = t62 * t81 - t86;
t78 = (-g(1) * t51 - g(2) * t53 + t61 * t87) * MDP(33) + (g(1) * t52 - g(2) * t54 + t63 * t87) * MDP(34);
t59 = g(1) * t70 - g(2) * t73;
t76 = (t59 * t62 + t87) * MDP(20) + (-t71 * MDP(26) + t68 * MDP(27) - t63 * MDP(33) + t61 * MDP(34) - MDP(19)) * (-g(3) * t62 + t59 * t64);
t72 = cos(qJ(3));
t69 = sin(qJ(3));
t58 = t62 * t79 - t84;
t57 = t62 * t80 + t83;
t56 = t62 * t83 + t80;
t55 = -t62 * t84 + t79;
t1 = [(-g(1) * (-t70 * pkin(1) + t73 * qJ(2)) - g(2) * (t73 * pkin(1) + t70 * qJ(2))) * MDP(6) + (-g(1) * t58 - g(2) * t56) * MDP(26) + (g(1) * t57 - g(2) * t55) * MDP(27) + (-g(1) * t54 - g(2) * t52) * MDP(33) + (g(1) * t53 - g(2) * t51) * MDP(34) + (MDP(2) - MDP(4)) * t59 + (-t69 * MDP(12) - t72 * MDP(13) - MDP(19) * t62 - MDP(20) * t64 + MDP(3) - MDP(5)) * (g(1) * t73 + g(2) * t70); -t59 * MDP(6); (g(3) * t69 - t59 * t72) * MDP(12) + (g(3) * t72 + t59 * t69) * MDP(13) + t76; t76; (-g(1) * t55 - g(2) * t57 + t68 * t87) * MDP(26) + (g(1) * t56 - g(2) * t58 + t71 * t87) * MDP(27) + t78; t78;];
taug  = t1;
