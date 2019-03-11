% Calculate Gravitation load on the joints for
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:14:34
% EndTime: 2019-03-09 07:14:34
% DurationCPUTime: 0.27s
% Computational Cost: add. (320->66), mult. (291->96), div. (0->0), fcn. (290->12), ass. (0->45)
t75 = pkin(11) + qJ(3);
t70 = sin(t75);
t106 = t70 * MDP(14);
t76 = qJ(4) + qJ(5);
t74 = qJ(6) + t76;
t68 = sin(t74);
t69 = cos(t74);
t72 = sin(t76);
t73 = cos(t76);
t79 = sin(qJ(4));
t81 = cos(qJ(4));
t105 = t81 * MDP(20) - t79 * MDP(21) + t73 * MDP(27) - t72 * MDP(28) + t69 * MDP(34) - t68 * MDP(35) + MDP(13);
t104 = g(3) * t70;
t80 = sin(qJ(1));
t103 = t80 * t68;
t102 = t80 * t69;
t101 = t80 * t72;
t100 = t80 * t73;
t99 = t80 * t79;
t98 = t80 * t81;
t82 = cos(qJ(1));
t97 = t82 * t68;
t96 = t82 * t69;
t95 = t82 * t72;
t94 = t82 * t73;
t93 = t82 * t79;
t92 = t82 * t81;
t71 = cos(t75);
t54 = t103 * t71 + t96;
t55 = -t102 * t71 + t97;
t56 = -t71 * t97 + t102;
t57 = t71 * t96 + t103;
t91 = (-g(1) * t56 + g(2) * t54 + t104 * t68) * MDP(34) + (g(1) * t57 - g(2) * t55 + t104 * t69) * MDP(35);
t58 = t101 * t71 + t94;
t59 = -t100 * t71 + t95;
t60 = -t71 * t95 + t100;
t61 = t71 * t94 + t101;
t84 = (-g(1) * t60 + g(2) * t58 + t104 * t72) * MDP(27) + (g(1) * t61 - g(2) * t59 + t104 * t73) * MDP(28) + t91;
t83 = g(1) * t82 + g(2) * t80;
t66 = g(1) * t80 - g(2) * t82;
t65 = t71 * t92 + t99;
t64 = -t71 * t93 + t98;
t63 = -t71 * t98 + t93;
t62 = t71 * t99 + t92;
t1 = [(-g(1) * (-t80 * pkin(1) + t82 * qJ(2)) - g(2) * (t82 * pkin(1) + t80 * qJ(2))) * MDP(7) + (-g(1) * t63 - g(2) * t65) * MDP(20) + (-g(1) * t62 - g(2) * t64) * MDP(21) + (-g(1) * t59 - g(2) * t61) * MDP(27) + (-g(1) * t58 - g(2) * t60) * MDP(28) + (-g(1) * t55 - g(2) * t57) * MDP(34) + (-g(1) * t54 - g(2) * t56) * MDP(35) + (MDP(3) - MDP(6)) * t83 + (t71 * MDP(13) - t106 + MDP(4) * cos(pkin(11)) - MDP(5) * sin(pkin(11)) + MDP(2)) * t66; -t66 * MDP(7); (-t105 * t71 + t106) * g(3) + (MDP(14) * t71 + t105 * t70) * t83; (-g(1) * t64 + g(2) * t62 + t104 * t79) * MDP(20) + (g(1) * t65 - g(2) * t63 + t104 * t81) * MDP(21) + t84; t84; t91;];
taug  = t1;
