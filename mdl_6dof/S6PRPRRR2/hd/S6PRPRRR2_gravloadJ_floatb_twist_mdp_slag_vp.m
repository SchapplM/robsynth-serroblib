% Calculate Gravitation load on the joints for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:29:56
% EndTime: 2019-03-08 20:29:58
% DurationCPUTime: 0.49s
% Computational Cost: add. (372->86), mult. (864->158), div. (0->0), fcn. (1103->14), ass. (0->39)
t100 = cos(pkin(12));
t81 = sin(pkin(12));
t88 = sin(qJ(2));
t91 = cos(qJ(2));
t96 = t91 * t100 - t88 * t81;
t83 = sin(pkin(6));
t111 = g(3) * t83;
t80 = qJ(5) + qJ(6);
t78 = sin(t80);
t90 = cos(qJ(4));
t110 = t78 * t90;
t79 = cos(t80);
t109 = t79 * t90;
t87 = sin(qJ(4));
t108 = t83 * t87;
t107 = t83 * t90;
t85 = cos(pkin(6));
t106 = t85 * t88;
t105 = t85 * t91;
t86 = sin(qJ(5));
t104 = t86 * t90;
t89 = cos(qJ(5));
t102 = t89 * t90;
t84 = cos(pkin(11));
t75 = -t88 * t100 - t91 * t81;
t73 = t75 * t85;
t82 = sin(pkin(11));
t98 = -t84 * t73 + t82 * t96;
t59 = -t84 * t108 + t90 * t98;
t97 = t82 * t73 + t84 * t96;
t61 = t82 * t108 + t90 * t97;
t93 = t96 * t85;
t63 = t82 * t75 + t84 * t93;
t66 = t84 * t75 - t82 * t93;
t72 = t75 * t83;
t69 = -t72 * t90 + t85 * t87;
t71 = t96 * t83;
t101 = (-g(1) * (-t61 * t78 - t66 * t79) - g(2) * (-t59 * t78 - t63 * t79) - g(3) * (-t69 * t78 - t71 * t79)) * MDP(25) + (-g(1) * (-t61 * t79 + t66 * t78) - g(2) * (-t59 * t79 + t63 * t78) - g(3) * (-t69 * t79 + t71 * t78)) * MDP(26);
t1 = [(-MDP(1) - MDP(5)) * g(3); (-g(1) * (t82 * t106 - t84 * t91) - g(2) * (-t84 * t106 - t82 * t91) + t88 * t111) * MDP(4) + (-g(1) * (t66 * t102 + t86 * t97) - g(2) * (t63 * t102 + t86 * t98) - g(3) * (t71 * t102 - t72 * t86)) * MDP(18) + (-g(1) * (-t66 * t104 + t89 * t97) - g(2) * (-t63 * t104 + t89 * t98) - g(3) * (-t71 * t104 - t72 * t89)) * MDP(19) + (-g(1) * (t66 * t109 + t78 * t97) - g(2) * (t63 * t109 + t78 * t98) - g(3) * (t71 * t109 - t72 * t78)) * MDP(25) + (-g(1) * (-t66 * t110 + t79 * t97) - g(2) * (-t63 * t110 + t79 * t98) - g(3) * (-t71 * t110 - t72 * t79)) * MDP(26) + (-t90 * MDP(11) + MDP(12) * t87) * (g(1) * t66 + g(2) * t63 + g(3) * t71) + (pkin(2) * MDP(5) + MDP(3)) * (-g(1) * (-t82 * t105 - t84 * t88) - g(2) * (t84 * t105 - t82 * t88) - t91 * t111); (-g(3) * t85 + (-g(1) * t82 + g(2) * t84) * t83) * MDP(5); (g(1) * t61 + g(2) * t59 + g(3) * t69) * MDP(12) + (-MDP(18) * t89 + MDP(19) * t86 - MDP(25) * t79 + MDP(26) * t78 - MDP(11)) * (g(1) * (t82 * t107 - t87 * t97) + g(2) * (-t84 * t107 - t87 * t98) + g(3) * (t72 * t87 + t85 * t90)); (-g(1) * (-t61 * t86 - t66 * t89) - g(2) * (-t59 * t86 - t63 * t89) - g(3) * (-t69 * t86 - t71 * t89)) * MDP(18) + (-g(1) * (-t61 * t89 + t66 * t86) - g(2) * (-t59 * t89 + t63 * t86) - g(3) * (-t69 * t89 + t71 * t86)) * MDP(19) + t101; t101;];
taug  = t1;
