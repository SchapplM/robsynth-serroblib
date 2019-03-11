% Calculate Gravitation load on the joints for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:53:26
% EndTime: 2019-03-09 13:53:27
% DurationCPUTime: 0.34s
% Computational Cost: add. (277->63), mult. (442->100), div. (0->0), fcn. (473->10), ass. (0->40)
t89 = sin(qJ(2));
t92 = cos(qJ(4));
t115 = t89 * t92;
t88 = sin(qJ(4));
t93 = cos(qJ(2));
t101 = t93 * t88 - t115;
t110 = qJ(4) + qJ(5);
t109 = cos(t110);
t106 = t89 * t109;
t86 = sin(t110);
t122 = -t93 * t86 + t106;
t116 = g(3) * t122;
t90 = sin(qJ(1));
t64 = t122 * t90;
t95 = t109 * t93 + t89 * t86;
t65 = t95 * t90;
t94 = cos(qJ(1));
t113 = t93 * t94;
t66 = -t106 * t94 + t113 * t86;
t67 = t95 * t94;
t87 = sin(qJ(6));
t91 = cos(qJ(6));
t108 = (MDP(34) * t91 - MDP(35) * t87 + MDP(27)) * (g(1) * t66 - g(2) * t64 + g(3) * t95) + (g(1) * t67 + g(2) * t65 + t116) * MDP(28);
t70 = t101 * t90;
t78 = t89 * t88 + t93 * t92;
t71 = t78 * t90;
t72 = t113 * t88 - t115 * t94;
t73 = t78 * t94;
t125 = t108 + (g(1) * t72 + g(2) * t70 + g(3) * t78) * MDP(20) + (g(1) * t73 + g(2) * t71 - g(3) * t101) * MDP(21);
t82 = g(1) * t94 + g(2) * t90;
t121 = MDP(9) + MDP(11);
t120 = MDP(10) - MDP(13);
t105 = t93 * pkin(2) + t89 * qJ(3);
t103 = t65 * t91 + t94 * t87;
t102 = t65 * t87 - t94 * t91;
t100 = pkin(1) + t105;
t68 = -g(3) * t93 + t82 * t89;
t63 = t67 * t91 - t90 * t87;
t62 = -t67 * t87 - t90 * t91;
t1 = [((-g(1) * pkin(7) - g(2) * t100) * t94 + (-g(2) * pkin(7) + g(1) * t100) * t90) * MDP(14) + (g(1) * t71 - g(2) * t73) * MDP(20) + (-g(1) * t70 + g(2) * t72) * MDP(21) + (g(1) * t65 - g(2) * t67) * MDP(27) + (g(1) * t64 + g(2) * t66) * MDP(28) + (g(1) * t103 - g(2) * t63) * MDP(34) + (-g(1) * t102 - g(2) * t62) * MDP(35) + (MDP(3) - MDP(12)) * t82 + (-t120 * t89 + t121 * t93 + MDP(2)) * (g(1) * t90 - g(2) * t94); (-g(3) * t105 + t82 * (pkin(2) * t89 - qJ(3) * t93)) * MDP(14) + t120 * (g(3) * t89 + t82 * t93) + t121 * t68 - t125; -t68 * MDP(14); t125; t108; (-g(1) * t62 + g(2) * t102 + t116 * t87) * MDP(34) + (g(1) * t63 + g(2) * t103 + t116 * t91) * MDP(35);];
taug  = t1;
