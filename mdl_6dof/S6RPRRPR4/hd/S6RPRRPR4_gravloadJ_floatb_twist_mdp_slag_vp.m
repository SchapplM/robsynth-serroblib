% Calculate Gravitation load on the joints for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:10:17
% EndTime: 2019-03-09 05:10:18
% DurationCPUTime: 0.42s
% Computational Cost: add. (337->72), mult. (293->101), div. (0->0), fcn. (256->12), ass. (0->42)
t86 = pkin(10) + qJ(3);
t81 = cos(t86);
t82 = qJ(4) + t86;
t76 = sin(t82);
t77 = cos(t82);
t98 = t77 * pkin(4) + t76 * qJ(5);
t111 = pkin(3) * t81 + t98;
t110 = MDP(21) - MDP(24);
t91 = sin(qJ(1));
t92 = cos(qJ(1));
t72 = g(1) * t92 + g(2) * t91;
t61 = -g(3) * t77 + t72 * t76;
t109 = pkin(4) * t76;
t108 = g(3) * t76;
t85 = pkin(11) + qJ(6);
t78 = sin(t85);
t106 = t78 * t92;
t80 = cos(t85);
t105 = t80 * t92;
t87 = sin(pkin(11));
t104 = t87 * t92;
t89 = cos(pkin(11));
t103 = t89 * t92;
t102 = t91 * t78;
t101 = t91 * t80;
t100 = t91 * t87;
t99 = t91 * t89;
t97 = qJ(5) * t77;
t79 = sin(t86);
t96 = -pkin(3) * t79 - t109;
t71 = g(1) * t91 - g(2) * t92;
t90 = cos(pkin(10));
t95 = pkin(2) * t90 + pkin(1) + t111;
t94 = t110 * (t72 * t77 + t108) + (MDP(22) * t89 - MDP(23) * t87 + MDP(31) * t80 - MDP(32) * t78 + MDP(20)) * t61;
t84 = -pkin(8) - pkin(7) - qJ(2);
t70 = t92 * t97;
t69 = t91 * t97;
t66 = t77 * t105 + t102;
t65 = -t77 * t106 + t101;
t64 = -t77 * t101 + t106;
t63 = t77 * t102 + t105;
t1 = [(-g(1) * (-t91 * pkin(1) + qJ(2) * t92) - g(2) * (pkin(1) * t92 + t91 * qJ(2))) * MDP(7) + (-g(1) * (-t77 * t99 + t104) - g(2) * (t77 * t103 + t100)) * MDP(22) + (-g(1) * (t77 * t100 + t103) - g(2) * (-t77 * t104 + t99)) * MDP(23) + ((g(1) * t84 - g(2) * t95) * t92 + (g(1) * t95 + g(2) * t84) * t91) * MDP(25) + (-g(1) * t64 - g(2) * t66) * MDP(31) + (-g(1) * t63 - g(2) * t65) * MDP(32) + (MDP(3) - MDP(6)) * t72 + (-t110 * t76 + t81 * MDP(13) - t79 * MDP(14) + t77 * MDP(20) + MDP(4) * t90 - MDP(5) * sin(pkin(10)) + MDP(2)) * t71; (-MDP(25) - MDP(7)) * t71; (-g(3) * t81 + t72 * t79) * MDP(13) + (g(3) * t79 + t72 * t81) * MDP(14) + (-g(1) * (t96 * t92 + t70) - g(2) * (t96 * t91 + t69) - g(3) * t111) * MDP(25) + t94; (-g(1) * (-t92 * t109 + t70) - g(2) * (-t91 * t109 + t69) - g(3) * t98) * MDP(25) + t94; -t61 * MDP(25); (-g(1) * t65 + g(2) * t63 + t78 * t108) * MDP(31) + (g(1) * t66 - g(2) * t64 + t80 * t108) * MDP(32);];
taug  = t1;
