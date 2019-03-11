% Calculate Gravitation load on the joints for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:23:04
% EndTime: 2019-03-10 01:23:05
% DurationCPUTime: 0.30s
% Computational Cost: add. (381->86), mult. (400->127), div. (0->0), fcn. (384->10), ass. (0->45)
t102 = cos(qJ(1));
t99 = sin(qJ(1));
t107 = g(1) * t102 + g(2) * t99;
t124 = MDP(10) - MDP(32);
t98 = sin(qJ(2));
t119 = g(3) * t98;
t101 = cos(qJ(2));
t114 = t101 * t99;
t96 = qJ(3) + qJ(4);
t94 = qJ(5) + t96;
t89 = sin(t94);
t90 = cos(t94);
t67 = t102 * t90 + t89 * t114;
t112 = t101 * t102;
t69 = -t89 * t112 + t99 * t90;
t123 = -g(1) * t69 + g(2) * t67 + t89 * t119;
t76 = -g(3) * t101 + t107 * t98;
t91 = sin(t96);
t85 = -pkin(4) * t91 - pkin(5) * t89;
t97 = sin(qJ(3));
t78 = t97 * pkin(3) - t85;
t118 = pkin(7) + t78;
t68 = t102 * t89 - t90 * t114;
t70 = t90 * t112 + t99 * t89;
t115 = t123 * MDP(30) + (g(1) * t70 - g(2) * t68 + t90 * t119) * MDP(31);
t92 = cos(t96);
t86 = pkin(4) * t92 + pkin(5) * t90;
t100 = cos(qJ(3));
t113 = t99 * t100;
t111 = t102 * t100;
t79 = t100 * pkin(3) + t86;
t72 = t102 * t92 + t91 * t114;
t73 = t102 * t91 - t92 * t114;
t74 = -t91 * t112 + t99 * t92;
t75 = t92 * t112 + t99 * t91;
t109 = (-g(1) * t74 + g(2) * t72 + t91 * t119) * MDP(23) + (g(1) * t75 - g(2) * t73 + t92 * t119) * MDP(24) + t115;
t71 = pkin(2) + t79;
t93 = -qJ(6) - pkin(10) - pkin(9) - pkin(8);
t105 = t101 * t71 - t98 * t93;
t104 = pkin(1) + t105;
t83 = t101 * t111 + t99 * t97;
t82 = -t97 * t112 + t113;
t81 = -t101 * t113 + t102 * t97;
t80 = t97 * t114 + t111;
t1 = [t107 * MDP(3) + (-g(1) * t81 - g(2) * t83) * MDP(16) + (-g(1) * t80 - g(2) * t82) * MDP(17) + (-g(1) * t73 - g(2) * t75) * MDP(23) + (-g(1) * t72 - g(2) * t74) * MDP(24) + (-g(1) * t68 - g(2) * t70) * MDP(30) + (-g(1) * t67 - g(2) * t69) * MDP(31) + ((g(1) * t104 - g(2) * t118) * t99 + (-g(1) * t118 - g(2) * t104) * t102) * MDP(33) + (t101 * MDP(9) - t124 * t98 + MDP(2)) * (g(1) * t99 - g(2) * t102); (-g(3) * t105 + t107 * (t101 * t93 + t71 * t98)) * MDP(33) + t124 * (t107 * t101 + t119) + (t100 * MDP(16) - t97 * MDP(17) + t92 * MDP(23) - t91 * MDP(24) + t90 * MDP(30) - t89 * MDP(31) + MDP(9)) * t76; (-g(1) * t82 + g(2) * t80 + t97 * t119) * MDP(16) + (g(1) * t83 - g(2) * t81 + t100 * t119) * MDP(17) + (-g(1) * (-t78 * t112 + t99 * t79) - g(2) * (-t102 * t79 - t78 * t114) + t78 * t119) * MDP(33) + t109; (-g(1) * (t85 * t112 + t99 * t86) - g(2) * (-t102 * t86 + t85 * t114) - t85 * t119) * MDP(33) + t109; t123 * MDP(33) * pkin(5) + t115; -t76 * MDP(33);];
taug  = t1;
