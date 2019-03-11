% Calculate Gravitation load on the joints for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:04:00
% EndTime: 2019-03-09 22:04:01
% DurationCPUTime: 0.29s
% Computational Cost: add. (370->67), mult. (331->90), div. (0->0), fcn. (284->10), ass. (0->39)
t90 = cos(qJ(2));
t85 = qJ(2) + qJ(3);
t82 = qJ(4) + t85;
t78 = sin(t82);
t79 = cos(t82);
t101 = t79 * pkin(4) + t78 * qJ(5);
t81 = cos(t85);
t99 = pkin(3) * t81 + t101;
t112 = t90 * pkin(2) + t99;
t111 = MDP(23) - MDP(26);
t110 = MDP(24) - MDP(27);
t108 = pkin(4) * t78;
t80 = sin(t85);
t97 = -pkin(3) * t80 - t108;
t106 = g(3) * t79;
t86 = sin(qJ(6));
t88 = sin(qJ(1));
t105 = t88 * t86;
t89 = cos(qJ(6));
t104 = t88 * t89;
t91 = cos(qJ(1));
t103 = t91 * t86;
t102 = t91 * t89;
t100 = qJ(5) * t79;
t87 = sin(qJ(2));
t98 = -t87 * pkin(2) + t97;
t73 = g(1) * t91 + g(2) * t88;
t59 = t73 * t78 - t106;
t95 = t111 * t59 + (-MDP(34) * t86 - MDP(35) * t89 + t110) * (g(3) * t78 + t73 * t79);
t94 = pkin(1) + t112;
t93 = (-g(3) * t81 + t73 * t80) * MDP(16) + (g(3) * t80 + t73 * t81) * MDP(17) + t95;
t84 = -pkin(9) - pkin(8) - pkin(7);
t72 = t91 * t100;
t71 = t88 * t100;
t68 = -t78 * t105 + t102;
t67 = t78 * t104 + t103;
t66 = t78 * t103 + t104;
t65 = t78 * t102 - t105;
t1 = [((g(1) * t84 - g(2) * t94) * t91 + (g(1) * t94 + g(2) * t84) * t88) * MDP(28) + (-g(1) * t68 - g(2) * t66) * MDP(34) + (g(1) * t67 - g(2) * t65) * MDP(35) + (MDP(3) - MDP(25)) * t73 + (-t87 * MDP(10) + MDP(16) * t81 - MDP(17) * t80 + t90 * MDP(9) - t110 * t78 + t111 * t79 + MDP(2)) * (g(1) * t88 - g(2) * t91); (-g(3) * t90 + t73 * t87) * MDP(9) + (g(3) * t87 + t73 * t90) * MDP(10) + (-g(1) * (t91 * t98 + t72) - g(2) * (t88 * t98 + t71) - g(3) * t112) * MDP(28) + t93; (-g(1) * (t91 * t97 + t72) - g(2) * (t88 * t97 + t71) - g(3) * t99) * MDP(28) + t93; (-g(1) * (-t91 * t108 + t72) - g(2) * (-t88 * t108 + t71) - g(3) * t101) * MDP(28) + t95; -t59 * MDP(28); (-g(1) * t65 - g(2) * t67 + t89 * t106) * MDP(34) + (g(1) * t66 - g(2) * t68 - t86 * t106) * MDP(35);];
taug  = t1;
