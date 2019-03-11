% Calculate Gravitation load on the joints for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:16:45
% EndTime: 2019-03-09 06:16:46
% DurationCPUTime: 0.32s
% Computational Cost: add. (268->72), mult. (293->103), div. (0->0), fcn. (267->10), ass. (0->43)
t91 = sin(qJ(1));
t93 = cos(qJ(1));
t75 = g(1) * t93 + g(2) * t91;
t117 = MDP(14) - MDP(29);
t85 = pkin(10) + qJ(3);
t78 = sin(t85);
t111 = g(3) * t78;
t86 = qJ(4) + qJ(5);
t80 = sin(t86);
t105 = t91 * t80;
t81 = cos(t86);
t107 = t81 * t93;
t79 = cos(t85);
t62 = t79 * t105 + t107;
t104 = t91 * t81;
t108 = t80 * t93;
t64 = -t79 * t108 + t104;
t116 = -g(1) * t64 + g(2) * t62 + t80 * t111;
t60 = -g(3) * t79 + t75 * t78;
t90 = sin(qJ(4));
t72 = pkin(4) * t90 + pkin(5) * t80;
t109 = t72 * t79;
t106 = t90 * t93;
t103 = t91 * t90;
t92 = cos(qJ(4));
t102 = t91 * t92;
t101 = t92 * t93;
t63 = -t79 * t104 + t108;
t65 = t79 * t107 + t105;
t100 = t116 * MDP(27) + (g(1) * t65 - g(2) * t63 + t81 * t111) * MDP(28);
t99 = t72 + pkin(7) + qJ(2);
t73 = t92 * pkin(4) + pkin(5) * t81;
t74 = g(1) * t91 - g(2) * t93;
t71 = pkin(3) + t73;
t84 = -qJ(6) - pkin(9) - pkin(8);
t97 = t71 * t79 - t78 * t84;
t88 = cos(pkin(10));
t95 = pkin(2) * t88 + pkin(1) + t97;
t70 = t79 * t101 + t103;
t69 = -t79 * t106 + t102;
t68 = -t79 * t102 + t106;
t67 = t79 * t103 + t101;
t1 = [(-g(1) * (-t91 * pkin(1) + qJ(2) * t93) - g(2) * (pkin(1) * t93 + t91 * qJ(2))) * MDP(7) + (-g(1) * t68 - g(2) * t70) * MDP(20) + (-g(1) * t67 - g(2) * t69) * MDP(21) + (-g(1) * t63 - g(2) * t65) * MDP(27) + (-g(1) * t62 - g(2) * t64) * MDP(28) + ((-g(1) * t99 - g(2) * t95) * t93 + (g(1) * t95 - g(2) * t99) * t91) * MDP(30) + (MDP(3) - MDP(6)) * t75 + (t79 * MDP(13) + MDP(4) * t88 - MDP(5) * sin(pkin(10)) - t117 * t78 + MDP(2)) * t74; (-MDP(30) - MDP(7)) * t74; (-g(3) * t97 + t75 * (t71 * t78 + t79 * t84)) * MDP(30) + t117 * (t75 * t79 + t111) + (t92 * MDP(20) - t90 * MDP(21) + t81 * MDP(27) - t80 * MDP(28) + MDP(13)) * t60; (-g(1) * t69 + g(2) * t67 + t90 * t111) * MDP(20) + (g(1) * t70 - g(2) * t68 + t92 * t111) * MDP(21) + (-g(1) * (-t93 * t109 + t91 * t73) - g(2) * (-t91 * t109 - t73 * t93) + t72 * t111) * MDP(30) + t100; t116 * MDP(30) * pkin(5) + t100; -t60 * MDP(30);];
taug  = t1;
