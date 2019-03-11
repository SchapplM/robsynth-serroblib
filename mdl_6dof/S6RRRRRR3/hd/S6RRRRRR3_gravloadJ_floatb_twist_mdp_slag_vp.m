% Calculate Gravitation load on the joints for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:42:39
% EndTime: 2019-03-10 03:42:39
% DurationCPUTime: 0.20s
% Computational Cost: add. (384->64), mult. (356->94), div. (0->0), fcn. (350->12), ass. (0->45)
t89 = qJ(2) + qJ(3);
t84 = sin(t89);
t116 = g(3) * t84;
t88 = qJ(4) + qJ(5);
t87 = qJ(6) + t88;
t81 = sin(t87);
t92 = sin(qJ(1));
t114 = t92 * t81;
t82 = cos(t87);
t113 = t92 * t82;
t83 = sin(t88);
t112 = t92 * t83;
t85 = cos(t88);
t111 = t92 * t85;
t90 = sin(qJ(4));
t110 = t92 * t90;
t93 = cos(qJ(4));
t109 = t92 * t93;
t95 = cos(qJ(1));
t108 = t95 * t81;
t107 = t95 * t82;
t106 = t95 * t83;
t105 = t95 * t85;
t104 = t95 * t90;
t103 = t95 * t93;
t86 = cos(t89);
t69 = t114 * t86 + t107;
t70 = -t113 * t86 + t108;
t71 = -t108 * t86 + t113;
t72 = t107 * t86 + t114;
t102 = (-g(1) * t71 + g(2) * t69 + t116 * t81) * MDP(37) + (g(1) * t72 - g(2) * t70 + t116 * t82) * MDP(38);
t73 = t112 * t86 + t105;
t74 = -t111 * t86 + t106;
t75 = -t106 * t86 + t111;
t76 = t105 * t86 + t112;
t101 = (-g(1) * t75 + g(2) * t73 + t116 * t83) * MDP(30) + (g(1) * t76 - g(2) * t74 + t116 * t85) * MDP(31) + t102;
t100 = g(1) * t95 + g(2) * t92;
t98 = (t100 * t86 + t116) * MDP(17) + (t93 * MDP(23) - t90 * MDP(24) + t85 * MDP(30) - t83 * MDP(31) + t82 * MDP(37) - t81 * MDP(38) + MDP(16)) * (-g(3) * t86 + t100 * t84);
t94 = cos(qJ(2));
t91 = sin(qJ(2));
t80 = t103 * t86 + t110;
t79 = -t104 * t86 + t109;
t78 = -t109 * t86 + t104;
t77 = t110 * t86 + t103;
t1 = [t100 * MDP(3) + (-g(1) * t78 - g(2) * t80) * MDP(23) + (-g(1) * t77 - g(2) * t79) * MDP(24) + (-g(1) * t74 - g(2) * t76) * MDP(30) + (-g(1) * t73 - g(2) * t75) * MDP(31) + (-g(1) * t70 - g(2) * t72) * MDP(37) + (-g(1) * t69 - g(2) * t71) * MDP(38) + (-MDP(10) * t91 + MDP(16) * t86 - MDP(17) * t84 + MDP(9) * t94 + MDP(2)) * (g(1) * t92 - g(2) * t95); (-g(3) * t94 + t100 * t91) * MDP(9) + (g(3) * t91 + t100 * t94) * MDP(10) + t98; t98; (-g(1) * t79 + g(2) * t77 + t116 * t90) * MDP(23) + (g(1) * t80 - g(2) * t78 + t93 * t116) * MDP(24) + t101; t101; t102;];
taug  = t1;
