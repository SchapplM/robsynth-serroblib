% Calculate Gravitation load on the joints for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:13:55
% EndTime: 2019-03-09 10:13:56
% DurationCPUTime: 0.29s
% Computational Cost: add. (277->65), mult. (263->85), div. (0->0), fcn. (222->10), ass. (0->38)
t88 = cos(qJ(2));
t79 = t88 * pkin(2);
t82 = qJ(2) + pkin(10);
t78 = qJ(4) + t82;
t75 = sin(t78);
t76 = cos(t78);
t97 = t76 * pkin(4) + t75 * qJ(5);
t107 = pkin(3) * cos(t82) + t79 + t97;
t106 = MDP(18) - MDP(21);
t105 = MDP(19) - MDP(22);
t104 = pkin(4) * t75;
t102 = g(3) * t76;
t84 = sin(qJ(6));
t86 = sin(qJ(1));
t101 = t86 * t84;
t87 = cos(qJ(6));
t100 = t86 * t87;
t89 = cos(qJ(1));
t99 = t89 * t84;
t98 = t89 * t87;
t83 = -qJ(3) - pkin(7);
t95 = qJ(5) * t76;
t85 = sin(qJ(2));
t94 = -pkin(3) * sin(t82) - t85 * pkin(2) - t104;
t70 = g(1) * t89 + g(2) * t86;
t69 = g(1) * t86 - g(2) * t89;
t57 = t70 * t75 - t102;
t93 = t106 * t57 + (-MDP(29) * t84 - MDP(30) * t87 + t105) * (g(3) * t75 + t70 * t76);
t92 = pkin(1) + t107;
t81 = -pkin(8) + t83;
t77 = t79 + pkin(1);
t68 = t89 * t95;
t67 = t86 * t95;
t64 = -t101 * t75 + t98;
t63 = t100 * t75 + t99;
t62 = t75 * t99 + t100;
t61 = t75 * t98 - t101;
t1 = [(-g(1) * (-t86 * t77 - t89 * t83) - g(2) * (t89 * t77 - t86 * t83)) * MDP(12) + ((g(1) * t81 - g(2) * t92) * t89 + (g(1) * t92 + g(2) * t81) * t86) * MDP(23) + (-g(1) * t64 - g(2) * t62) * MDP(29) + (g(1) * t63 - g(2) * t61) * MDP(30) + (MDP(3) - MDP(11) - MDP(20)) * t70 + (-t85 * MDP(10) + t88 * MDP(9) - t105 * t75 + t106 * t76 + MDP(2)) * t69; (g(3) * t85 + t70 * t88) * MDP(10) + (-g(1) * (t89 * t94 + t68) - g(2) * (t86 * t94 + t67) - g(3) * t107) * MDP(23) + t93 + (MDP(12) * pkin(2) + MDP(9)) * (-g(3) * t88 + t70 * t85); (-MDP(12) - MDP(23)) * t69; (-g(1) * (-t104 * t89 + t68) - g(2) * (-t104 * t86 + t67) - g(3) * t97) * MDP(23) + t93; -t57 * MDP(23); (-g(1) * t61 - g(2) * t63 + t102 * t87) * MDP(29) + (g(1) * t62 - g(2) * t64 - t102 * t84) * MDP(30);];
taug  = t1;
