% Calculate Gravitation load on the joints for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RRPPRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:31:58
% EndTime: 2019-03-09 08:31:59
% DurationCPUTime: 0.47s
% Computational Cost: add. (216->75), mult. (285->98), div. (0->0), fcn. (242->8), ass. (0->43)
t82 = qJ(2) + pkin(9);
t78 = sin(t82);
t79 = cos(t82);
t118 = t79 * pkin(3) + t78 * qJ(4);
t117 = -MDP(14) + MDP(24);
t87 = sin(qJ(1));
t90 = cos(qJ(1));
t115 = -g(1) * t90 - g(2) * t87;
t86 = sin(qJ(2));
t113 = pkin(2) * t86;
t108 = g(3) * t79;
t85 = sin(qJ(5));
t107 = t79 * t85;
t84 = -qJ(3) - pkin(7);
t106 = t84 * t90;
t105 = t85 * t90;
t104 = t87 * t85;
t88 = cos(qJ(5));
t103 = t87 * t88;
t102 = t88 * t90;
t101 = pkin(5) * t88 + pkin(4) - t84;
t100 = qJ(4) * t79;
t99 = -MDP(16) - MDP(25);
t89 = cos(qJ(2));
t80 = t89 * pkin(2);
t97 = t80 + t118;
t77 = t80 + pkin(1);
t72 = t90 * t77;
t96 = g(2) * (t118 * t90 + t72);
t95 = -pkin(3) * t78 - t113;
t69 = g(1) * t87 - g(2) * t90;
t83 = -qJ(6) - pkin(8);
t94 = t78 * t85 * pkin(5) - t79 * t83;
t62 = t102 * t78 - t104;
t64 = t103 * t78 + t105;
t93 = -t77 - t118;
t60 = -g(3) * t78 + t115 * t79;
t68 = t90 * t100;
t66 = t87 * t100;
t65 = -t104 * t78 + t102;
t63 = t105 * t78 + t103;
t59 = -t115 * t78 - t108;
t1 = [(-g(1) * (-t87 * t77 - t106) - g(2) * (-t87 * t84 + t72)) * MDP(12) + (g(1) * t106 - t96 + (-g(1) * t93 + g(2) * t84) * t87) * MDP(16) + (-g(1) * t65 - g(2) * t63) * MDP(22) + (g(1) * t64 - g(2) * t62) * MDP(23) + (-t96 + (-g(1) * t101 - g(2) * t94) * t90 + (-g(1) * (t93 - t94) - g(2) * t101) * t87) * MDP(25) - (MDP(3) - MDP(11) - MDP(13)) * t115 + (-t86 * MDP(10) + t78 * MDP(15) + t89 * MDP(9) + t117 * t79 + MDP(2)) * t69; (g(3) * t86 - t115 * t89) * MDP(10) + (-g(1) * (t90 * t95 + t68) - g(2) * (t87 * t95 + t66) - g(3) * t97) * MDP(16) + (-g(1) * t68 - g(2) * t66 - g(3) * (t94 + t97) + t115 * (pkin(5) * t107 - t113 + (-pkin(3) + t83) * t78)) * MDP(25) + (MDP(12) * pkin(2) + MDP(9)) * (-g(3) * t89 - t115 * t86) + t117 * t59 + (MDP(22) * t85 + MDP(23) * t88 + MDP(15)) * t60; (-MDP(12) + t99) * t69; t99 * t59; (g(1) * t63 - g(2) * t65 - g(3) * t107) * MDP(23) + (MDP(25) * pkin(5) + MDP(22)) * (-g(1) * t62 - g(2) * t64 + t108 * t88); t60 * MDP(25);];
taug  = t1;
