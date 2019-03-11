% Calculate Gravitation load on the joints for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPPRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:23:03
% EndTime: 2019-03-08 19:23:04
% DurationCPUTime: 0.34s
% Computational Cost: add. (233->77), mult. (591->133), div. (0->0), fcn. (725->12), ass. (0->39)
t105 = cos(qJ(2));
t97 = sin(pkin(6));
t120 = t105 * t97;
t102 = sin(qJ(2));
t122 = t102 * t97;
t124 = pkin(2) * t120 + qJ(3) * t122;
t101 = sin(qJ(5));
t123 = t101 * t97;
t104 = cos(qJ(5));
t121 = t104 * t97;
t119 = cos(pkin(6));
t100 = sin(qJ(6));
t118 = t100 * t104;
t103 = cos(qJ(6));
t117 = t103 * t104;
t116 = MDP(10) + MDP(7);
t110 = t105 * t119;
t96 = sin(pkin(10));
t99 = cos(pkin(10));
t86 = t96 * t102 - t110 * t99;
t111 = t102 * t119;
t87 = t96 * t105 + t111 * t99;
t95 = sin(pkin(11));
t98 = cos(pkin(11));
t115 = -t86 * t98 + t87 * t95;
t71 = t86 * t95 + t87 * t98;
t88 = t99 * t102 + t110 * t96;
t89 = t99 * t105 - t111 * t96;
t114 = -t88 * t98 + t89 * t95;
t75 = t88 * t95 + t89 * t98;
t113 = -t86 * pkin(2) + t87 * qJ(3);
t112 = -t88 * pkin(2) + t89 * qJ(3);
t106 = -g(1) * t88 - g(2) * t86 + g(3) * t120;
t83 = -t120 * t95 + t122 * t98;
t82 = (t102 * t95 + t105 * t98) * t97;
t77 = -t101 * t119 + t83 * t104;
t65 = t75 * t104 - t123 * t96;
t63 = t71 * t104 + t123 * t99;
t1 = [(-MDP(1) - t116) * g(3); (-g(1) * t112 - g(2) * t113 - g(3) * t124) * MDP(7) + (-g(1) * t75 - g(2) * t71 - g(3) * t83) * MDP(9) + (-g(1) * (-t88 * pkin(3) + t112) - g(2) * (-t86 * pkin(3) + t113) - g(3) * (pkin(3) * t120 + t124)) * MDP(10) + (-g(1) * (-t100 * t75 + t114 * t117) - g(2) * (-t100 * t71 + t115 * t117) - g(3) * (-t83 * t100 + t117 * t82)) * MDP(23) + (-g(1) * (-t103 * t75 - t114 * t118) - g(2) * (-t103 * t71 - t115 * t118) - g(3) * (-t83 * t103 - t118 * t82)) * MDP(24) - (MDP(3) + MDP(5)) * t106 + (MDP(4) - MDP(6)) * (g(1) * t89 + g(2) * t87 + g(3) * t122) + (-MDP(16) * t104 + MDP(17) * t101 - MDP(8)) * (g(1) * t114 + g(2) * t115 + g(3) * t82); t116 * t106; (g(3) * t119 + (g(1) * t96 - g(2) * t99) * t97) * MDP(10); (g(1) * t65 + g(2) * t63 + g(3) * t77) * MDP(17) + (-MDP(23) * t103 + MDP(24) * t100 - MDP(16)) * (g(1) * (-t75 * t101 - t121 * t96) + g(2) * (-t71 * t101 + t121 * t99) + g(3) * (-t83 * t101 - t104 * t119)); (-g(1) * (-t65 * t100 + t103 * t114) - g(2) * (-t63 * t100 + t103 * t115) - g(3) * (-t77 * t100 + t82 * t103)) * MDP(23) + (-g(1) * (-t100 * t114 - t65 * t103) - g(2) * (-t100 * t115 - t63 * t103) - g(3) * (-t82 * t100 - t77 * t103)) * MDP(24);];
taug  = t1;
