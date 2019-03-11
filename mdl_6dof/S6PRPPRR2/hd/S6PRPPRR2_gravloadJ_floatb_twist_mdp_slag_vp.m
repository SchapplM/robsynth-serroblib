% Calculate Gravitation load on the joints for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPPRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S6PRPPRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:20:02
% EndTime: 2019-03-08 19:20:03
% DurationCPUTime: 0.45s
% Computational Cost: add. (231->76), mult. (597->137), div. (0->0), fcn. (736->12), ass. (0->40)
t108 = cos(pkin(11));
t87 = sin(pkin(11));
t94 = sin(qJ(2));
t97 = cos(qJ(2));
t81 = -t94 * t108 - t97 * t87;
t102 = t97 * t108 - t94 * t87;
t91 = cos(pkin(6));
t109 = t81 * t91;
t88 = sin(pkin(10));
t90 = cos(pkin(10));
t72 = t102 * t90 + t88 * t109;
t67 = -t102 * t88 + t90 * t109;
t120 = t88 * t94;
t89 = sin(pkin(6));
t93 = sin(qJ(5));
t119 = t89 * t93;
t96 = cos(qJ(5));
t118 = t89 * t96;
t117 = t89 * t97;
t115 = t91 * t94;
t114 = t91 * t97;
t92 = sin(qJ(6));
t113 = t92 * t93;
t95 = cos(qJ(6));
t112 = t93 * t95;
t107 = MDP(5) + MDP(8);
t106 = t90 * t114;
t103 = -t88 * t114 - t90 * t94;
t99 = t91 * t102;
t98 = -g(1) * t103 - g(3) * t117;
t82 = pkin(2) * t106;
t79 = t81 * t89;
t78 = t102 * t89;
t74 = -t78 * t93 + t91 * t96;
t71 = t90 * t81 - t88 * t99;
t68 = t88 * t81 + t90 * t99;
t65 = t90 * t118 + t68 * t93;
t63 = t88 * t118 - t71 * t93;
t61 = g(1) * t71 + g(2) * t68 + g(3) * t78;
t1 = [(-MDP(1) - t107) * g(3); (-g(2) * (t106 - t120) + t98) * MDP(3) + (-g(1) * (t88 * t115 - t90 * t97) - g(2) * (-t90 * t115 - t88 * t97) + g(3) * t89 * t94) * MDP(4) + (-g(2) * t82 + (g(2) * t120 + t98) * pkin(2)) * MDP(5) + t61 * MDP(6) + (-g(1) * (t103 * pkin(2) + t71 * pkin(3) + qJ(4) * t72) - g(2) * (-pkin(2) * t120 + t68 * pkin(3) - t67 * qJ(4) + t82) - g(3) * (pkin(2) * t117 + t78 * pkin(3) - t79 * qJ(4))) * MDP(8) + (-g(1) * (t112 * t72 + t71 * t92) - g(2) * (-t67 * t112 + t68 * t92) - g(3) * (-t79 * t112 + t78 * t92)) * MDP(21) + (-g(1) * (-t113 * t72 + t71 * t95) - g(2) * (t67 * t113 + t68 * t95) - g(3) * (t79 * t113 + t78 * t95)) * MDP(22) + (t93 * MDP(14) + MDP(15) * t96 + MDP(7)) * (-g(1) * t72 + g(2) * t67 + g(3) * t79); t107 * (-g(3) * t91 + (-g(1) * t88 + g(2) * t90) * t89); t61 * MDP(8); (g(1) * t63 - g(2) * t65 + g(3) * t74) * MDP(15) + (-MDP(21) * t95 + MDP(22) * t92 - MDP(14)) * (g(1) * (-t88 * t119 - t71 * t96) + g(2) * (t90 * t119 - t68 * t96) + g(3) * (-t78 * t96 - t91 * t93)); (-g(1) * (-t63 * t92 + t72 * t95) - g(2) * (t65 * t92 - t67 * t95) - g(3) * (-t74 * t92 - t79 * t95)) * MDP(21) + (-g(1) * (-t63 * t95 - t72 * t92) - g(2) * (t65 * t95 + t67 * t92) - g(3) * (-t74 * t95 + t79 * t92)) * MDP(22);];
taug  = t1;
