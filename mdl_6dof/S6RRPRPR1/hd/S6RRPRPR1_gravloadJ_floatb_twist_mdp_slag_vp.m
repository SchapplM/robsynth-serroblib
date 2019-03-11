% Calculate Gravitation load on the joints for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:09:57
% EndTime: 2019-03-09 10:09:58
% DurationCPUTime: 0.32s
% Computational Cost: add. (333->74), mult. (297->101), div. (0->0), fcn. (258->12), ass. (0->44)
t91 = qJ(2) + pkin(10);
t86 = qJ(4) + t91;
t81 = sin(t86);
t82 = cos(t86);
t106 = t82 * pkin(4) + t81 * qJ(5);
t97 = cos(qJ(2));
t87 = t97 * pkin(2);
t119 = pkin(3) * cos(t91) + t87 + t106;
t118 = MDP(19) - MDP(22);
t96 = sin(qJ(1));
t98 = cos(qJ(1));
t77 = g(1) * t98 + g(2) * t96;
t65 = -g(3) * t82 + t77 * t81;
t117 = pkin(4) * t81;
t116 = g(3) * t81;
t90 = pkin(11) + qJ(6);
t84 = sin(t90);
t114 = t96 * t84;
t85 = cos(t90);
t113 = t96 * t85;
t92 = sin(pkin(11));
t112 = t96 * t92;
t93 = cos(pkin(11));
t111 = t96 * t93;
t110 = t98 * t84;
t109 = t98 * t85;
t108 = t98 * t92;
t107 = t98 * t93;
t94 = -qJ(3) - pkin(7);
t104 = qJ(5) * t82;
t95 = sin(qJ(2));
t103 = -pkin(3) * sin(t91) - t95 * pkin(2) - t117;
t76 = g(1) * t96 - g(2) * t98;
t102 = pkin(1) + t119;
t101 = t118 * (t77 * t82 + t116) + (MDP(20) * t93 - MDP(21) * t92 + MDP(29) * t85 - MDP(30) * t84 + MDP(18)) * t65;
t89 = -pkin(8) + t94;
t83 = t87 + pkin(1);
t75 = t98 * t104;
t74 = t96 * t104;
t70 = t82 * t109 + t114;
t69 = -t82 * t110 + t113;
t68 = -t82 * t113 + t110;
t67 = t82 * t114 + t109;
t1 = [(-g(1) * (-t96 * t83 - t98 * t94) - g(2) * (t98 * t83 - t96 * t94)) * MDP(12) + (-g(1) * (-t82 * t111 + t108) - g(2) * (t82 * t107 + t112)) * MDP(20) + (-g(1) * (t82 * t112 + t107) - g(2) * (-t82 * t108 + t111)) * MDP(21) + ((g(1) * t89 - g(2) * t102) * t98 + (g(1) * t102 + g(2) * t89) * t96) * MDP(23) + (-g(1) * t68 - g(2) * t70) * MDP(29) + (-g(1) * t67 - g(2) * t69) * MDP(30) + (MDP(3) - MDP(11)) * t77 + (-t95 * MDP(10) + t82 * MDP(18) + t97 * MDP(9) - t118 * t81 + MDP(2)) * t76; (g(3) * t95 + t77 * t97) * MDP(10) + (-g(1) * (t103 * t98 + t75) - g(2) * (t103 * t96 + t74) - g(3) * t119) * MDP(23) + t101 + (MDP(12) * pkin(2) + MDP(9)) * (-g(3) * t97 + t77 * t95); (-MDP(12) - MDP(23)) * t76; (-g(1) * (-t98 * t117 + t75) - g(2) * (-t96 * t117 + t74) - g(3) * t106) * MDP(23) + t101; -t65 * MDP(23); (-g(1) * t69 + g(2) * t67 + t84 * t116) * MDP(29) + (g(1) * t70 - g(2) * t68 + t85 * t116) * MDP(30);];
taug  = t1;
