% Calculate Gravitation load on the joints for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:53:43
% EndTime: 2019-03-09 12:53:44
% DurationCPUTime: 0.54s
% Computational Cost: add. (321->95), mult. (482->127), div. (0->0), fcn. (468->8), ass. (0->50)
t152 = MDP(27) + MDP(29);
t151 = MDP(28) - MDP(31);
t116 = sin(qJ(2));
t119 = cos(qJ(2));
t155 = t119 * pkin(2) + t116 * qJ(3);
t117 = sin(qJ(1));
t120 = cos(qJ(1));
t100 = g(1) * t120 + g(2) * t117;
t115 = sin(qJ(4));
t118 = cos(qJ(4));
t136 = t118 * t120;
t93 = -t115 * t117 + t116 * t136;
t137 = t117 * t118;
t95 = t115 * t120 + t116 * t137;
t154 = -g(1) * t93 - g(2) * t95;
t153 = MDP(10) - MDP(13);
t150 = MDP(9) - MDP(12) + MDP(30);
t147 = pkin(4) * t115;
t146 = pkin(4) * t118;
t145 = g(1) * t117;
t141 = g(3) * t119;
t139 = t116 * t117;
t138 = t116 * t120;
t135 = t119 * t120;
t121 = -pkin(9) - pkin(8);
t134 = t119 * t121;
t133 = t115 * t138;
t132 = g(3) * t155;
t114 = qJ(4) + qJ(5);
t106 = sin(t114);
t107 = cos(t114);
t87 = t106 * t117 - t107 * t138;
t89 = t106 * t120 + t107 * t139;
t81 = g(1) * t87 - g(2) * t89 + t107 * t141;
t88 = t106 * t138 + t107 * t117;
t90 = -t106 * t139 + t107 * t120;
t131 = t152 * t81 - t151 * (-g(1) * t88 + g(2) * t90 + t106 * t141);
t130 = t120 * pkin(1) + pkin(2) * t135 + t117 * pkin(7) + qJ(3) * t138;
t127 = -pkin(5) * t107 - qJ(6) * t106;
t126 = -pkin(1) - t155;
t125 = pkin(5) * t106 - qJ(6) * t107 + t147;
t122 = -g(1) * (-t87 * pkin(5) + qJ(6) * t88) - g(2) * (t89 * pkin(5) - qJ(6) * t90);
t111 = t120 * pkin(7);
t105 = pkin(3) + t146;
t103 = qJ(3) * t135;
t101 = t117 * t119 * qJ(3);
t96 = -t115 * t139 + t136;
t94 = t133 + t137;
t91 = t100 * t116 - t141;
t1 = [(-g(1) * t111 - g(2) * t130 - t126 * t145) * MDP(14) + (-g(1) * t96 - g(2) * t94) * MDP(20) + (g(1) * t95 - g(2) * t93) * MDP(21) + (-g(1) * (t90 * pkin(5) + t89 * qJ(6) + t120 * t105 + t111) - g(2) * (pkin(4) * t133 + t88 * pkin(5) + t87 * qJ(6) - t120 * t134 + t130) + (-g(1) * (-t116 * t147 + t126 + t134) - g(2) * t105) * t117) * MDP(32) + t152 * (-g(1) * t90 - g(2) * t88) + t151 * (g(1) * t89 + g(2) * t87) + (MDP(3) - MDP(11)) * t100 + (-t153 * t116 + t150 * t119 + MDP(2)) * (-g(2) * t120 + t145); (-g(1) * (-pkin(2) * t138 + t103) - g(2) * (-pkin(2) * t139 + t101) - t132) * MDP(14) + (-g(1) * t103 - g(2) * t101 - t132 + (g(3) * t121 - t100 * t125) * t119 + (-g(3) * t125 + t100 * (pkin(2) - t121)) * t116) * MDP(32) + t150 * t91 + (-t115 * MDP(20) - t118 * MDP(21) - t152 * t106 - t151 * t107 + t153) * (g(3) * t116 + t100 * t119); (-MDP(14) - MDP(32)) * t91; (t118 * t141 + t154) * MDP(20) + (g(1) * t94 - g(2) * t96 - t115 * t141) * MDP(21) + (t154 * pkin(4) - (t127 - t146) * t141 + t122) * MDP(32) + t131; (-t127 * t141 + t122) * MDP(32) + t131; -t81 * MDP(32);];
taug  = t1;
