% Calculate Gravitation load on the joints for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:45:13
% EndTime: 2019-03-09 11:45:14
% DurationCPUTime: 0.37s
% Computational Cost: add. (401->77), mult. (396->102), div. (0->0), fcn. (365->10), ass. (0->42)
t102 = sin(qJ(5));
t105 = cos(qJ(5));
t133 = pkin(5) * t105 + qJ(6) * t102;
t132 = MDP(19) - MDP(28);
t129 = MDP(25) + MDP(27);
t128 = MDP(26) - MDP(29);
t104 = sin(qJ(1));
t107 = cos(qJ(1));
t87 = g(1) * t107 + g(2) * t104;
t100 = qJ(2) + pkin(10);
t96 = qJ(4) + t100;
t91 = sin(t96);
t131 = t87 * t91;
t92 = cos(t96);
t130 = t92 * pkin(4) + t91 * pkin(9);
t127 = pkin(9) * t92;
t126 = g(3) * t91;
t101 = -qJ(3) - pkin(7);
t106 = cos(qJ(2));
t97 = t106 * pkin(2);
t122 = pkin(3) * cos(t100) + t97;
t120 = t104 * t102;
t119 = t104 * t105;
t118 = t107 * t102;
t117 = t107 * t105;
t116 = t133 * t92 + t130;
t86 = g(1) * t104 - g(2) * t107;
t114 = pkin(1) + t122 + t130;
t112 = t132 * (t87 * t92 + t126) + (-t128 * t102 + t129 * t105 + MDP(18)) * (-g(3) * t92 + t131);
t76 = t92 * t120 + t117;
t78 = t92 * t118 - t119;
t63 = g(1) * t78 + g(2) * t76 + t102 * t126;
t108 = (pkin(4) + t133) * t131;
t103 = sin(qJ(2));
t99 = -pkin(8) + t101;
t95 = t97 + pkin(1);
t85 = t107 * t127;
t83 = t104 * t127;
t81 = -pkin(3) * sin(t100) - t103 * pkin(2);
t79 = t92 * t117 + t120;
t77 = t92 * t119 - t118;
t1 = [(-g(1) * (-t107 * t101 - t104 * t95) - g(2) * (-t104 * t101 + t107 * t95)) * MDP(12) + (-g(1) * (-t77 * pkin(5) - t76 * qJ(6)) - g(2) * (t79 * pkin(5) + t78 * qJ(6)) + (g(1) * t99 - g(2) * t114) * t107 + (g(1) * t114 + g(2) * t99) * t104) * MDP(30) + (MDP(3) - MDP(11)) * t87 + t129 * (g(1) * t77 - g(2) * t79) - t128 * (g(1) * t76 - g(2) * t78) + (-t103 * MDP(10) + t92 * MDP(18) + t106 * MDP(9) - t132 * t91 + MDP(2)) * t86; (g(3) * t103 + t87 * t106) * MDP(10) + (-g(1) * (t107 * t81 + t85) - g(2) * (t104 * t81 + t83) - g(3) * (t116 + t122) + t108) * MDP(30) + t112 + (MDP(12) * pkin(2) + MDP(9)) * (-g(3) * t106 + t87 * t103); (-MDP(12) - MDP(30)) * t86; (-g(1) * t85 - g(2) * t83 - g(3) * t116 + t108) * MDP(30) + t112; (-g(1) * (-t78 * pkin(5) + t79 * qJ(6)) - g(2) * (-t76 * pkin(5) + t77 * qJ(6)) - (-pkin(5) * t102 + qJ(6) * t105) * t126) * MDP(30) + t129 * t63 + t128 * (g(1) * t79 + g(2) * t77 + t105 * t126); -t63 * MDP(30);];
taug  = t1;
