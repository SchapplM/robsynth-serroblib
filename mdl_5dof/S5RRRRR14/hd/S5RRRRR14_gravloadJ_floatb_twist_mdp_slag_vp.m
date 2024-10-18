% Calculate Gravitation load on the joints for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR14_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR14_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR14_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:37
% EndTime: 2024-09-27 18:43:38
% DurationCPUTime: 0.17s
% Computational Cost: add. (557->63), mult. (304->93), div. (0->0), fcn. (276->20), ass. (0->49)
t127 = g(3) * sin(pkin(5));
t104 = qJ(3) + qJ(4);
t121 = pkin(5) - t104;
t119 = -qJ(5) + t121;
t113 = sin(t119);
t105 = qJ(1) + qJ(2);
t100 = sin(t105);
t102 = cos(t105);
t98 = pkin(5) + t104;
t95 = qJ(5) + t98;
t123 = sin(t95) / 0.2e1;
t81 = t123 - t113 / 0.2e1;
t103 = qJ(5) + t104;
t97 = cos(t103);
t116 = t100 * t81 - t102 * t97;
t117 = t100 * t97 + t102 * t81;
t89 = cos(t95) / 0.2e1;
t93 = cos(t119);
t82 = t93 / 0.2e1 + t89;
t96 = sin(t103);
t73 = t100 * t96 - t102 * t82;
t74 = -t100 * t82 - t102 * t96;
t126 = (-g(1) * t74 + g(2) * t73 - g(3) * (t123 + t113 / 0.2e1)) * MDP(26) + (-g(1) * t116 + g(2) * t117 - g(3) * (t89 - t93 / 0.2e1)) * MDP(27);
t107 = cos(pkin(5));
t108 = sin(qJ(3));
t125 = t107 * t108;
t110 = cos(qJ(3));
t124 = t107 * t110;
t122 = sin(t98) / 0.2e1;
t101 = cos(t104);
t118 = sin(t121);
t85 = t122 - t118 / 0.2e1;
t114 = t100 * t101 + t102 * t85;
t115 = t100 * t85 - t102 * t101;
t92 = cos(t98) / 0.2e1;
t94 = cos(t121);
t86 = t94 / 0.2e1 + t92;
t99 = sin(t104);
t75 = t100 * t99 - t102 * t86;
t76 = -t100 * t86 - t102 * t99;
t120 = (-g(1) * t76 + g(2) * t75 - g(3) * (t122 + t118 / 0.2e1)) * MDP(19) + (-g(1) * t115 + g(2) * t114 - g(3) * (t92 - t94 / 0.2e1)) * MDP(20) + t126;
t77 = t100 * t108 - t102 * t124;
t78 = -t100 * t110 - t102 * t125;
t79 = -t100 * t124 - t102 * t108;
t80 = -t100 * t125 + t102 * t110;
t112 = (g(1) * t100 - g(2) * t102) * MDP(5) + (g(1) * t102 + g(2) * t100) * MDP(6) + (-g(1) * t78 - g(2) * t80) * MDP(12) + (-g(1) * t77 - g(2) * t79) * MDP(13) + (g(1) * t114 + g(2) * t115) * MDP(19) + (-g(1) * t75 - g(2) * t76) * MDP(20) + (g(1) * t117 + g(2) * t116) * MDP(26) + (-g(1) * t73 - g(2) * t74) * MDP(27);
t111 = cos(qJ(1));
t109 = sin(qJ(1));
t1 = [(g(1) * t109 - g(2) * t111) * MDP(2) + (g(1) * t111 + g(2) * t109) * MDP(3) + t112; t112; (-g(1) * t79 + g(2) * t77 - t110 * t127) * MDP(12) + (g(1) * t80 - g(2) * t78 + t108 * t127) * MDP(13) + t120; t120; t126;];
taug = t1;
