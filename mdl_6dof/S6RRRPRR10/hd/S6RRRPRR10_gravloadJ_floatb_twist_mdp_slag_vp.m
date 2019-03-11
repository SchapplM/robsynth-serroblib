% Calculate Gravitation load on the joints for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:21:58
% EndTime: 2019-03-09 19:21:59
% DurationCPUTime: 0.62s
% Computational Cost: add. (316->86), mult. (647->130), div. (0->0), fcn. (729->10), ass. (0->41)
t105 = sin(qJ(5));
t106 = sin(qJ(3));
t109 = cos(qJ(5));
t110 = cos(qJ(3));
t117 = t105 * t106 + t109 * t110;
t118 = t105 * t110 - t106 * t109;
t108 = sin(qJ(1));
t111 = cos(qJ(2));
t112 = cos(qJ(1));
t131 = t112 * t106;
t97 = -t108 * t110 + t111 * t131;
t130 = t112 * t110;
t98 = t108 * t106 + t111 * t130;
t121 = t98 * t105 - t97 * t109;
t132 = t108 * t111;
t95 = t106 * t132 + t130;
t96 = t110 * t132 - t131;
t122 = t95 * t105 + t96 * t109;
t128 = t96 * t105 - t95 * t109;
t104 = qJ(5) + qJ(6);
t102 = sin(t104);
t103 = cos(t104);
t119 = t102 * t106 + t103 * t110;
t120 = t102 * t110 - t103 * t106;
t123 = t98 * t102 - t97 * t103;
t124 = t95 * t102 + t96 * t103;
t129 = t96 * t102 - t95 * t103;
t107 = sin(qJ(2));
t137 = g(3) * t107;
t81 = t97 * t102 + t98 * t103;
t135 = (g(1) * t123 + g(2) * t129 + t120 * t137) * MDP(34) + (g(1) * t81 + g(2) * t124 + t119 * t137) * MDP(35);
t84 = t97 * t105 + t98 * t109;
t162 = t135 + (g(1) * t121 + g(2) * t128 + t118 * t137) * MDP(27) + (g(1) * t84 + g(2) * t122 + t117 * t137) * MDP(28);
t157 = MDP(10) - MDP(19);
t151 = MDP(16) + MDP(18);
t150 = MDP(17) - MDP(20);
t126 = g(1) * t112 + g(2) * t108;
t116 = t111 * pkin(2) + t107 * pkin(8) + pkin(1);
t115 = pkin(3) * t110 + qJ(4) * t106 + pkin(2);
t82 = g(1) * t97 + g(2) * t95 + t106 * t137;
t1 = [t126 * MDP(3) + (-g(1) * (-t96 * pkin(3) - t95 * qJ(4)) - g(2) * (t98 * pkin(3) + t97 * qJ(4)) + (-g(1) * pkin(7) - g(2) * t116) * t112 + (-g(2) * pkin(7) + g(1) * t116) * t108) * MDP(21) + (g(1) * t122 - g(2) * t84) * MDP(27) + (-g(1) * t128 + g(2) * t121) * MDP(28) + (g(1) * t124 - g(2) * t81) * MDP(34) + (-g(1) * t129 + g(2) * t123) * MDP(35) + t151 * (g(1) * t96 - g(2) * t98) - t150 * (g(1) * t95 - g(2) * t97) + (t111 * MDP(9) - t157 * t107 + MDP(2)) * (g(1) * t108 - g(2) * t112); ((-t126 * pkin(8) - g(3) * t115) * t111 + (-g(3) * pkin(8) + t115 * t126) * t107) * MDP(21) + t157 * (t126 * t111 + t137) + (-t117 * MDP(27) + t118 * MDP(28) - t119 * MDP(34) + t120 * MDP(35) + t150 * t106 - t151 * t110 - MDP(9)) * (g(3) * t111 - t126 * t107); (-g(1) * (-t97 * pkin(3) + t98 * qJ(4)) - g(2) * (-t95 * pkin(3) + t96 * qJ(4)) - (-pkin(3) * t106 + qJ(4) * t110) * t137) * MDP(21) + t151 * t82 + t150 * (g(1) * t98 + g(2) * t96 + t110 * t137) - t162; -t82 * MDP(21); t162; t135;];
taug  = t1;
