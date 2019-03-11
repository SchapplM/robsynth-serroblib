% Calculate Gravitation load on the joints for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRPRP9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:29:38
% EndTime: 2019-03-09 03:29:39
% DurationCPUTime: 0.64s
% Computational Cost: add. (250->87), mult. (374->117), div. (0->0), fcn. (349->8), ass. (0->40)
t87 = sin(qJ(1));
t112 = g(1) * t87;
t89 = cos(qJ(1));
t118 = g(2) * t89 - t112;
t117 = MDP(23) + MDP(25);
t116 = MDP(24) - MDP(27);
t86 = sin(qJ(3));
t88 = cos(qJ(3));
t68 = -g(3) * t86 - t118 * t88;
t114 = MDP(13) - MDP(16) - MDP(26);
t113 = -pkin(1) - pkin(7);
t109 = g(3) * t88;
t83 = sin(pkin(9));
t108 = t83 * t89;
t106 = t86 * t89;
t82 = pkin(9) + qJ(5);
t75 = sin(t82);
t105 = t87 * t75;
t76 = cos(t82);
t104 = t87 * t76;
t103 = t87 * t83;
t84 = cos(pkin(9));
t102 = t87 * t84;
t85 = -pkin(8) - qJ(4);
t101 = t88 * t85;
t100 = t89 * t76;
t99 = pkin(1) * t89 + qJ(2) * t87;
t98 = qJ(4) * t88;
t97 = MDP(17) + MDP(28);
t96 = pkin(7) * t89 + t99;
t93 = pkin(3) * t86 - t98;
t74 = pkin(4) * t84 + pkin(3);
t92 = pkin(5) * t76 + qJ(6) * t75 + t74;
t63 = t105 * t86 - t100;
t65 = t106 * t75 + t104;
t59 = g(1) * t63 - g(2) * t65 + t109 * t75;
t78 = t89 * qJ(2);
t66 = t100 * t86 - t105;
t64 = t104 * t86 + t75 * t89;
t1 = [(-g(1) * (-t87 * pkin(1) + t78) - g(2) * t99) * MDP(6) + (-g(1) * (t106 * t84 - t103) - g(2) * (t102 * t86 + t108)) * MDP(14) + (-g(1) * (-t106 * t83 - t102) - g(2) * (-t103 * t86 + t84 * t89)) * MDP(15) + (-g(1) * (pkin(3) * t106 - t89 * t98 + t78) - g(2) * t96 + (-g(1) * t113 - g(2) * t93) * t87) * MDP(17) + (-g(1) * (t66 * pkin(5) + t65 * qJ(6) + t89 * t101 + t74 * t106 + t78) - g(2) * (pkin(4) * t108 + t64 * pkin(5) + t63 * qJ(6) + t96) + (-g(1) * (-pkin(4) * t83 + t113) - g(2) * (t74 * t86 + t101)) * t87) * MDP(28) + t116 * (g(1) * t65 + g(2) * t63) - (MDP(2) - MDP(4)) * t118 + t117 * (-g(1) * t66 - g(2) * t64) + (-MDP(12) * t86 - t114 * t88 + MDP(3) - MDP(5)) * (g(1) * t89 + g(2) * t87); -(-MDP(6) - t97) * t118; (g(3) * t93 + t118 * (pkin(3) * t88 + qJ(4) * t86)) * MDP(17) + (-g(3) * (-t86 * t92 - t101) + t118 * (-t85 * t86 + t88 * t92)) * MDP(28) + t114 * (-g(2) * t106 + t112 * t86 + t109) + (-t84 * MDP(14) + t83 * MDP(15) + t116 * t75 - t117 * t76 - MDP(12)) * t68; t97 * t68; (-g(1) * (-pkin(5) * t63 + qJ(6) * t64) - g(2) * (pkin(5) * t65 - qJ(6) * t66) - (-pkin(5) * t75 + qJ(6) * t76) * t109) * MDP(28) + t116 * (g(1) * t64 - g(2) * t66 + t109 * t76) + t117 * t59; -t59 * MDP(28);];
taug  = t1;
