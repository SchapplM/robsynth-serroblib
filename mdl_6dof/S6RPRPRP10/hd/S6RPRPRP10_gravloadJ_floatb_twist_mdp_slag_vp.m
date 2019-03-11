% Calculate Gravitation load on the joints for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRPRP10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:32:32
% EndTime: 2019-03-09 03:32:34
% DurationCPUTime: 0.46s
% Computational Cost: add. (161->81), mult. (354->107), div. (0->0), fcn. (326->6), ass. (0->36)
t85 = sin(qJ(5));
t88 = cos(qJ(5));
t94 = pkin(5) * t85 - qJ(6) * t88;
t116 = MDP(13) - MDP(16);
t115 = MDP(23) + MDP(25);
t114 = MDP(24) - MDP(27);
t113 = -MDP(12) + MDP(15) - MDP(26);
t112 = -pkin(1) - pkin(7);
t111 = -pkin(3) - pkin(8);
t87 = sin(qJ(1));
t109 = g(1) * t87;
t90 = cos(qJ(1));
t108 = g(2) * t90;
t86 = sin(qJ(3));
t107 = g(3) * t86;
t105 = t86 * pkin(3);
t104 = t86 * t90;
t89 = cos(qJ(3));
t103 = t87 * t89;
t102 = t89 * t90;
t99 = qJ(4) * t86;
t101 = pkin(3) * t103 + t87 * t99;
t100 = t90 * pkin(1) + t87 * qJ(2);
t79 = t89 * qJ(4);
t97 = MDP(17) + MDP(28);
t96 = t90 * pkin(7) + t87 * t105 + t100;
t69 = -t108 + t109;
t80 = t90 * qJ(2);
t93 = pkin(3) * t104 - t90 * t79 + t80;
t63 = -t88 * t102 + t85 * t87;
t65 = t88 * t103 + t85 * t90;
t92 = -g(1) * t65 - g(2) * t63 + t88 * t107;
t66 = -t85 * t103 + t88 * t90;
t64 = t85 * t102 + t87 * t88;
t62 = g(1) * t103 - g(2) * t102 - t107;
t1 = [(-g(1) * (-t87 * pkin(1) + t80) - g(2) * t100) * MDP(6) + (-g(1) * (t112 * t87 + t93) - g(2) * (-t87 * t79 + t96)) * MDP(17) + (-g(1) * (-t64 * pkin(5) + pkin(8) * t104 - t63 * qJ(6) + t93) - g(2) * (pkin(4) * t90 + t66 * pkin(5) + t65 * qJ(6) + t96) + (-g(1) * (-pkin(4) + t112) - g(2) * (t86 * pkin(8) - t79)) * t87) * MDP(28) - t114 * (g(1) * t63 - g(2) * t65) + t115 * (g(1) * t64 - g(2) * t66) + (MDP(2) - MDP(4) + MDP(14)) * t69 + (t113 * t86 - t116 * t89 + MDP(3) - MDP(5)) * (g(1) * t90 + g(2) * t87); (-MDP(6) - t97) * t69; (-g(1) * t101 - g(3) * (t79 - t105) - (-pkin(3) * t89 - t99) * t108) * MDP(17) + (-g(1) * (pkin(8) * t103 + t101) - g(3) * (t94 * t89 + t79) + (-g(3) * t111 - t94 * t109) * t86 - (t111 * t89 + (-qJ(4) - t94) * t86) * t108) * MDP(28) + t113 * t62 + (-t114 * t88 - t115 * t85 + t116) * (g(3) * t89 + t69 * t86); t97 * t62; (-g(1) * (-pkin(5) * t65 + qJ(6) * t66) - g(2) * (-pkin(5) * t63 + qJ(6) * t64) - (pkin(5) * t88 + qJ(6) * t85) * t107) * MDP(28) + t114 * (g(1) * t66 + g(2) * t64 + t85 * t107) - t115 * t92; t92 * MDP(28);];
taug  = t1;
