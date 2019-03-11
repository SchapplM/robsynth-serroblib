% Calculate Gravitation load on the joints for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRPRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:19:46
% EndTime: 2019-03-09 03:19:47
% DurationCPUTime: 0.47s
% Computational Cost: add. (217->70), mult. (274->92), div. (0->0), fcn. (233->8), ass. (0->42)
t83 = pkin(9) + qJ(3);
t80 = sin(t83);
t81 = cos(t83);
t101 = t81 * pkin(3) + t80 * qJ(4);
t85 = cos(pkin(9));
t121 = pkin(2) * t85 + pkin(1) + t101;
t91 = cos(qJ(1));
t120 = g(2) * t91;
t113 = g(1) * t91;
t89 = sin(qJ(1));
t73 = g(2) * t89 + t113;
t119 = MDP(14) - MDP(17);
t61 = g(3) * t80 + t73 * t81;
t117 = MDP(13) - MDP(16) + MDP(26);
t116 = pkin(3) * t80;
t88 = sin(qJ(5));
t115 = pkin(5) * t88;
t109 = g(3) * t81;
t86 = -qJ(6) - pkin(8);
t107 = t81 * t86;
t106 = t81 * t88;
t105 = t88 * t91;
t104 = t89 * t88;
t90 = cos(qJ(5));
t103 = t89 * t90;
t102 = t90 * t91;
t87 = -pkin(7) - qJ(2);
t100 = pkin(5) * t90 + pkin(4) - t87;
t99 = qJ(4) * t81;
t98 = -MDP(18) - MDP(27);
t97 = pkin(5) * t106;
t95 = t121 * t120;
t72 = g(1) * t89 - t120;
t94 = t80 * t115 - t107;
t64 = t80 * t102 - t104;
t66 = t80 * t103 + t105;
t70 = t91 * t99;
t68 = t89 * t99;
t67 = -t80 * t104 + t102;
t65 = t80 * t105 + t103;
t60 = t73 * t80 - t109;
t1 = [(-g(1) * (-t89 * pkin(1) + qJ(2) * t91) - g(2) * (pkin(1) * t91 + t89 * qJ(2))) * MDP(7) + (t87 * t113 - t95 + (g(1) * t121 + g(2) * t87) * t89) * MDP(18) + (-g(1) * t67 - g(2) * t65) * MDP(24) + (g(1) * t66 - g(2) * t64) * MDP(25) + (-t95 + (-g(1) * t100 - g(2) * t94) * t91 + (-g(1) * (-t121 - t94) - g(2) * t100) * t89) * MDP(27) + (MDP(3) - MDP(6) - MDP(15)) * t73 + (t85 * MDP(4) - MDP(5) * sin(pkin(9)) + t117 * t81 - t119 * t80 + MDP(2)) * t72; (-MDP(7) + t98) * t72; (-g(1) * (-t91 * t116 + t70) - g(2) * (-t89 * t116 + t68) - g(3) * t101) * MDP(18) + (-g(1) * (t91 * t97 + t70) - g(2) * (t89 * t97 + t68) - g(3) * (t101 - t107) + (-g(3) * t115 + t73 * (pkin(3) - t86)) * t80) * MDP(27) + t117 * t60 + (-t88 * MDP(24) - MDP(25) * t90 + t119) * t61; t98 * t60; (g(1) * t65 - g(2) * t67 - g(3) * t106) * MDP(25) + (pkin(5) * MDP(27) + MDP(24)) * (-g(1) * t64 - g(2) * t66 + t90 * t109); -t61 * MDP(27);];
taug  = t1;
