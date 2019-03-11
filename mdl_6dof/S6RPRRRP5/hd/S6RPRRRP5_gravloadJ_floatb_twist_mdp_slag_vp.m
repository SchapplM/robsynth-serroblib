% Calculate Gravitation load on the joints for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:12:25
% EndTime: 2019-03-09 06:12:26
% DurationCPUTime: 0.40s
% Computational Cost: add. (405->76), mult. (392->102), div. (0->0), fcn. (363->10), ass. (0->40)
t101 = cos(qJ(5));
t99 = sin(qJ(5));
t127 = pkin(5) * t101 + qJ(6) * t99;
t126 = MDP(21) - MDP(30);
t123 = MDP(27) + MDP(29);
t122 = MDP(28) - MDP(31);
t100 = sin(qJ(1));
t102 = cos(qJ(1));
t83 = g(1) * t102 + g(2) * t100;
t96 = pkin(10) + qJ(3);
t93 = qJ(4) + t96;
t87 = sin(t93);
t125 = t83 * t87;
t88 = cos(t93);
t124 = t88 * pkin(4) + t87 * pkin(9);
t91 = sin(t96);
t121 = pkin(3) * t91;
t120 = pkin(9) * t88;
t119 = g(3) * t87;
t114 = t100 * t99;
t113 = t102 * t99;
t112 = t100 * t101;
t111 = t101 * t102;
t110 = t127 * t88 + t124;
t82 = g(1) * t100 - g(2) * t102;
t92 = cos(t96);
t86 = pkin(3) * t92;
t98 = cos(pkin(10));
t108 = pkin(2) * t98 + pkin(1) + t124 + t86;
t106 = t126 * (t83 * t88 + t119) + (t123 * t101 - t122 * t99 + MDP(20)) * (-g(3) * t88 + t125);
t73 = t88 * t114 + t111;
t75 = t88 * t113 - t112;
t60 = g(1) * t75 + g(2) * t73 + t99 * t119;
t103 = (pkin(4) + t127) * t125;
t95 = -pkin(8) - pkin(7) - qJ(2);
t81 = t102 * t120;
t79 = t100 * t120;
t76 = t88 * t111 + t114;
t74 = t88 * t112 - t113;
t1 = [(-g(1) * (-t100 * pkin(1) + qJ(2) * t102) - g(2) * (pkin(1) * t102 + t100 * qJ(2))) * MDP(7) + (-g(1) * (-t74 * pkin(5) - t73 * qJ(6)) - g(2) * (t76 * pkin(5) + t75 * qJ(6)) + (g(1) * t95 - g(2) * t108) * t102 + (g(1) * t108 + g(2) * t95) * t100) * MDP(32) + (MDP(3) - MDP(6)) * t83 + t123 * (g(1) * t74 - g(2) * t76) - t122 * (g(1) * t73 - g(2) * t75) + (-t126 * t87 + t92 * MDP(13) - t91 * MDP(14) + t88 * MDP(20) + MDP(4) * t98 - MDP(5) * sin(pkin(10)) + MDP(2)) * t82; (-MDP(32) - MDP(7)) * t82; (-g(3) * t92 + t83 * t91) * MDP(13) + (g(3) * t91 + t83 * t92) * MDP(14) + (-g(1) * (-t102 * t121 + t81) - g(2) * (-t100 * t121 + t79) - g(3) * (t86 + t110) + t103) * MDP(32) + t106; (-g(1) * t81 - g(2) * t79 - g(3) * t110 + t103) * MDP(32) + t106; (-g(1) * (-pkin(5) * t75 + qJ(6) * t76) - g(2) * (-pkin(5) * t73 + qJ(6) * t74) - (-pkin(5) * t99 + qJ(6) * t101) * t119) * MDP(32) + t123 * t60 + t122 * (g(1) * t76 + g(2) * t74 + t101 * t119); -t60 * MDP(32);];
taug  = t1;
