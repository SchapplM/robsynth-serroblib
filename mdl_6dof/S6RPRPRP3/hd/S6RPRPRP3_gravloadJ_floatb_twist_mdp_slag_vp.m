% Calculate Gravitation load on the joints for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:09:37
% EndTime: 2019-03-09 03:09:38
% DurationCPUTime: 0.48s
% Computational Cost: add. (363->82), mult. (360->118), div. (0->0), fcn. (337->10), ass. (0->39)
t114 = MDP(21) + MDP(23);
t113 = MDP(22) - MDP(25);
t116 = MDP(11) - MDP(14) - MDP(24);
t84 = -pkin(8) - qJ(4);
t85 = sin(qJ(3));
t102 = t85 * t84;
t83 = cos(pkin(10));
t73 = t83 * pkin(4) + pkin(3);
t87 = cos(qJ(3));
t115 = t87 * t73 - t102;
t81 = qJ(1) + pkin(9);
t75 = sin(t81);
t77 = cos(t81);
t96 = g(1) * t77 + g(2) * t75;
t65 = -g(3) * t87 + t96 * t85;
t112 = g(1) * t75;
t109 = g(3) * t85;
t107 = t75 * t87;
t82 = sin(pkin(10));
t106 = t77 * t82;
t105 = t77 * t87;
t104 = t82 * t87;
t103 = t83 * t87;
t100 = -MDP(15) - MDP(26);
t88 = cos(qJ(1));
t99 = t88 * pkin(1) + t77 * pkin(2) + t75 * pkin(7);
t86 = sin(qJ(1));
t98 = -t86 * pkin(1) + t77 * pkin(7);
t93 = t87 * pkin(3) + t85 * qJ(4);
t80 = pkin(10) + qJ(5);
t74 = sin(t80);
t76 = cos(t80);
t91 = pkin(5) * t76 + qJ(6) * t74 + t73;
t61 = t74 * t107 + t77 * t76;
t63 = t74 * t105 - t75 * t76;
t57 = g(1) * t63 + g(2) * t61 + t74 * t109;
t64 = t76 * t105 + t75 * t74;
t62 = t76 * t107 - t77 * t74;
t1 = [(g(1) * t88 + g(2) * t86) * MDP(3) + (-g(1) * (-t75 * t103 + t106) - g(2) * (t77 * t103 + t75 * t82)) * MDP(12) + (-g(1) * (t75 * t104 + t77 * t83) - g(2) * (-t77 * t104 + t75 * t83)) * MDP(13) + (-g(1) * t98 - g(2) * (t93 * t77 + t99) - (-pkin(2) - t93) * t112) * MDP(15) + (-g(1) * (pkin(4) * t106 - t62 * pkin(5) - t61 * qJ(6) + t98) - g(2) * (t64 * pkin(5) + t63 * qJ(6) + t115 * t77 + t99) + (-g(1) * (-pkin(2) - t115) - g(2) * pkin(4) * t82) * t75) * MDP(26) - t113 * (g(1) * t61 - g(2) * t63) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t86 - g(2) * t88) + t114 * (g(1) * t62 - g(2) * t64) + (t87 * MDP(10) - t116 * t85) * (-g(2) * t77 + t112); (-MDP(4) + t100) * g(3); (-g(3) * t93 + t96 * (pkin(3) * t85 - qJ(4) * t87)) * MDP(15) + (-g(3) * (t91 * t87 - t102) + t96 * (t84 * t87 + t91 * t85)) * MDP(26) + t116 * (t96 * t87 + t109) + (t83 * MDP(12) - t82 * MDP(13) - t113 * t74 + t114 * t76 + MDP(10)) * t65; t100 * t65; (-g(1) * (-t63 * pkin(5) + t64 * qJ(6)) - g(2) * (-t61 * pkin(5) + t62 * qJ(6)) - (-pkin(5) * t74 + qJ(6) * t76) * t109) * MDP(26) + t113 * (g(1) * t64 + g(2) * t62 + t76 * t109) + t114 * t57; -t57 * MDP(26);];
taug  = t1;
