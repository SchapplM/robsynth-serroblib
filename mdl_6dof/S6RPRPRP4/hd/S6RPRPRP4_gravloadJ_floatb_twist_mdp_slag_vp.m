% Calculate Gravitation load on the joints for
% S6RPRPRP4
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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:13:07
% EndTime: 2019-03-09 03:13:08
% DurationCPUTime: 0.44s
% Computational Cost: add. (272->76), mult. (340->106), div. (0->0), fcn. (314->8), ass. (0->38)
t117 = MDP(11) - MDP(14);
t114 = MDP(21) + MDP(23);
t113 = MDP(22) - MDP(25);
t116 = MDP(10) - MDP(13) + MDP(24);
t85 = sin(qJ(3));
t79 = t85 * qJ(4);
t88 = cos(qJ(3));
t102 = t88 * pkin(3) + t79;
t83 = qJ(1) + pkin(9);
t77 = sin(t83);
t78 = cos(t83);
t115 = -g(1) * t78 - g(2) * t77;
t112 = pkin(3) * t85;
t111 = g(1) * t77;
t107 = g(3) * t88;
t106 = t88 * pkin(8);
t105 = t78 * t88;
t84 = sin(qJ(5));
t104 = t84 * t85;
t87 = cos(qJ(5));
t103 = t85 * t87;
t101 = qJ(4) * t88;
t100 = -MDP(15) - MDP(26);
t86 = sin(qJ(1));
t99 = -t86 * pkin(1) + t78 * pkin(7);
t89 = cos(qJ(1));
t98 = t89 * pkin(1) + pkin(3) * t105 + t77 * pkin(7) + (pkin(2) + t79) * t78;
t93 = pkin(5) * t84 - qJ(6) * t87;
t92 = -pkin(2) - t102;
t63 = -t103 * t78 + t77 * t84;
t65 = t103 * t77 + t78 * t84;
t57 = g(1) * t63 - g(2) * t65 + t87 * t107;
t72 = t78 * t101;
t70 = t77 * t101;
t66 = -t104 * t77 + t78 * t87;
t64 = t104 * t78 + t77 * t87;
t61 = -t115 * t85 - t107;
t1 = [(g(1) * t89 + g(2) * t86) * MDP(3) + t115 * MDP(12) + (-g(1) * t99 - g(2) * t98 - t111 * t92) * MDP(15) + (-g(1) * (t78 * pkin(4) + t66 * pkin(5) + t65 * qJ(6) + t99) - g(2) * (t64 * pkin(5) + pkin(8) * t105 + t63 * qJ(6) + t98) + (-g(1) * (t92 - t106) - g(2) * pkin(4)) * t77) * MDP(26) + t113 * (g(1) * t65 + g(2) * t63) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t86 - g(2) * t89) + t114 * (-g(1) * t66 - g(2) * t64) + (t116 * t88 - t117 * t85) * (-g(2) * t78 + t111); (-MDP(4) + t100) * g(3); (-g(1) * (-t112 * t78 + t72) - g(2) * (-t112 * t77 + t70) - g(3) * t102) * MDP(15) + (-g(1) * t72 - g(2) * t70 - g(3) * (t85 * t93 + t102 + t106) + t115 * (t93 * t88 + (-pkin(3) - pkin(8)) * t85)) * MDP(26) + t116 * t61 + (-t113 * t87 - t114 * t84 + t117) * (g(3) * t85 - t115 * t88); t100 * t61; (-g(1) * (-t63 * pkin(5) + t64 * qJ(6)) - g(2) * (t65 * pkin(5) - t66 * qJ(6)) - (-pkin(5) * t87 - qJ(6) * t84) * t107) * MDP(26) - t113 * (-g(1) * t64 + g(2) * t66 + t84 * t107) + t114 * t57; -t57 * MDP(26);];
taug  = t1;
