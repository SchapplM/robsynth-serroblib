% Calculate Gravitation load on the joints for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRRRP8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:24:34
% EndTime: 2019-03-09 06:24:35
% DurationCPUTime: 0.41s
% Computational Cost: add. (277->76), mult. (387->100), div. (0->0), fcn. (359->8), ass. (0->40)
t90 = sin(qJ(5));
t124 = qJ(6) * t90 + pkin(4);
t123 = MDP(20) - MDP(29);
t121 = MDP(26) + MDP(28);
t120 = MDP(27) - MDP(30);
t89 = qJ(3) + qJ(4);
t84 = cos(t89);
t81 = t84 * pkin(9);
t91 = sin(qJ(3));
t122 = -pkin(3) * t91 + t81;
t94 = cos(qJ(3));
t117 = pkin(3) * t94;
t83 = sin(t89);
t116 = pkin(9) * t83;
t95 = cos(qJ(1));
t115 = g(2) * t95;
t114 = g(3) * t83;
t113 = g(3) * t84;
t92 = sin(qJ(1));
t111 = t90 * t92;
t110 = t90 * t95;
t93 = cos(qJ(5));
t109 = t92 * t93;
t108 = t95 * t93;
t107 = t95 * pkin(1) + t92 * qJ(2);
t105 = t92 * t116 + (pkin(5) * t109 + t124 * t92) * t84;
t77 = g(1) * t92 - t115;
t103 = -pkin(5) * t93 - t124;
t102 = t123 * (t77 * t83 + t113) + (t120 * t90 - t121 * t93 - MDP(19)) * (t77 * t84 - t114);
t101 = pkin(4) * t83 - t122;
t71 = t83 * t111 - t108;
t73 = t83 * t110 + t109;
t58 = g(1) * t71 - g(2) * t73 + t90 * t113;
t99 = t103 * t114;
t97 = t103 * t84 - t116;
t96 = -pkin(8) - pkin(7);
t86 = t95 * qJ(2);
t74 = t83 * t108 - t111;
t72 = t83 * t109 + t110;
t1 = [(-g(1) * (-t92 * pkin(1) + t86) - g(2) * t107) * MDP(6) + (-g(1) * (t74 * pkin(5) + t73 * qJ(6) + t86) - g(2) * (t72 * pkin(5) + t71 * qJ(6) + t107) + (-g(1) * t101 + g(2) * t96) * t95 + (-g(1) * (-pkin(1) + t96) - g(2) * t101) * t92) * MDP(31) + (MDP(2) - MDP(4)) * t77 + t121 * (-g(1) * t74 - g(2) * t72) + t120 * (g(1) * t73 + g(2) * t71) + (-t91 * MDP(12) - t94 * MDP(13) - t83 * MDP(19) - t123 * t84 + MDP(3) - MDP(5)) * (g(1) * t95 + g(2) * t92); (-MDP(31) - MDP(6)) * t77; (g(3) * t91 - t77 * t94) * MDP(12) + (g(3) * t94 + t77 * t91) * MDP(13) + (-g(1) * (t92 * t117 + t105) - g(3) * t122 - t99 - (t97 - t117) * t115) * MDP(31) + t102; (-g(1) * t105 - g(3) * t81 - t97 * t115 - t99) * MDP(31) + t102; (-g(1) * (-pkin(5) * t71 + qJ(6) * t72) - g(2) * (pkin(5) * t73 - qJ(6) * t74) - (-pkin(5) * t90 + qJ(6) * t93) * t113) * MDP(31) + t121 * t58 + t120 * (g(1) * t72 - g(2) * t74 + t93 * t113); -t58 * MDP(31);];
taug  = t1;
