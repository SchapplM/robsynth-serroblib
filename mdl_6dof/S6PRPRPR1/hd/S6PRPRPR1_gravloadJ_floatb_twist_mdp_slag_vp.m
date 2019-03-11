% Calculate Gravitation load on the joints for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:27:58
% EndTime: 2019-03-08 19:28:00
% DurationCPUTime: 0.57s
% Computational Cost: add. (270->87), mult. (614->157), div. (0->0), fcn. (750->14), ass. (0->47)
t110 = cos(pkin(11));
t88 = sin(pkin(11));
t96 = sin(qJ(2));
t99 = cos(qJ(2));
t105 = t99 * t110 - t96 * t88;
t87 = qJ(4) + pkin(12);
t86 = cos(t87);
t94 = sin(qJ(6));
t121 = t86 * t94;
t97 = cos(qJ(6));
t120 = t86 * t97;
t89 = sin(pkin(10));
t90 = sin(pkin(6));
t119 = t89 * t90;
t118 = t89 * t96;
t91 = cos(pkin(10));
t117 = t90 * t91;
t95 = sin(qJ(4));
t116 = t90 * t95;
t98 = cos(qJ(4));
t115 = t90 * t98;
t114 = t90 * t99;
t92 = cos(pkin(6));
t113 = t92 * t96;
t112 = t92 * t99;
t109 = MDP(14) + MDP(5);
t108 = t91 * t112;
t80 = -t96 * t110 - t99 * t88;
t78 = t80 * t92;
t67 = t105 * t89 - t78 * t91;
t68 = -t105 * t91 - t89 * t78;
t106 = -t89 * t112 - t91 * t96;
t102 = t92 * t105;
t66 = t91 * t102 + t89 * t80;
t69 = -t89 * t102 + t80 * t91;
t76 = t105 * t90;
t103 = g(1) * t69 + g(2) * t66 + g(3) * t76;
t101 = -g(1) * t106 - g(3) * t114;
t93 = -qJ(5) - pkin(8);
t85 = sin(t87);
t84 = pkin(4) * t98 + pkin(3);
t81 = pkin(2) * t108;
t77 = t80 * t90;
t72 = -t77 * t86 + t85 * t92;
t64 = t85 * t119 - t68 * t86;
t62 = -t85 * t117 + t67 * t86;
t1 = [(-MDP(1) - t109) * g(3); (-g(2) * (t108 - t118) + t101) * MDP(3) + (-g(1) * (t89 * t113 - t91 * t99) - g(2) * (-t91 * t113 - t89 * t99) + g(3) * t90 * t96) * MDP(4) + (-g(2) * t81 + (g(2) * t118 + t101) * pkin(2)) * MDP(5) + (g(1) * t68 - g(2) * t67 + g(3) * t77) * MDP(13) + (-g(1) * (t106 * pkin(2) + t68 * t93 + t69 * t84) - g(2) * (-pkin(2) * t118 + t66 * t84 - t67 * t93 + t81) - g(3) * (pkin(2) * t114 + t76 * t84 + t77 * t93)) * MDP(14) + (-g(1) * (t69 * t120 - t68 * t94) - g(2) * (t66 * t120 + t67 * t94) - g(3) * (t76 * t120 - t77 * t94)) * MDP(20) + (-g(1) * (-t69 * t121 - t68 * t97) - g(2) * (-t66 * t121 + t67 * t97) - g(3) * (-t76 * t121 - t77 * t97)) * MDP(21) + (-MDP(11) * t98 + MDP(12) * t95) * t103; t109 * (-g(3) * t92 + (-g(1) * t89 + g(2) * t91) * t90); (-g(1) * (-t89 * t116 + t68 * t98) - g(2) * (t91 * t116 - t67 * t98) - g(3) * (t77 * t98 - t92 * t95)) * MDP(12) + (-MDP(20) * t97 + MDP(21) * t94) * (g(1) * (t86 * t119 + t68 * t85) + g(2) * (-t86 * t117 - t67 * t85) + g(3) * (t77 * t85 + t86 * t92)) + (pkin(4) * MDP(14) + MDP(11)) * (-g(1) * (t89 * t115 + t68 * t95) - g(2) * (-t91 * t115 - t67 * t95) - g(3) * (t77 * t95 + t92 * t98)); t103 * MDP(14); (-g(1) * (-t64 * t94 - t69 * t97) - g(2) * (-t62 * t94 - t66 * t97) - g(3) * (-t72 * t94 - t76 * t97)) * MDP(20) + (-g(1) * (-t64 * t97 + t69 * t94) - g(2) * (-t62 * t97 + t66 * t94) - g(3) * (-t72 * t97 + t76 * t94)) * MDP(21);];
taug  = t1;
