% Calculate Gravitation load on the joints for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:16:44
% EndTime: 2019-03-08 20:16:45
% DurationCPUTime: 0.58s
% Computational Cost: add. (217->81), mult. (530->130), div. (0->0), fcn. (607->10), ass. (0->45)
t114 = MDP(14) - MDP(22);
t79 = sin(pkin(10));
t81 = cos(pkin(10));
t88 = cos(qJ(2));
t85 = sin(qJ(2));
t99 = cos(pkin(6));
t96 = t85 * t99;
t68 = t79 * t88 + t81 * t96;
t70 = -t79 * t96 + t81 * t88;
t113 = -g(1) * t70 - g(2) * t68;
t80 = sin(pkin(6));
t110 = g(3) * t80;
t95 = t88 * t99;
t67 = t79 * t85 - t81 * t95;
t83 = sin(qJ(5));
t109 = t67 * t83;
t69 = t79 * t95 + t81 * t85;
t108 = t69 * t83;
t84 = sin(qJ(4));
t107 = t80 * t84;
t106 = t80 * t85;
t87 = cos(qJ(4));
t105 = t80 * t87;
t104 = t80 * t88;
t103 = t83 * t84;
t102 = t83 * t85;
t86 = cos(qJ(5));
t101 = t84 * t86;
t100 = t85 * t86;
t98 = MDP(23) + MDP(7);
t97 = g(3) * (pkin(2) * t104 + qJ(3) * t106);
t78 = t86 * pkin(5) + pkin(4);
t82 = -qJ(6) - pkin(9);
t94 = t78 * t84 + t82 * t87;
t59 = t79 * t107 - t69 * t87;
t61 = t81 * t107 + t67 * t87;
t71 = t87 * t104 + t99 * t84;
t92 = g(1) * t59 - g(2) * t61 + g(3) * t71;
t58 = -g(1) * t69 - g(2) * t67 + g(3) * t104;
t72 = -t84 * t104 + t99 * t87;
t66 = t69 * pkin(2);
t65 = t67 * pkin(2);
t62 = t81 * t105 - t67 * t84;
t60 = t79 * t105 + t69 * t84;
t1 = [(-MDP(1) - t98) * g(3); (-g(1) * (t70 * qJ(3) - t66) - g(2) * (t68 * qJ(3) - t65) - t97) * MDP(7) + (-g(1) * (t70 * t101 - t108) - g(2) * (t68 * t101 - t109) - (t84 * t100 + t83 * t88) * t110) * MDP(20) + (-g(1) * (-t70 * t103 - t69 * t86) - g(2) * (-t68 * t103 - t67 * t86) - (-t84 * t102 + t86 * t88) * t110) * MDP(21) + (-g(1) * (-pkin(5) * t108 - t69 * pkin(8) - t66) - g(2) * (-pkin(5) * t109 - t67 * pkin(8) - t65) - t97 - ((pkin(5) * t83 + pkin(8)) * t88 + t94 * t85) * t110 + t113 * (qJ(3) + t94)) * MDP(23) + (-MDP(3) + MDP(5)) * t58 + (-t84 * MDP(13) - t114 * t87 + MDP(4) - MDP(6)) * (g(3) * t106 - t113); t98 * t58; (-g(1) * (-t59 * t78 - t60 * t82) - g(2) * (t61 * t78 + t62 * t82) - g(3) * (-t71 * t78 - t72 * t82)) * MDP(23) + t114 * (g(1) * t60 - g(2) * t62 + g(3) * t72) + (MDP(20) * t86 - MDP(21) * t83 + MDP(13)) * t92; (-g(1) * (-t60 * t86 - t70 * t83) - g(2) * (t62 * t86 - t68 * t83) - g(3) * (-t80 * t102 - t72 * t86)) * MDP(21) + (pkin(5) * MDP(23) + MDP(20)) * (-g(1) * (-t60 * t83 + t70 * t86) - g(2) * (t62 * t83 + t68 * t86) - g(3) * (t80 * t100 - t72 * t83)); -t92 * MDP(23);];
taug  = t1;
