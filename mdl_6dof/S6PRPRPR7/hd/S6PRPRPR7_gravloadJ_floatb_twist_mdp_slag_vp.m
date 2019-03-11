% Calculate Gravitation load on the joints for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:53:48
% EndTime: 2019-03-08 19:53:49
% DurationCPUTime: 0.57s
% Computational Cost: add. (205->77), mult. (517->123), div. (0->0), fcn. (595->10), ass. (0->40)
t107 = MDP(13) - MDP(16);
t106 = -MDP(14) + MDP(17);
t77 = sin(pkin(10));
t79 = cos(pkin(10));
t85 = cos(qJ(2));
t82 = sin(qJ(2));
t95 = cos(pkin(6));
t92 = t82 * t95;
t67 = t77 * t85 + t79 * t92;
t69 = -t77 * t92 + t79 * t85;
t105 = -g(1) * t69 - g(2) * t67;
t78 = sin(pkin(6));
t102 = g(3) * t78;
t81 = sin(qJ(4));
t101 = t78 * t81;
t100 = t78 * t82;
t84 = cos(qJ(4));
t99 = t78 * t84;
t98 = t78 * t85;
t80 = sin(qJ(6));
t97 = t80 * t84;
t83 = cos(qJ(6));
t96 = t83 * t84;
t94 = MDP(18) + MDP(7);
t93 = g(3) * (pkin(2) * t98 + qJ(3) * t100);
t91 = t85 * t95;
t90 = pkin(4) * t81 - qJ(5) * t84;
t68 = t77 * t91 + t79 * t82;
t58 = t77 * t101 - t68 * t84;
t66 = t77 * t82 - t79 * t91;
t60 = t79 * t101 + t66 * t84;
t70 = t95 * t81 + t84 * t98;
t88 = g(1) * t58 - g(2) * t60 + g(3) * t70;
t56 = -g(1) * t68 - g(2) * t66 + g(3) * t98;
t71 = -t81 * t98 + t95 * t84;
t65 = t68 * pkin(2);
t64 = t66 * pkin(2);
t61 = -t66 * t81 + t79 * t99;
t59 = t68 * t81 + t77 * t99;
t1 = [(-MDP(1) - t94) * g(3); (-g(1) * (t69 * qJ(3) - t65) - g(2) * (t67 * qJ(3) - t64) - t93) * MDP(7) + (-g(1) * (-t68 * pkin(8) - t65) - g(2) * (-t66 * pkin(8) - t64) - t93 - (pkin(8) * t85 + t90 * t82) * t102 + t105 * (qJ(3) + t90)) * MDP(18) + (-g(1) * (-t68 * t83 - t69 * t97) - g(2) * (-t66 * t83 - t67 * t97) - (-t82 * t97 + t83 * t85) * t102) * MDP(24) + (-g(1) * (t68 * t80 - t69 * t96) - g(2) * (t66 * t80 - t67 * t96) - (-t80 * t85 - t82 * t96) * t102) * MDP(25) + (MDP(5) - MDP(3) - MDP(15)) * t56 + (t106 * t84 - t107 * t81 + MDP(4) - MDP(6)) * (g(3) * t100 - t105); t94 * t56; (-g(1) * (-t58 * pkin(4) + t59 * qJ(5)) - g(2) * (t60 * pkin(4) - t61 * qJ(5)) - g(3) * (-t70 * pkin(4) + t71 * qJ(5))) * MDP(18) + (-MDP(24) * t80 - MDP(25) * t83 - t106) * (g(1) * t59 - g(2) * t61 + g(3) * t71) + t107 * t88; -t88 * MDP(18); (-g(1) * (t58 * t83 - t69 * t80) - g(2) * (-t60 * t83 - t67 * t80) - g(3) * (-t80 * t100 + t70 * t83)) * MDP(24) + (-g(1) * (-t58 * t80 - t69 * t83) - g(2) * (t60 * t80 - t67 * t83) - g(3) * (-t83 * t100 - t70 * t80)) * MDP(25);];
taug  = t1;
