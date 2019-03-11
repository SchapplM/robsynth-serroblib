% Calculate Gravitation load on the joints for
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:05:30
% EndTime: 2019-03-09 14:05:32
% DurationCPUTime: 0.32s
% Computational Cost: add. (370->79), mult. (350->118), div. (0->0), fcn. (346->12), ass. (0->45)
t85 = sin(qJ(1));
t87 = cos(qJ(1));
t93 = g(1) * t87 + g(2) * t85;
t109 = MDP(10) - MDP(13);
t84 = sin(qJ(2));
t86 = cos(qJ(2));
t70 = -g(3) * t86 + t84 * t93;
t106 = g(3) * t84;
t104 = t85 * t86;
t81 = pkin(11) + qJ(4);
t80 = qJ(5) + t81;
t77 = qJ(6) + t80;
t73 = sin(t77);
t103 = t87 * t73;
t74 = cos(t77);
t102 = t87 * t74;
t75 = sin(t80);
t101 = t87 * t75;
t76 = cos(t80);
t100 = t87 * t76;
t78 = sin(t81);
t99 = t87 * t78;
t79 = cos(t81);
t98 = t87 * t79;
t82 = sin(pkin(11));
t97 = t87 * t82;
t83 = cos(pkin(11));
t96 = t87 * t83;
t58 = t104 * t73 + t102;
t59 = -t104 * t74 + t103;
t60 = -t103 * t86 + t85 * t74;
t61 = t102 * t86 + t85 * t73;
t95 = (-g(1) * t60 + g(2) * t58 + t106 * t73) * MDP(34) + (g(1) * t61 - g(2) * t59 + t106 * t74) * MDP(35);
t62 = t104 * t75 + t100;
t63 = -t104 * t76 + t101;
t64 = -t101 * t86 + t85 * t76;
t65 = t100 * t86 + t85 * t75;
t94 = (-g(1) * t64 + g(2) * t62 + t106 * t75) * MDP(27) + (g(1) * t65 - g(2) * t63 + t106 * t76) * MDP(28) + t95;
t91 = t86 * pkin(2) + t84 * qJ(3);
t89 = pkin(1) + t91;
t69 = t85 * t78 + t86 * t98;
t68 = t85 * t79 - t86 * t99;
t67 = -t104 * t79 + t99;
t66 = t104 * t78 + t98;
t1 = [t93 * MDP(3) + (-g(1) * (-t104 * t83 + t97) - g(2) * (t85 * t82 + t86 * t96)) * MDP(11) + (-g(1) * (t104 * t82 + t96) - g(2) * (t85 * t83 - t86 * t97)) * MDP(12) + ((-g(1) * pkin(7) - g(2) * t89) * t87 + (-g(2) * pkin(7) + g(1) * t89) * t85) * MDP(14) + (-g(1) * t67 - g(2) * t69) * MDP(20) + (-g(1) * t66 - g(2) * t68) * MDP(21) + (-g(1) * t63 - g(2) * t65) * MDP(27) + (-g(1) * t62 - g(2) * t64) * MDP(28) + (-g(1) * t59 - g(2) * t61) * MDP(34) + (-g(1) * t58 - g(2) * t60) * MDP(35) + (MDP(9) * t86 - t109 * t84 + MDP(2)) * (g(1) * t85 - g(2) * t87); (-g(3) * t91 + t93 * (pkin(2) * t84 - qJ(3) * t86)) * MDP(14) + t109 * (t86 * t93 + t106) + (MDP(11) * t83 - MDP(12) * t82 + MDP(20) * t79 - MDP(21) * t78 + MDP(27) * t76 - MDP(28) * t75 + MDP(34) * t74 - MDP(35) * t73 + MDP(9)) * t70; -t70 * MDP(14); (-g(1) * t68 + g(2) * t66 + t106 * t78) * MDP(20) + (g(1) * t69 - g(2) * t67 + t106 * t79) * MDP(21) + t94; t94; t95;];
taug  = t1;
