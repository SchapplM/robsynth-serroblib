% Calculate Gravitation load on the joints for
% S6RRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:08:40
% EndTime: 2019-03-09 18:08:41
% DurationCPUTime: 0.26s
% Computational Cost: add. (283->62), mult. (274->91), div. (0->0), fcn. (251->12), ass. (0->41)
t83 = qJ(2) + qJ(3);
t75 = pkin(11) + t83;
t72 = sin(t75);
t106 = g(3) * t72;
t82 = qJ(5) + qJ(6);
t76 = sin(t82);
t86 = sin(qJ(1));
t104 = t86 * t76;
t78 = cos(t82);
t103 = t86 * t78;
t84 = sin(qJ(5));
t102 = t86 * t84;
t87 = cos(qJ(5));
t101 = t86 * t87;
t89 = cos(qJ(1));
t100 = t89 * t76;
t99 = t89 * t78;
t98 = t89 * t84;
t97 = t89 * t87;
t73 = cos(t75);
t60 = t73 * t104 + t99;
t61 = -t73 * t103 + t100;
t62 = -t73 * t100 + t103;
t63 = t73 * t99 + t104;
t96 = (-g(1) * t62 + g(2) * t60 + t76 * t106) * MDP(32) + (g(1) * t63 - g(2) * t61 + t78 * t106) * MDP(33);
t79 = cos(t83);
t88 = cos(qJ(2));
t95 = t88 * pkin(2) + pkin(3) * t79;
t94 = g(1) * t89 + g(2) * t86;
t70 = g(1) * t86 - g(2) * t89;
t77 = sin(t83);
t90 = -g(3) * t79 + t94 * t77;
t93 = t90 * MDP(16) + (g(3) * t77 + t94 * t79) * MDP(17) + (t87 * MDP(25) - t84 * MDP(26) + t78 * MDP(32) - t76 * MDP(33)) * (-g(3) * t73 + t94 * t72);
t85 = sin(qJ(2));
t81 = -qJ(4) - pkin(8) - pkin(7);
t68 = pkin(1) + t95;
t67 = t73 * t97 + t102;
t66 = -t73 * t98 + t101;
t65 = -t73 * t101 + t98;
t64 = t73 * t102 + t97;
t1 = [(-g(1) * (-t86 * t68 - t89 * t81) - g(2) * (t89 * t68 - t86 * t81)) * MDP(19) + (-g(1) * t65 - g(2) * t67) * MDP(25) + (-g(1) * t64 - g(2) * t66) * MDP(26) + (-g(1) * t61 - g(2) * t63) * MDP(32) + (-g(1) * t60 - g(2) * t62) * MDP(33) + (MDP(3) - MDP(18)) * t94 + (-t85 * MDP(10) + MDP(16) * t79 - MDP(17) * t77 + t88 * MDP(9) + MDP(2)) * t70; (-g(3) * t88 + t94 * t85) * MDP(9) + (g(3) * t85 + t94 * t88) * MDP(10) + (-g(3) * t95 - t94 * (-t85 * pkin(2) - pkin(3) * t77)) * MDP(19) + t93; t90 * MDP(19) * pkin(3) + t93; -t70 * MDP(19); (-g(1) * t66 + g(2) * t64 + t84 * t106) * MDP(25) + (g(1) * t67 - g(2) * t65 + t87 * t106) * MDP(26) + t96; t96;];
taug  = t1;
