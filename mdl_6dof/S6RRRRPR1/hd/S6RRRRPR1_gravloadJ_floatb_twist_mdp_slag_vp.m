% Calculate Gravitation load on the joints for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:54:55
% EndTime: 2019-03-09 21:54:56
% DurationCPUTime: 0.20s
% Computational Cost: add. (294->56), mult. (257->80), div. (0->0), fcn. (218->12), ass. (0->35)
t72 = qJ(2) + qJ(3);
t70 = qJ(4) + t72;
t64 = pkin(11) + t70;
t60 = sin(t64);
t91 = g(3) * t60;
t73 = sin(qJ(6));
t78 = cos(qJ(1));
t89 = t73 * t78;
t75 = sin(qJ(1));
t88 = t75 * t73;
t76 = cos(qJ(6));
t87 = t75 * t76;
t86 = t76 * t78;
t66 = cos(t70);
t68 = cos(t72);
t85 = pkin(3) * t68 + pkin(4) * t66;
t77 = cos(qJ(2));
t84 = t77 * pkin(2) + t85;
t61 = cos(t64);
t65 = sin(t70);
t81 = g(1) * t78 + g(2) * t75;
t79 = -g(3) * t66 + t65 * t81;
t83 = t79 * MDP(23) + (g(3) * t65 + t66 * t81) * MDP(24) + (t76 * MDP(32) - t73 * MDP(33)) * (-g(3) * t61 + t60 * t81);
t67 = sin(t72);
t82 = -pkin(3) * t67 - pkin(4) * t65;
t58 = g(1) * t75 - g(2) * t78;
t80 = (-g(3) * t68 + t67 * t81) * MDP(16) + (g(3) * t67 + t68 * t81) * MDP(17) + t83;
t74 = sin(qJ(2));
t69 = -qJ(5) - pkin(9) - pkin(8) - pkin(7);
t55 = pkin(1) + t84;
t54 = t61 * t86 + t88;
t53 = -t61 * t89 + t87;
t52 = -t61 * t87 + t89;
t51 = t61 * t88 + t86;
t1 = [(-g(1) * (-t75 * t55 - t69 * t78) - g(2) * (t55 * t78 - t75 * t69)) * MDP(26) + (-g(1) * t52 - g(2) * t54) * MDP(32) + (-g(1) * t51 - g(2) * t53) * MDP(33) + (MDP(3) - MDP(25)) * t81 + (-t74 * MDP(10) + MDP(16) * t68 - MDP(17) * t67 + MDP(23) * t66 - MDP(24) * t65 + t77 * MDP(9) + MDP(2)) * t58; (-g(3) * t77 + t74 * t81) * MDP(9) + (g(3) * t74 + t77 * t81) * MDP(10) + (-g(3) * t84 - t81 * (-pkin(2) * t74 + t82)) * MDP(26) + t80; (-g(3) * t85 - t81 * t82) * MDP(26) + t80; MDP(26) * pkin(4) * t79 + t83; -t58 * MDP(26); (-g(1) * t53 + g(2) * t51 + t73 * t91) * MDP(32) + (g(1) * t54 - g(2) * t52 + t76 * t91) * MDP(33);];
taug  = t1;
