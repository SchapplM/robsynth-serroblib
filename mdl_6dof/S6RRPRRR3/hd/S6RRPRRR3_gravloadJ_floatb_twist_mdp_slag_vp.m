% Calculate Gravitation load on the joints for
% S6RRPRRR3
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
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:25:19
% EndTime: 2019-03-09 13:25:20
% DurationCPUTime: 0.27s
% Computational Cost: add. (313->68), mult. (292->98), div. (0->0), fcn. (289->12), ass. (0->50)
t113 = MDP(12) * pkin(2) + MDP(9);
t80 = qJ(4) + qJ(5);
t78 = qJ(6) + t80;
t71 = sin(t78);
t72 = cos(t78);
t76 = sin(t80);
t77 = cos(t80);
t82 = sin(qJ(4));
t85 = cos(qJ(4));
t112 = -MDP(18) * t85 + MDP(19) * t82 - MDP(25) * t77 + MDP(26) * t76 - MDP(32) * t72 + MDP(33) * t71;
t79 = qJ(2) + pkin(11);
t74 = sin(t79);
t111 = g(3) * t74;
t84 = sin(qJ(1));
t110 = t84 * t71;
t109 = t84 * t72;
t108 = t84 * t76;
t107 = t84 * t77;
t106 = t84 * t82;
t105 = t84 * t85;
t87 = cos(qJ(1));
t104 = t87 * t71;
t103 = t87 * t72;
t102 = t87 * t76;
t101 = t87 * t77;
t100 = t87 * t82;
t99 = t87 * t85;
t75 = cos(t79);
t57 = t75 * t110 + t103;
t58 = -t75 * t109 + t104;
t59 = -t75 * t104 + t109;
t60 = t75 * t103 + t110;
t98 = (-g(1) * t59 + g(2) * t57 + t71 * t111) * MDP(32) + (g(1) * t60 - g(2) * t58 + t72 * t111) * MDP(33);
t83 = sin(qJ(2));
t90 = t83 * MDP(10);
t61 = t75 * t108 + t101;
t62 = -t75 * t107 + t102;
t63 = -t75 * t102 + t107;
t64 = t75 * t101 + t108;
t89 = (-g(1) * t63 + g(2) * t61 + t76 * t111) * MDP(25) + (g(1) * t64 - g(2) * t62 + t77 * t111) * MDP(26) + t98;
t88 = g(1) * t87 + g(2) * t84;
t69 = g(1) * t84 - g(2) * t87;
t86 = cos(qJ(2));
t81 = -qJ(3) - pkin(7);
t73 = t86 * pkin(2) + pkin(1);
t68 = t75 * t99 + t106;
t67 = -t75 * t100 + t105;
t66 = -t75 * t105 + t100;
t65 = t75 * t106 + t99;
t1 = [(-g(1) * (-t84 * t73 - t87 * t81) - g(2) * (t87 * t73 - t84 * t81)) * MDP(12) + (-g(1) * t66 - g(2) * t68) * MDP(18) + (-g(1) * t65 - g(2) * t67) * MDP(19) + (-g(1) * t62 - g(2) * t64) * MDP(25) + (-g(1) * t61 - g(2) * t63) * MDP(26) + (-g(1) * t58 - g(2) * t60) * MDP(32) + (-g(1) * t57 - g(2) * t59) * MDP(33) + (MDP(3) - MDP(11)) * t88 + (t86 * MDP(9) + MDP(2) - t90) * t69; (t112 * t75 - t113 * t86 + t90) * g(3) + (MDP(10) * t86 - t112 * t74 + t113 * t83) * t88; -t69 * MDP(12); (-g(1) * t67 + g(2) * t65 + t82 * t111) * MDP(18) + (g(1) * t68 - g(2) * t66 + t85 * t111) * MDP(19) + t89; t89; t98;];
taug  = t1;
