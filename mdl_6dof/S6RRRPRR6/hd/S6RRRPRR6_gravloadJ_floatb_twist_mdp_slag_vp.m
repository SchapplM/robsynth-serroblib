% Calculate Gravitation load on the joints for
% S6RRRPRR6
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
%   see S6RRRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:29:31
% EndTime: 2019-03-09 18:29:32
% DurationCPUTime: 0.26s
% Computational Cost: add. (349->72), mult. (339->105), div. (0->0), fcn. (331->10), ass. (0->43)
t85 = sin(qJ(1));
t88 = cos(qJ(1));
t94 = g(1) * t88 + g(2) * t85;
t113 = MDP(10) - MDP(18);
t84 = sin(qJ(2));
t87 = cos(qJ(2));
t68 = -g(3) * t87 + t94 * t84;
t107 = g(3) * t84;
t105 = t85 * t87;
t81 = qJ(3) + pkin(11) + qJ(5);
t79 = qJ(6) + t81;
t75 = sin(t79);
t104 = t88 * t75;
t76 = cos(t79);
t103 = t88 * t76;
t77 = sin(t81);
t102 = t88 * t77;
t78 = cos(t81);
t101 = t88 * t78;
t83 = sin(qJ(3));
t100 = t88 * t83;
t86 = cos(qJ(3));
t99 = t88 * t86;
t60 = t75 * t105 + t103;
t61 = -t76 * t105 + t104;
t62 = -t87 * t104 + t85 * t76;
t63 = t87 * t103 + t85 * t75;
t98 = (-g(1) * t62 + g(2) * t60 + t75 * t107) * MDP(32) + (g(1) * t63 - g(2) * t61 + t76 * t107) * MDP(33);
t96 = pkin(3) * t83 + pkin(7);
t64 = t77 * t105 + t101;
t65 = -t78 * t105 + t102;
t66 = -t87 * t102 + t85 * t78;
t67 = t87 * t101 + t85 * t77;
t95 = (-g(1) * t66 + g(2) * t64 + t77 * t107) * MDP(25) + (g(1) * t67 - g(2) * t65 + t78 * t107) * MDP(26) + t98;
t80 = t86 * pkin(3) + pkin(2);
t82 = -qJ(4) - pkin(8);
t92 = t87 * t80 - t84 * t82;
t90 = pkin(1) + t92;
t72 = -t87 * t100 + t85 * t86;
t70 = t83 * t105 + t99;
t73 = t85 * t83 + t87 * t99;
t71 = -t86 * t105 + t100;
t1 = [t94 * MDP(3) + (-g(1) * t71 - g(2) * t73) * MDP(16) + (-g(1) * t70 - g(2) * t72) * MDP(17) + ((-g(1) * t96 - g(2) * t90) * t88 + (g(1) * t90 - g(2) * t96) * t85) * MDP(19) + (-g(1) * t65 - g(2) * t67) * MDP(25) + (-g(1) * t64 - g(2) * t66) * MDP(26) + (-g(1) * t61 - g(2) * t63) * MDP(32) + (-g(1) * t60 - g(2) * t62) * MDP(33) + (t87 * MDP(9) - t113 * t84 + MDP(2)) * (g(1) * t85 - g(2) * t88); (-g(3) * t92 + t94 * (t80 * t84 + t82 * t87)) * MDP(19) + t113 * (t94 * t87 + t107) + (t86 * MDP(16) - t83 * MDP(17) + t78 * MDP(25) - t77 * MDP(26) + t76 * MDP(32) - t75 * MDP(33) + MDP(9)) * t68; (g(1) * t73 - g(2) * t71 + t86 * t107) * MDP(17) + t95 + (pkin(3) * MDP(19) + MDP(16)) * (-g(1) * t72 + g(2) * t70 + t83 * t107); -t68 * MDP(19); t95; t98;];
taug  = t1;
