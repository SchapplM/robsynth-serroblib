% Calculate Gravitation load on the joints for
% S6RRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:35:48
% EndTime: 2019-03-10 03:35:49
% DurationCPUTime: 0.22s
% Computational Cost: add. (358->54), mult. (318->78), div. (0->0), fcn. (296->12), ass. (0->37)
t73 = qJ(2) + qJ(3);
t71 = qJ(4) + t73;
t65 = sin(t71);
t96 = g(3) * t65;
t72 = qJ(5) + qJ(6);
t67 = sin(t72);
t76 = sin(qJ(1));
t94 = t76 * t67;
t69 = cos(t72);
t93 = t76 * t69;
t74 = sin(qJ(5));
t92 = t76 * t74;
t77 = cos(qJ(5));
t91 = t76 * t77;
t79 = cos(qJ(1));
t90 = t79 * t67;
t89 = t79 * t69;
t88 = t79 * t74;
t87 = t79 * t77;
t66 = cos(t71);
t57 = t66 * t94 + t89;
t58 = -t66 * t93 + t90;
t59 = -t66 * t90 + t93;
t60 = t66 * t89 + t94;
t86 = (-g(1) * t59 + g(2) * t57 + t67 * t96) * MDP(37) + (g(1) * t60 - g(2) * t58 + t69 * t96) * MDP(38);
t85 = g(1) * t79 + g(2) * t76;
t83 = (t85 * t66 + t96) * MDP(24) + (t77 * MDP(30) - t74 * MDP(31) + t69 * MDP(37) - t67 * MDP(38) + MDP(23)) * (-g(3) * t66 + t85 * t65);
t68 = sin(t73);
t70 = cos(t73);
t82 = (-g(3) * t70 + t68 * t85) * MDP(16) + (g(3) * t68 + t70 * t85) * MDP(17) + t83;
t78 = cos(qJ(2));
t75 = sin(qJ(2));
t64 = t66 * t87 + t92;
t63 = -t66 * t88 + t91;
t62 = -t66 * t91 + t88;
t61 = t66 * t92 + t87;
t1 = [t85 * MDP(3) + (-g(1) * t62 - g(2) * t64) * MDP(30) + (-g(1) * t61 - g(2) * t63) * MDP(31) + (-g(1) * t58 - g(2) * t60) * MDP(37) + (-g(1) * t57 - g(2) * t59) * MDP(38) + (-t75 * MDP(10) + MDP(16) * t70 - MDP(17) * t68 + MDP(23) * t66 - MDP(24) * t65 + t78 * MDP(9) + MDP(2)) * (g(1) * t76 - g(2) * t79); (-g(3) * t78 + t75 * t85) * MDP(9) + (g(3) * t75 + t78 * t85) * MDP(10) + t82; t82; t83; (-g(1) * t63 + g(2) * t61 + t74 * t96) * MDP(30) + (g(1) * t64 - g(2) * t62 + t77 * t96) * MDP(31) + t86; t86;];
taug  = t1;
