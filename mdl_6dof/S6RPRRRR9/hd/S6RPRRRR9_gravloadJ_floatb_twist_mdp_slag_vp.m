% Calculate Gravitation load on the joints for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% MDP [34x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(34,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [34 1]), ...
  'S6RPRRRR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [34x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:25:19
% EndTime: 2019-03-09 07:25:20
% DurationCPUTime: 0.21s
% Computational Cost: add. (242->64), mult. (284->94), div. (0->0), fcn. (284->10), ass. (0->43)
t77 = cos(qJ(3));
t102 = t77 * MDP(13);
t72 = qJ(4) + qJ(5);
t70 = qJ(6) + t72;
t66 = sin(t70);
t67 = cos(t70);
t68 = sin(t72);
t69 = cos(t72);
t73 = sin(qJ(4));
t76 = cos(qJ(4));
t101 = t76 * MDP(19) - t73 * MDP(20) + t69 * MDP(26) - t68 * MDP(27) + t67 * MDP(33) - t66 * MDP(34) + MDP(12);
t100 = g(3) * t77;
t75 = sin(qJ(1));
t99 = t75 * t66;
t98 = t75 * t67;
t97 = t75 * t68;
t96 = t75 * t69;
t95 = t75 * t73;
t94 = t75 * t76;
t78 = cos(qJ(1));
t93 = t78 * t66;
t92 = t78 * t67;
t91 = t78 * t68;
t90 = t78 * t69;
t89 = t78 * t73;
t88 = t78 * t76;
t74 = sin(qJ(3));
t52 = -t74 * t99 + t92;
t53 = t74 * t98 + t93;
t54 = t74 * t93 + t98;
t55 = t74 * t92 - t99;
t87 = (-g(1) * t52 - g(2) * t54 + t66 * t100) * MDP(33) + (g(1) * t53 - g(2) * t55 + t67 * t100) * MDP(34);
t56 = -t74 * t97 + t90;
t57 = t74 * t96 + t91;
t58 = t74 * t91 + t96;
t59 = t74 * t90 - t97;
t80 = (-g(1) * t56 - g(2) * t58 + t68 * t100) * MDP(26) + (g(1) * t57 - g(2) * t59 + t69 * t100) * MDP(27) + t87;
t64 = g(1) * t75 - g(2) * t78;
t63 = t74 * t88 - t95;
t62 = t74 * t89 + t94;
t61 = t74 * t94 + t89;
t60 = -t74 * t95 + t88;
t1 = [(-g(1) * (-t75 * pkin(1) + t78 * qJ(2)) - g(2) * (t78 * pkin(1) + t75 * qJ(2))) * MDP(6) + (-g(1) * t63 - g(2) * t61) * MDP(19) + (g(1) * t62 - g(2) * t60) * MDP(20) + (-g(1) * t59 - g(2) * t57) * MDP(26) + (g(1) * t58 - g(2) * t56) * MDP(27) + (-g(1) * t55 - g(2) * t53) * MDP(33) + (g(1) * t54 - g(2) * t52) * MDP(34) + (MDP(2) - MDP(4)) * t64 + (-t74 * MDP(12) + MDP(3) - MDP(5) - t102) * (g(1) * t78 + g(2) * t75); -t64 * MDP(6); (t101 * t74 + t102) * g(3) + (MDP(13) * t74 - t101 * t77) * t64; (-g(1) * t60 - g(2) * t62 + t73 * t100) * MDP(19) + (g(1) * t61 - g(2) * t63 + t76 * t100) * MDP(20) + t80; t80; t87;];
taug  = t1;
