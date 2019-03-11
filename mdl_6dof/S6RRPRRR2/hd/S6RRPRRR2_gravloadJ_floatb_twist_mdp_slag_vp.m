% Calculate Gravitation load on the joints for
% S6RRPRRR2
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
%   see S6RRPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:19:10
% EndTime: 2019-03-09 13:19:10
% DurationCPUTime: 0.21s
% Computational Cost: add. (281->56), mult. (260->80), div. (0->0), fcn. (241->10), ass. (0->36)
t66 = qJ(2) + pkin(11) + qJ(4);
t63 = sin(t66);
t92 = g(3) * t63;
t69 = qJ(5) + qJ(6);
t67 = sin(t69);
t73 = sin(qJ(1));
t90 = t73 * t67;
t68 = cos(t69);
t89 = t73 * t68;
t71 = sin(qJ(5));
t88 = t73 * t71;
t74 = cos(qJ(5));
t87 = t73 * t74;
t76 = cos(qJ(1));
t86 = t76 * t67;
t85 = t76 * t68;
t84 = t76 * t71;
t83 = t76 * t74;
t64 = cos(t66);
t53 = t64 * t90 + t85;
t54 = -t64 * t89 + t86;
t55 = -t64 * t86 + t89;
t56 = t64 * t85 + t90;
t82 = (-g(1) * t55 + g(2) * t53 + t67 * t92) * MDP(32) + (g(1) * t56 - g(2) * t54 + t68 * t92) * MDP(33);
t81 = g(1) * t76 + g(2) * t73;
t61 = g(1) * t73 - g(2) * t76;
t80 = (t81 * t64 + t92) * MDP(19) + (t74 * MDP(25) - t71 * MDP(26) + t68 * MDP(32) - t67 * MDP(33) + MDP(18)) * (-g(3) * t64 + t81 * t63);
t75 = cos(qJ(2));
t72 = sin(qJ(2));
t70 = -qJ(3) - pkin(7);
t65 = t75 * pkin(2) + pkin(1);
t60 = t64 * t83 + t88;
t59 = -t64 * t84 + t87;
t58 = -t64 * t87 + t84;
t57 = t64 * t88 + t83;
t1 = [(-g(1) * (-t73 * t65 - t76 * t70) - g(2) * (t76 * t65 - t73 * t70)) * MDP(12) + (-g(1) * t58 - g(2) * t60) * MDP(25) + (-g(1) * t57 - g(2) * t59) * MDP(26) + (-g(1) * t54 - g(2) * t56) * MDP(32) + (-g(1) * t53 - g(2) * t55) * MDP(33) + (MDP(3) - MDP(11)) * t81 + (-MDP(10) * t72 + MDP(18) * t64 - MDP(19) * t63 + MDP(9) * t75 + MDP(2)) * t61; (g(3) * t72 + t81 * t75) * MDP(10) + t80 + (MDP(12) * pkin(2) + MDP(9)) * (-g(3) * t75 + t81 * t72); -t61 * MDP(12); t80; (-g(1) * t59 + g(2) * t57 + t71 * t92) * MDP(25) + (g(1) * t60 - g(2) * t58 + t74 * t92) * MDP(26) + t82; t82;];
taug  = t1;
