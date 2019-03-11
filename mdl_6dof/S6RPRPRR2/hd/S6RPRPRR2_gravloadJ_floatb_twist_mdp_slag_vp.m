% Calculate Gravitation load on the joints for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:38:58
% EndTime: 2019-03-09 03:38:58
% DurationCPUTime: 0.19s
% Computational Cost: add. (225->58), mult. (199->84), div. (0->0), fcn. (185->12), ass. (0->41)
t91 = MDP(13) * pkin(3) + MDP(10);
t63 = qJ(5) + qJ(6);
t59 = sin(t63);
t60 = cos(t63);
t65 = sin(qJ(5));
t68 = cos(qJ(5));
t90 = -MDP(19) * t68 + MDP(20) * t65 - MDP(26) * t60 + MDP(27) * t59;
t61 = qJ(3) + pkin(11);
t55 = sin(t61);
t89 = g(3) * t55;
t62 = qJ(1) + pkin(10);
t56 = sin(t62);
t88 = t56 * t59;
t87 = t56 * t60;
t86 = t56 * t65;
t85 = t56 * t68;
t58 = cos(t62);
t84 = t58 * t59;
t83 = t58 * t60;
t82 = t58 * t65;
t81 = t58 * t68;
t57 = cos(t61);
t46 = t57 * t88 + t83;
t47 = -t57 * t87 + t84;
t48 = -t57 * t84 + t87;
t49 = t57 * t83 + t88;
t80 = (-g(1) * t48 + g(2) * t46 + t59 * t89) * MDP(26) + (g(1) * t49 - g(2) * t47 + t60 * t89) * MDP(27);
t66 = sin(qJ(3));
t74 = t66 * MDP(11);
t73 = g(1) * t58 + g(2) * t56;
t72 = g(1) * t56 - g(2) * t58;
t70 = cos(qJ(1));
t69 = cos(qJ(3));
t67 = sin(qJ(1));
t64 = -qJ(4) - pkin(7);
t54 = pkin(3) * t69 + pkin(2);
t53 = t57 * t81 + t86;
t52 = -t57 * t82 + t85;
t51 = -t57 * t85 + t82;
t50 = t57 * t86 + t81;
t1 = [(g(1) * t70 + g(2) * t67) * MDP(3) - t73 * MDP(12) + (-g(1) * (-pkin(1) * t67 - t54 * t56 - t58 * t64) - g(2) * (pkin(1) * t70 + t58 * t54 - t56 * t64)) * MDP(13) + (-g(1) * t51 - g(2) * t53) * MDP(19) + (-g(1) * t50 - g(2) * t52) * MDP(20) + (-g(1) * t47 - g(2) * t49) * MDP(26) + (-g(1) * t46 - g(2) * t48) * MDP(27) + (t69 * MDP(10) - t74) * t72 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t67 - g(2) * t70); (-MDP(13) - MDP(4)) * g(3); (t90 * t57 - t91 * t69 + t74) * g(3) + (MDP(11) * t69 - t90 * t55 + t91 * t66) * t73; -t72 * MDP(13); (-g(1) * t52 + g(2) * t50 + t65 * t89) * MDP(19) + (g(1) * t53 - g(2) * t51 + t68 * t89) * MDP(20) + t80; t80;];
taug  = t1;
