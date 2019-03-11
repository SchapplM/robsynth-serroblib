% Calculate Gravitation load on the joints for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:58:36
% EndTime: 2019-03-09 06:58:36
% DurationCPUTime: 0.16s
% Computational Cost: add. (280->52), mult. (242->81), div. (0->0), fcn. (228->12), ass. (0->32)
t66 = qJ(3) + qJ(4);
t61 = sin(t66);
t85 = g(3) * t61;
t65 = qJ(5) + qJ(6);
t60 = sin(t65);
t63 = cos(t66);
t83 = t60 * t63;
t62 = cos(t65);
t82 = t62 * t63;
t67 = sin(qJ(5));
t81 = t63 * t67;
t70 = cos(qJ(5));
t80 = t63 * t70;
t64 = qJ(1) + pkin(11);
t58 = sin(t64);
t59 = cos(t64);
t50 = t58 * t83 + t59 * t62;
t51 = -t58 * t82 + t59 * t60;
t52 = t58 * t62 - t59 * t83;
t53 = t58 * t60 + t59 * t82;
t79 = (-g(1) * t52 + g(2) * t50 + t60 * t85) * MDP(31) + (g(1) * t53 - g(2) * t51 + t62 * t85) * MDP(32);
t78 = g(1) * t59 + g(2) * t58;
t75 = (t78 * t63 + t85) * MDP(18) + (t70 * MDP(24) - t67 * MDP(25) + t62 * MDP(31) - t60 * MDP(32) + MDP(17)) * (-g(3) * t63 + t78 * t61);
t72 = cos(qJ(1));
t71 = cos(qJ(3));
t69 = sin(qJ(1));
t68 = sin(qJ(3));
t57 = t58 * t67 + t59 * t80;
t56 = t58 * t70 - t59 * t81;
t55 = -t58 * t80 + t59 * t67;
t54 = t58 * t81 + t59 * t70;
t1 = [(g(1) * t72 + g(2) * t69) * MDP(3) + (-g(1) * t55 - g(2) * t57) * MDP(24) + (-g(1) * t54 - g(2) * t56) * MDP(25) + (-g(1) * t51 - g(2) * t53) * MDP(31) + (-g(1) * t50 - g(2) * t52) * MDP(32) + (MDP(4) * pkin(1) + MDP(2)) * (g(1) * t69 - g(2) * t72) + (t71 * MDP(10) - t68 * MDP(11) + MDP(17) * t63 - MDP(18) * t61) * (g(1) * t58 - g(2) * t59); -g(3) * MDP(4); (-g(3) * t71 + t78 * t68) * MDP(10) + (g(3) * t68 + t78 * t71) * MDP(11) + t75; t75; (-g(1) * t56 + g(2) * t54 + t67 * t85) * MDP(24) + (g(1) * t57 - g(2) * t55 + t70 * t85) * MDP(25) + t79; t79;];
taug  = t1;
