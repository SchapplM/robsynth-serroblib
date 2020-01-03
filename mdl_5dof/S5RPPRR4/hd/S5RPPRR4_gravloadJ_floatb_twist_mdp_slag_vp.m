% Calculate Gravitation load on the joints for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:31:23
% EndTime: 2020-01-03 11:31:24
% DurationCPUTime: 0.21s
% Computational Cost: add. (154->55), mult. (178->85), div. (0->0), fcn. (175->10), ass. (0->40)
t68 = sin(pkin(8));
t90 = g(1) * t68;
t66 = pkin(9) + qJ(4);
t62 = qJ(5) + t66;
t58 = sin(t62);
t71 = sin(qJ(1));
t89 = t71 * t58;
t59 = cos(t62);
t88 = t71 * t59;
t60 = sin(t66);
t87 = t71 * t60;
t61 = cos(t66);
t86 = t71 * t61;
t67 = sin(pkin(9));
t85 = t71 * t67;
t69 = cos(pkin(9));
t84 = t71 * t69;
t72 = cos(qJ(1));
t83 = t72 * t58;
t82 = t72 * t59;
t81 = t72 * t60;
t80 = t72 * t61;
t79 = t72 * t67;
t78 = t72 * t69;
t70 = cos(pkin(8));
t47 = -t70 * t89 - t82;
t48 = t70 * t88 - t83;
t49 = t70 * t83 - t88;
t50 = t70 * t82 + t89;
t77 = (-g(2) * t47 - g(3) * t49 + t58 * t90) * MDP(24) + (g(2) * t48 - g(3) * t50 + t59 * t90) * MDP(25);
t76 = t72 * pkin(1) + t71 * qJ(2);
t75 = t71 * pkin(1) - t72 * qJ(2);
t57 = g(2) * t72 + g(3) * t71;
t74 = g(2) * t71 - g(3) * t72;
t73 = pkin(2) * t70 + qJ(3) * t68;
t54 = t70 * t80 + t87;
t53 = t70 * t81 - t86;
t52 = t70 * t86 - t81;
t51 = -t70 * t87 - t80;
t1 = [(-g(2) * t76 - g(3) * t75) * MDP(7) + (-g(2) * (t70 * t78 + t85) - g(3) * (t70 * t84 - t79)) * MDP(8) + (-g(2) * (-t70 * t79 + t84) - g(3) * (-t70 * t85 - t78)) * MDP(9) + (-g(2) * (t72 * t73 + t76) - g(3) * (t71 * t73 + t75)) * MDP(11) + (-g(2) * t54 - g(3) * t52) * MDP(17) + (g(2) * t53 - g(3) * t51) * MDP(18) + (-g(2) * t50 - g(3) * t48) * MDP(24) + (g(2) * t49 - g(3) * t47) * MDP(25) + (MDP(3) - MDP(6)) * t74 + (-MDP(2) - t70 * MDP(4) + (MDP(5) - MDP(10)) * t68) * t57; (MDP(11) + MDP(7)) * t57; (g(1) * t70 - t68 * t74) * MDP(11); (-g(2) * t51 - g(3) * t53 + t60 * t90) * MDP(17) + (g(2) * t52 - g(3) * t54 + t61 * t90) * MDP(18) + t77; t77;];
taug = t1;
