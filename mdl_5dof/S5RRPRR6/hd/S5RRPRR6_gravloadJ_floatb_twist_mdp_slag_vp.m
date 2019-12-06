% Calculate Gravitation load on the joints for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:36:11
% EndTime: 2019-12-05 18:36:11
% DurationCPUTime: 0.16s
% Computational Cost: add. (215->48), mult. (201->76), div. (0->0), fcn. (196->10), ass. (0->31)
t68 = sin(pkin(9));
t82 = g(1) * t68;
t67 = qJ(1) + qJ(2);
t63 = sin(t67);
t69 = cos(pkin(9));
t81 = t63 * t69;
t65 = cos(t67);
t80 = t65 * t69;
t70 = sin(qJ(4));
t79 = t69 * t70;
t72 = cos(qJ(4));
t78 = t69 * t72;
t66 = qJ(4) + qJ(5);
t62 = sin(t66);
t64 = cos(t66);
t46 = t62 * t81 + t65 * t64;
t47 = -t65 * t62 + t64 * t81;
t48 = t62 * t80 - t63 * t64;
t49 = -t63 * t62 - t64 * t80;
t77 = (-g(2) * t46 + g(3) * t48 + t62 * t82) * MDP(23) + (-g(2) * t47 - g(3) * t49 + t64 * t82) * MDP(24);
t76 = -t63 * pkin(2) + t65 * qJ(3);
t60 = g(2) * t65 + g(3) * t63;
t75 = -t65 * pkin(2) - t63 * qJ(3);
t52 = t63 * t79 + t65 * t72;
t53 = t63 * t78 - t65 * t70;
t54 = -t63 * t72 + t65 * t79;
t55 = -t63 * t70 - t65 * t78;
t74 = (-g(2) * t48 - g(3) * t46) * MDP(24) + (-g(2) * t49 + g(3) * t47) * MDP(23) + (-g(2) * t54 - g(3) * t52) * MDP(17) + (-g(2) * t55 + g(3) * t53) * MDP(16) + (MDP(9) - MDP(6)) * (g(2) * t63 - g(3) * t65) + (t69 * MDP(7) - t68 * MDP(8) + MDP(5)) * t60;
t73 = cos(qJ(1));
t71 = sin(qJ(1));
t1 = [(g(2) * t73 + g(3) * t71) * MDP(2) + (-g(2) * t71 + g(3) * t73) * MDP(3) + (-g(2) * (-t73 * pkin(1) + t75) - g(3) * (-t71 * pkin(1) + t76)) * MDP(10) + t74; (-g(2) * t75 - g(3) * t76) * MDP(10) + t74; -t60 * MDP(10); (-g(2) * t52 + g(3) * t54 + t70 * t82) * MDP(16) + (-g(2) * t53 - g(3) * t55 + t72 * t82) * MDP(17) + t77; t77;];
taug = t1;
