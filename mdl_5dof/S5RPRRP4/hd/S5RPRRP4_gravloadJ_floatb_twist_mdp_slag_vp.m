% Calculate Gravitation load on the joints for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:06:53
% EndTime: 2019-12-05 18:06:55
% DurationCPUTime: 0.26s
% Computational Cost: add. (141->60), mult. (207->89), div. (0->0), fcn. (194->8), ass. (0->36)
t68 = cos(pkin(8));
t66 = qJ(3) + qJ(4);
t62 = cos(t66);
t72 = cos(qJ(1));
t78 = t72 * t62;
t61 = sin(t66);
t70 = sin(qJ(1));
t83 = t70 * t61;
t46 = t68 * t83 + t78;
t79 = t72 * t61;
t82 = t70 * t62;
t48 = t68 * t79 - t82;
t67 = sin(pkin(8));
t87 = g(1) * t67;
t88 = -g(2) * t46 + g(3) * t48 + t61 * t87;
t69 = sin(qJ(3));
t56 = pkin(3) * t69 + pkin(4) * t61;
t84 = t56 * t68;
t81 = t70 * t69;
t71 = cos(qJ(3));
t80 = t70 * t71;
t77 = t72 * t69;
t76 = t72 * t71;
t47 = t68 * t82 - t79;
t49 = -t68 * t78 - t83;
t75 = t88 * MDP(20) + (-g(2) * t47 - g(3) * t49 + t62 * t87) * MDP(21);
t57 = t71 * pkin(3) + pkin(4) * t62;
t59 = g(2) * t72 + g(3) * t70;
t58 = g(2) * t70 - g(3) * t72;
t73 = (pkin(2) + t57) * t68 - (-qJ(5) - pkin(7) - pkin(6)) * t67 + pkin(1);
t63 = t72 * qJ(2);
t53 = -t68 * t76 - t81;
t52 = t68 * t77 - t80;
t51 = t68 * t80 - t77;
t50 = t68 * t81 + t76;
t1 = [(-g(2) * (-t72 * pkin(1) - t70 * qJ(2)) - g(3) * (-pkin(1) * t70 + t63)) * MDP(7) + (-g(2) * t53 + g(3) * t51) * MDP(13) + (-g(2) * t52 - g(3) * t50) * MDP(14) + (-g(2) * t49 + g(3) * t47) * MDP(20) + (-g(2) * t48 - g(3) * t46) * MDP(21) + (-g(3) * t63 + (g(2) * t73 - g(3) * t56) * t72 + (-g(2) * (-qJ(2) - t56) + g(3) * t73) * t70) * MDP(23) + (-MDP(3) + MDP(6)) * t58 + (t68 * MDP(4) + MDP(2) + (-MDP(5) + MDP(22)) * t67) * t59; (-MDP(23) - MDP(7)) * t59; (-g(2) * t50 + g(3) * t52 + t69 * t87) * MDP(13) + (-g(2) * t51 - g(3) * t53 + t71 * t87) * MDP(14) + (t56 * t87 - g(2) * (t57 * t72 + t70 * t84) - g(3) * (t57 * t70 - t72 * t84)) * MDP(23) + t75; MDP(23) * pkin(4) * t88 + t75; (g(1) * t68 + t58 * t67) * MDP(23);];
taug = t1;
