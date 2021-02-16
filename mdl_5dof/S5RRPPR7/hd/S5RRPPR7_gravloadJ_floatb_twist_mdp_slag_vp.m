% Calculate Gravitation load on the joints for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:59:02
% EndTime: 2021-01-15 19:59:03
% DurationCPUTime: 0.26s
% Computational Cost: add. (145->58), mult. (210->78), div. (0->0), fcn. (181->10), ass. (0->33)
t84 = MDP(11) - MDP(16);
t83 = MDP(12) - MDP(17);
t63 = qJ(2) + pkin(8);
t60 = cos(t63);
t56 = g(3) * t60;
t66 = -qJ(3) - pkin(6);
t69 = sin(qJ(1));
t80 = t69 * t66;
t67 = sin(qJ(5));
t79 = t69 * t67;
t70 = cos(qJ(5));
t78 = t69 * t70;
t72 = cos(qJ(1));
t77 = t72 * t67;
t76 = t72 * t70;
t75 = g(1) * t72 + g(2) * t69;
t54 = g(1) * t69 - g(2) * t72;
t71 = cos(qJ(2));
t68 = sin(qJ(2));
t65 = cos(pkin(8));
t64 = sin(pkin(8));
t61 = t71 * pkin(2);
t59 = sin(t63);
t58 = t61 + pkin(1);
t57 = t66 * t72;
t53 = -t64 * pkin(3) + qJ(4) * t65;
t52 = pkin(3) * t65 + qJ(4) * t64 + pkin(2);
t50 = -t59 * t79 + t76;
t49 = t59 * t78 + t77;
t48 = t59 * t77 + t78;
t47 = t59 * t76 - t79;
t41 = t52 * t71 + t53 * t68 + pkin(1);
t1 = [(-g(1) * (-t69 * t58 - t57) - g(2) * (t72 * t58 - t80)) * MDP(14) + (-g(1) * (-t41 * t69 - t57) - g(2) * (t41 * t72 - t80)) * MDP(18) + (-g(1) * t50 - g(2) * t48) * MDP(24) + (g(1) * t49 - g(2) * t47) * MDP(25) + (MDP(3) - MDP(13) - MDP(15)) * t75 + (-t68 * MDP(10) + t71 * MDP(9) - t83 * t59 + t84 * t60 + MDP(2)) * t54; (g(3) * t68 + t75 * t71) * MDP(10) + (-g(3) * (t60 * pkin(3) + t59 * qJ(4) + t61) - t75 * (-t52 * t68 + t53 * t71)) * MDP(18) + (MDP(14) * pkin(2) + MDP(9)) * (-g(3) * t71 + t75 * t68) + t84 * (t75 * t59 - t56) + (-MDP(24) * t67 - MDP(25) * t70 + t83) * (g(3) * t59 + t75 * t60); (-MDP(14) - MDP(18)) * t54; (t56 - t75 * (t64 * t71 + t65 * t68)) * MDP(18); (-g(1) * t47 - g(2) * t49 + t70 * t56) * MDP(24) + (g(1) * t48 - g(2) * t50 - t67 * t56) * MDP(25);];
taug = t1;
