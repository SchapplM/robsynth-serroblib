% Calculate Gravitation load on the joints for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:24:49
% EndTime: 2021-01-15 22:24:51
% DurationCPUTime: 0.29s
% Computational Cost: add. (260->61), mult. (240->75), div. (0->0), fcn. (189->8), ass. (0->34)
t82 = MDP(18) + MDP(22);
t81 = MDP(19) - MDP(24);
t68 = sin(qJ(1));
t70 = cos(qJ(1));
t54 = g(1) * t70 + g(2) * t68;
t66 = qJ(2) + qJ(3);
t61 = sin(t66);
t80 = t54 * t61;
t60 = pkin(8) + t66;
t57 = sin(t60);
t58 = cos(t60);
t72 = t58 * pkin(4) + t57 * qJ(5);
t79 = MDP(21) + MDP(25);
t78 = pkin(4) * t57;
t62 = cos(t66);
t77 = g(3) * t62;
t59 = pkin(3) * t62;
t69 = cos(qJ(2));
t63 = t69 * pkin(2);
t76 = t59 + t63;
t75 = qJ(5) * t58;
t74 = t59 + t72;
t67 = sin(qJ(2));
t50 = -t67 * pkin(2) - pkin(3) * t61;
t73 = t50 - t78;
t53 = g(1) * t68 - g(2) * t70;
t42 = -g(3) * t58 + t54 * t57;
t71 = (-t77 + t80) * MDP(16) + (g(3) * t61 + t54 * t62) * MDP(17) + t81 * (g(3) * t57 + t54 * t58) + t82 * t42;
t65 = -qJ(4) - pkin(7) - pkin(6);
t52 = t70 * t75;
t51 = t68 * t75;
t49 = pkin(1) + t76;
t48 = t70 * t49;
t1 = [(-g(1) * (-t68 * t49 - t70 * t65) - g(2) * (-t68 * t65 + t48)) * MDP(21) + (-g(2) * t48 + (g(1) * t65 - g(2) * t72) * t70 + (-g(1) * (-t49 - t72) + g(2) * t65) * t68) * MDP(25) + (MDP(3) - MDP(20) - MDP(23)) * t54 + (-t67 * MDP(10) + MDP(16) * t62 - MDP(17) * t61 + t69 * MDP(9) - t81 * t57 + t82 * t58 + MDP(2)) * t53; (-g(3) * t69 + t54 * t67) * MDP(9) + (g(3) * t67 + t54 * t69) * MDP(10) + (-g(3) * t76 - t54 * t50) * MDP(21) + (-g(1) * (t73 * t70 + t52) - g(2) * (t73 * t68 + t51) - g(3) * (t63 + t74)) * MDP(25) + t71; (-g(1) * (-t70 * t78 + t52) - g(2) * (-t68 * t78 + t51) - g(3) * t74) * MDP(25) + (-MDP(21) * t77 + t79 * t80) * pkin(3) + t71; -t79 * t53; -t42 * MDP(25);];
taug = t1;
