% Calculate Gravitation load on the joints for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:32
% EndTime: 2019-12-31 20:01:33
% DurationCPUTime: 0.38s
% Computational Cost: add. (195->65), mult. (285->88), div. (0->0), fcn. (269->8), ass. (0->32)
t87 = MDP(18) + MDP(20);
t86 = MDP(19) - MDP(22);
t62 = qJ(2) + pkin(8);
t59 = cos(t62);
t58 = sin(t62);
t82 = t58 * pkin(7);
t88 = t59 * pkin(3) + t82;
t66 = sin(qJ(1));
t69 = cos(qJ(1));
t53 = g(1) * t69 + g(2) * t66;
t83 = g(3) * t58;
t64 = sin(qJ(4));
t80 = t66 * t64;
t67 = cos(qJ(4));
t79 = t66 * t67;
t63 = -qJ(3) - pkin(6);
t78 = t69 * t63;
t77 = t69 * t64;
t76 = t69 * t67;
t52 = g(1) * t66 - g(2) * t69;
t74 = pkin(4) * t67 + qJ(5) * t64 + pkin(3);
t48 = t59 * t80 + t76;
t50 = t59 * t77 - t79;
t44 = g(1) * t50 + g(2) * t48 + t64 * t83;
t68 = cos(qJ(2));
t65 = sin(qJ(2));
t60 = t68 * pkin(2);
t57 = t60 + pkin(1);
t54 = t69 * t57;
t51 = t59 * t76 + t80;
t49 = t59 * t79 - t77;
t1 = [(-g(1) * (-t66 * t57 - t78) - g(2) * (-t66 * t63 + t54)) * MDP(12) + (-g(1) * (-t49 * pkin(4) - t48 * qJ(5) - t78) - g(2) * (t51 * pkin(4) + t50 * qJ(5) + t88 * t69 + t54) + (-g(1) * (-t57 - t88) + g(2) * t63) * t66) * MDP(23) - t86 * (g(1) * t48 - g(2) * t50) + (MDP(3) - MDP(11)) * t53 + t87 * (g(1) * t49 - g(2) * t51) + (-t65 * MDP(10) + t58 * MDP(21) + t68 * MDP(9) + MDP(2)) * t52; (g(3) * t65 + t53 * t68) * MDP(10) + (-t53 * t59 - t83) * MDP(21) + (-g(3) * (t74 * t59 + t60 + t82) + t53 * (pkin(2) * t65 - pkin(7) * t59 + t74 * t58)) * MDP(23) + (pkin(2) * MDP(12) + MDP(9)) * (-g(3) * t68 + t53 * t65) + (-t86 * t64 + t87 * t67) * (-g(3) * t59 + t53 * t58); (-MDP(12) - MDP(23)) * t52; (-g(1) * (-t50 * pkin(4) + t51 * qJ(5)) - g(2) * (-t48 * pkin(4) + t49 * qJ(5)) - (-pkin(4) * t64 + qJ(5) * t67) * t83) * MDP(23) + t86 * (g(1) * t51 + g(2) * t49 + t67 * t83) + t87 * t44; -t44 * MDP(23);];
taug = t1;
