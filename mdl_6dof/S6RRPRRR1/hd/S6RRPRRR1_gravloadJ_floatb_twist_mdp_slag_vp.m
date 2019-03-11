% Calculate Gravitation load on the joints for
% S6RRPRRR1
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
%   see S6RRPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:14:29
% EndTime: 2019-03-09 13:14:30
% DurationCPUTime: 0.18s
% Computational Cost: add. (279->46), mult. (226->64), div. (0->0), fcn. (195->10), ass. (0->28)
t52 = qJ(2) + pkin(11) + qJ(4);
t50 = qJ(5) + t52;
t46 = sin(t50);
t70 = g(3) * t46;
t54 = sin(qJ(6));
t59 = cos(qJ(1));
t68 = t54 * t59;
t56 = sin(qJ(1));
t67 = t56 * t54;
t57 = cos(qJ(6));
t66 = t56 * t57;
t65 = t57 * t59;
t47 = cos(t50);
t63 = g(1) * t59 + g(2) * t56;
t64 = (t47 * t63 + t70) * MDP(26) + (t57 * MDP(32) - t54 * MDP(33) + MDP(25)) * (-g(3) * t47 + t46 * t63);
t44 = g(1) * t56 - g(2) * t59;
t48 = sin(t52);
t49 = cos(t52);
t62 = (-g(3) * t49 + t48 * t63) * MDP(18) + (g(3) * t48 + t49 * t63) * MDP(19) + t64;
t58 = cos(qJ(2));
t55 = sin(qJ(2));
t53 = -qJ(3) - pkin(7);
t51 = pkin(2) * t58 + pkin(1);
t43 = t47 * t65 + t67;
t42 = -t47 * t68 + t66;
t41 = -t47 * t66 + t68;
t40 = t47 * t67 + t65;
t1 = [(-g(1) * (-t56 * t51 - t53 * t59) - g(2) * (t51 * t59 - t56 * t53)) * MDP(12) + (-g(1) * t41 - g(2) * t43) * MDP(32) + (-g(1) * t40 - g(2) * t42) * MDP(33) + (MDP(3) - MDP(11)) * t63 + (-MDP(10) * t55 + MDP(18) * t49 - MDP(19) * t48 + MDP(25) * t47 - MDP(26) * t46 + MDP(9) * t58 + MDP(2)) * t44; (g(3) * t55 + t58 * t63) * MDP(10) + t62 + (MDP(12) * pkin(2) + MDP(9)) * (-g(3) * t58 + t55 * t63); -t44 * MDP(12); t62; t64; (-g(1) * t42 + g(2) * t40 + t54 * t70) * MDP(32) + (g(1) * t43 - g(2) * t41 + t57 * t70) * MDP(33);];
taug  = t1;
