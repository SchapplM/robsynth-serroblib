% Calculate Gravitation load on the joints for
% S5RPRRP6
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
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 18:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 18:08:28
% EndTime: 2021-01-15 18:08:29
% DurationCPUTime: 0.22s
% Computational Cost: add. (174->46), mult. (224->66), div. (0->0), fcn. (203->8), ass. (0->29)
t58 = qJ(1) + pkin(8);
t56 = sin(t58);
t57 = cos(t58);
t72 = g(1) * t57 + g(2) * t56;
t86 = MDP(11) - MDP(21);
t85 = MDP(17) + MDP(19);
t84 = MDP(18) + MDP(20);
t61 = sin(qJ(3));
t64 = cos(qJ(3));
t48 = -g(3) * t64 + t72 * t61;
t78 = g(3) * t61;
t60 = sin(qJ(4));
t76 = t60 * t64;
t63 = cos(qJ(4));
t75 = t63 * t64;
t73 = pkin(4) * t60 + pkin(6);
t62 = sin(qJ(1));
t65 = cos(qJ(1));
t70 = g(1) * t62 - g(2) * t65;
t55 = pkin(4) * t63 + pkin(3);
t59 = -qJ(5) - pkin(7);
t69 = t55 * t64 - t59 * t61;
t67 = pkin(2) + t69;
t52 = t56 * t63 - t57 * t76;
t50 = t56 * t76 + t57 * t63;
t66 = t70 * pkin(1);
t53 = t56 * t60 + t57 * t75;
t51 = -t56 * t75 + t57 * t60;
t1 = [t70 * MDP(2) + (g(1) * t65 + g(2) * t62) * MDP(3) + MDP(4) * t66 + (t66 + (-g(1) * t73 - g(2) * t67) * t57 + (g(1) * t67 - g(2) * t73) * t56) * MDP(22) + t85 * (-g(1) * t51 - g(2) * t53) + t84 * (-g(1) * t50 - g(2) * t52) + (t64 * MDP(10) - t86 * t61) * (g(1) * t56 - g(2) * t57); (-MDP(22) - MDP(4)) * g(3); (-g(3) * t69 + t72 * (t55 * t61 + t59 * t64)) * MDP(22) + t86 * (t72 * t64 + t78) + (-t84 * t60 + t85 * t63 + MDP(10)) * t48; t84 * (g(1) * t53 - g(2) * t51 + t63 * t78) + (pkin(4) * MDP(22) + t85) * (-g(1) * t52 + g(2) * t50 + t60 * t78); -t48 * MDP(22);];
taug = t1;
