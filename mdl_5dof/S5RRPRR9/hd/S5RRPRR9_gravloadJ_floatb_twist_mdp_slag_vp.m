% Calculate Gravitation load on the joints for
% S5RRPRR9
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
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:47:34
% EndTime: 2021-01-15 21:47:35
% DurationCPUTime: 0.22s
% Computational Cost: add. (174->55), mult. (212->78), div. (0->0), fcn. (199->10), ass. (0->37)
t58 = qJ(4) + qJ(5);
t55 = sin(t58);
t56 = cos(t58);
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t83 = t63 * MDP(20) - t60 * MDP(21) + t56 * MDP(27) - t55 * MDP(28) + MDP(11);
t57 = qJ(2) + pkin(9);
t53 = sin(t57);
t61 = sin(qJ(2));
t82 = t61 * MDP(10) + MDP(12) * t53;
t81 = g(3) * t53;
t65 = cos(qJ(1));
t80 = t55 * t65;
t79 = t56 * t65;
t78 = t60 * t65;
t62 = sin(qJ(1));
t77 = t62 * t55;
t76 = t62 * t56;
t75 = t62 * t60;
t74 = t62 * t63;
t73 = t63 * t65;
t54 = cos(t57);
t42 = t54 * t77 + t79;
t43 = -t54 * t76 + t80;
t44 = -t54 * t80 + t76;
t45 = t54 * t79 + t77;
t72 = (-g(1) * t44 + g(2) * t42 + t55 * t81) * MDP(27) + (g(1) * t45 - g(2) * t43 + t56 * t81) * MDP(28);
t67 = g(1) * t65 + g(2) * t62;
t50 = g(1) * t62 - g(2) * t65;
t64 = cos(qJ(2));
t59 = -qJ(3) - pkin(6);
t52 = pkin(2) * t64 + pkin(1);
t49 = t54 * t73 + t75;
t48 = -t54 * t78 + t74;
t47 = -t54 * t74 + t78;
t46 = t54 * t75 + t73;
t1 = [(-g(1) * (-t62 * t52 - t59 * t65) - g(2) * (t52 * t65 - t62 * t59)) * MDP(14) + (-g(1) * t47 - g(2) * t49) * MDP(20) + (-g(1) * t46 - g(2) * t48) * MDP(21) + (-g(1) * t43 - g(2) * t45) * MDP(27) + (-g(1) * t42 - g(2) * t44) * MDP(28) + (MDP(3) - MDP(13)) * t67 + (MDP(11) * t54 + MDP(9) * t64 + MDP(2) - t82) * t50; (MDP(14) * pkin(2) + MDP(9)) * (-g(3) * t64 + t67 * t61) + (-t83 * t54 + t82) * g(3) + (t64 * MDP(10) + MDP(12) * t54 + t83 * t53) * t67; -t50 * MDP(14); (-g(1) * t48 + g(2) * t46 + t60 * t81) * MDP(20) + (g(1) * t49 - g(2) * t47 + t63 * t81) * MDP(21) + t72; t72;];
taug = t1;
