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
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:17:25
% EndTime: 2022-01-20 11:17:26
% DurationCPUTime: 0.14s
% Computational Cost: add. (207->47), mult. (191->75), div. (0->0), fcn. (188->10), ass. (0->30)
t83 = g(3) * sin(pkin(9));
t68 = qJ(1) + qJ(2);
t64 = sin(t68);
t70 = cos(pkin(9));
t82 = t64 * t70;
t66 = cos(t68);
t81 = t66 * t70;
t71 = sin(qJ(4));
t80 = t70 * t71;
t73 = cos(qJ(4));
t79 = t70 * t73;
t67 = qJ(4) + qJ(5);
t63 = sin(t67);
t65 = cos(t67);
t46 = t63 * t82 + t66 * t65;
t47 = t66 * t63 - t65 * t82;
t48 = -t63 * t81 + t64 * t65;
t49 = t64 * t63 + t65 * t81;
t78 = (-g(1) * t48 + g(2) * t46 + t63 * t83) * MDP(22) + (g(1) * t49 - g(2) * t47 + t65 * t83) * MDP(23);
t77 = t66 * pkin(2) + t64 * qJ(3);
t76 = -t64 * pkin(2) + t66 * qJ(3);
t58 = g(1) * t64 - g(2) * t66;
t51 = t64 * t80 + t66 * t73;
t52 = -t64 * t79 + t66 * t71;
t53 = t64 * t73 - t66 * t80;
t54 = t64 * t71 + t66 * t79;
t75 = (-g(1) * t46 - g(2) * t48) * MDP(23) + (-g(1) * t47 - g(2) * t49) * MDP(22) + (-g(1) * t51 - g(2) * t53) * MDP(16) + (-g(1) * t52 - g(2) * t54) * MDP(15) + (MDP(6) - MDP(8)) * (g(1) * t66 + g(2) * t64) + (t70 * MDP(7) + MDP(5)) * t58;
t74 = cos(qJ(1));
t72 = sin(qJ(1));
t1 = [(g(1) * t72 - g(2) * t74) * MDP(2) + (g(1) * t74 + g(2) * t72) * MDP(3) + (-g(1) * (-t72 * pkin(1) + t76) - g(2) * (t74 * pkin(1) + t77)) * MDP(9) + t75; (-g(1) * t76 - g(2) * t77) * MDP(9) + t75; -t58 * MDP(9); (-g(1) * t53 + g(2) * t51 + t71 * t83) * MDP(15) + (g(1) * t54 - g(2) * t52 + t73 * t83) * MDP(16) + t78; t78;];
taug = t1;
