% Calculate Gravitation load on the joints for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:27
% EndTime: 2019-12-05 16:04:29
% DurationCPUTime: 0.22s
% Computational Cost: add. (120->56), mult. (309->102), div. (0->0), fcn. (360->10), ass. (0->29)
t52 = sin(pkin(5));
t72 = g(3) * t52;
t56 = sin(qJ(4));
t71 = t52 * t56;
t59 = cos(qJ(4));
t70 = t52 * t59;
t60 = cos(qJ(2));
t69 = t52 * t60;
t54 = cos(pkin(5));
t57 = sin(qJ(2));
t68 = t54 * t57;
t67 = t54 * t60;
t55 = sin(qJ(5));
t66 = t55 * t56;
t65 = t55 * t57;
t58 = cos(qJ(5));
t64 = t56 * t58;
t63 = t57 * t58;
t51 = sin(pkin(9));
t53 = cos(pkin(9));
t43 = t51 * t57 - t53 * t67;
t45 = t51 * t67 + t53 * t57;
t37 = -g(1) * t45 - g(2) * t43 + g(3) * t69;
t48 = t54 * t59 - t56 * t69;
t46 = -t51 * t68 + t53 * t60;
t44 = t51 * t60 + t53 * t68;
t41 = -t43 * t56 + t53 * t70;
t39 = t45 * t56 + t51 * t70;
t1 = [(-MDP(1) - MDP(7)) * g(3); (-g(1) * (-t45 * pkin(2) + t46 * qJ(3)) - g(2) * (-t43 * pkin(2) + t44 * qJ(3)) - (pkin(2) * t60 + qJ(3) * t57) * t72) * MDP(7) + (-g(1) * (-t45 * t55 + t46 * t64) - g(2) * (-t43 * t55 + t44 * t64) - (t55 * t60 + t56 * t63) * t72) * MDP(20) + (-g(1) * (-t45 * t58 - t46 * t66) - g(2) * (-t43 * t58 - t44 * t66) - (-t56 * t65 + t58 * t60) * t72) * MDP(21) + (-MDP(3) + MDP(5)) * t37 + (-t56 * MDP(13) - t59 * MDP(14) + MDP(4) - MDP(6)) * (g(1) * t46 + g(2) * t44 + t57 * t72); t37 * MDP(7); (g(1) * t39 - g(2) * t41 + g(3) * t48) * MDP(14) + (-MDP(20) * t58 + MDP(21) * t55 - MDP(13)) * (g(1) * (t45 * t59 - t51 * t71) + g(2) * (t43 * t59 + t53 * t71) + g(3) * (-t54 * t56 - t59 * t69)); (-g(1) * (-t39 * t55 + t46 * t58) - g(2) * (t41 * t55 + t44 * t58) - g(3) * (-t48 * t55 + t52 * t63)) * MDP(20) + (-g(1) * (-t39 * t58 - t46 * t55) - g(2) * (t41 * t58 - t44 * t55) - g(3) * (-t48 * t58 - t52 * t65)) * MDP(21);];
taug = t1;
