% Calculate Gravitation load on the joints for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:12
% EndTime: 2019-12-31 18:22:12
% DurationCPUTime: 0.25s
% Computational Cost: add. (157->51), mult. (177->81), div. (0->0), fcn. (158->10), ass. (0->29)
t49 = qJ(1) + pkin(8);
t45 = sin(t49);
t47 = cos(t49);
t63 = g(1) * t47 + g(2) * t45;
t72 = MDP(11) - MDP(14);
t52 = sin(qJ(3));
t54 = cos(qJ(3));
t41 = -g(3) * t54 + t63 * t52;
t69 = g(3) * t52;
t67 = t45 * t54;
t66 = t47 * t54;
t50 = sin(pkin(9));
t65 = t50 * t54;
t51 = cos(pkin(9));
t64 = t51 * t54;
t53 = sin(qJ(1));
t55 = cos(qJ(1));
t61 = g(1) * t53 - g(2) * t55;
t60 = t54 * pkin(3) + t52 * qJ(4);
t58 = pkin(2) + t60;
t57 = t61 * pkin(1);
t48 = pkin(9) + qJ(5);
t46 = cos(t48);
t44 = sin(t48);
t40 = t45 * t44 + t46 * t66;
t39 = -t44 * t66 + t45 * t46;
t38 = t47 * t44 - t46 * t67;
t37 = t44 * t67 + t47 * t46;
t1 = [t61 * MDP(2) + (g(1) * t55 + g(2) * t53) * MDP(3) + MDP(4) * t57 + (-g(1) * (-t45 * t64 + t47 * t50) - g(2) * (t45 * t50 + t47 * t64)) * MDP(12) + (-g(1) * (t45 * t65 + t47 * t51) - g(2) * (t45 * t51 - t47 * t65)) * MDP(13) + (t57 + (-g(1) * pkin(6) - g(2) * t58) * t47 + (-g(2) * pkin(6) + g(1) * t58) * t45) * MDP(15) + (-g(1) * t38 - g(2) * t40) * MDP(21) + (-g(1) * t37 - g(2) * t39) * MDP(22) + (t54 * MDP(10) - t72 * t52) * (g(1) * t45 - g(2) * t47); (-MDP(15) - MDP(4)) * g(3); (-g(3) * t60 + t63 * (pkin(3) * t52 - qJ(4) * t54)) * MDP(15) + t72 * (t63 * t54 + t69) + (MDP(12) * t51 - MDP(13) * t50 + MDP(21) * t46 - MDP(22) * t44 + MDP(10)) * t41; -t41 * MDP(15); (-g(1) * t39 + g(2) * t37 + t44 * t69) * MDP(21) + (g(1) * t40 - g(2) * t38 + t46 * t69) * MDP(22);];
taug = t1;
