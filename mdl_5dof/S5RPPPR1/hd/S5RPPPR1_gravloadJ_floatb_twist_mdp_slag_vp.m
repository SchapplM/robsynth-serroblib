% Calculate Gravitation load on the joints for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:33
% EndTime: 2020-01-03 11:20:34
% DurationCPUTime: 0.18s
% Computational Cost: add. (129->47), mult. (128->77), div. (0->0), fcn. (113->10), ass. (0->28)
t52 = sin(pkin(8));
t68 = g(1) * t52;
t50 = qJ(1) + pkin(7);
t44 = sin(t50);
t54 = cos(pkin(8));
t67 = t44 * t54;
t46 = cos(t50);
t66 = t46 * t54;
t51 = sin(pkin(9));
t65 = t51 * t54;
t53 = cos(pkin(9));
t64 = t53 * t54;
t63 = MDP(12) + MDP(8);
t56 = cos(qJ(1));
t62 = t56 * pkin(1) + t46 * pkin(2) + t44 * qJ(3);
t55 = sin(qJ(1));
t61 = t55 * pkin(1) + t44 * pkin(2) - t46 * qJ(3);
t60 = g(2) * t46 + g(3) * t44;
t59 = -g(2) * t44 + g(3) * t46;
t57 = pkin(3) * t54 + qJ(4) * t52;
t49 = pkin(9) + qJ(5);
t45 = cos(t49);
t43 = sin(t49);
t37 = t44 * t43 + t45 * t66;
t36 = t43 * t66 - t44 * t45;
t35 = -t46 * t43 + t45 * t67;
t34 = -t43 * t67 - t46 * t45;
t1 = [(g(2) * t55 - g(3) * t56) * MDP(3) + t59 * MDP(7) + (-g(2) * t62 - g(3) * t61) * MDP(8) + (-g(2) * (t44 * t51 + t46 * t64) - g(3) * (t44 * t64 - t46 * t51)) * MDP(9) + (-g(2) * (t44 * t53 - t46 * t65) - g(3) * (-t44 * t65 - t46 * t53)) * MDP(10) + (-g(2) * (t57 * t46 + t62) - g(3) * (t57 * t44 + t61)) * MDP(12) + (-g(2) * t37 - g(3) * t35) * MDP(18) + (g(2) * t36 - g(3) * t34) * MDP(19) + (-t54 * MDP(5) + (MDP(6) - MDP(11)) * t52) * t60 + (MDP(4) * pkin(1) + MDP(2)) * (-g(2) * t56 - g(3) * t55); (-MDP(4) - t63) * g(1); t63 * t60; (g(1) * t54 + t59 * t52) * MDP(12); (-g(2) * t34 - g(3) * t36 + t43 * t68) * MDP(18) + (g(2) * t35 - g(3) * t37 + t45 * t68) * MDP(19);];
taug = t1;
