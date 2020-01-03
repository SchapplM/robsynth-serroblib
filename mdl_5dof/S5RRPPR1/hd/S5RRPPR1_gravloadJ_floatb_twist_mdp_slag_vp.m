% Calculate Gravitation load on the joints for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:56:00
% EndTime: 2020-01-03 11:56:00
% DurationCPUTime: 0.07s
% Computational Cost: add. (155->36), mult. (113->49), div. (0->0), fcn. (80->10), ass. (0->22)
t54 = qJ(1) + qJ(2);
t48 = pkin(8) + t54;
t42 = sin(t48);
t43 = cos(t48);
t50 = cos(t54);
t45 = pkin(2) * t50;
t64 = t43 * pkin(3) + t42 * qJ(4) + t45;
t49 = sin(t54);
t44 = pkin(2) * t49;
t63 = t42 * pkin(3) - t43 * qJ(4) + t44;
t62 = g(2) * t43 + g(3) * t42;
t61 = g(2) * t42 - g(3) * t43;
t60 = -g(2) * t50 - g(3) * t49;
t53 = pkin(9) + qJ(5);
t46 = sin(t53);
t47 = cos(t53);
t59 = -t61 * MDP(10) + (g(2) * t49 - g(3) * t50) * MDP(6) + t60 * MDP(5) + (-t47 * MDP(17) + t46 * MDP(18) - cos(pkin(9)) * MDP(8) + sin(pkin(9)) * MDP(9)) * t62;
t58 = cos(qJ(1));
t57 = sin(qJ(1));
t52 = t58 * pkin(1);
t51 = t57 * pkin(1);
t1 = [(-g(2) * t58 - g(3) * t57) * MDP(2) + (g(2) * t57 - g(3) * t58) * MDP(3) + (-g(2) * (t45 + t52) - g(3) * (t44 + t51)) * MDP(7) + (-g(2) * (t52 + t64) - g(3) * (t51 + t63)) * MDP(11) + t59; (-g(2) * t64 - g(3) * t63) * MDP(11) + t60 * MDP(7) * pkin(2) + t59; (-MDP(11) - MDP(7)) * g(1); t62 * MDP(11); (-g(1) * t47 + t61 * t46) * MDP(17) + (g(1) * t46 + t61 * t47) * MDP(18);];
taug = t1;
