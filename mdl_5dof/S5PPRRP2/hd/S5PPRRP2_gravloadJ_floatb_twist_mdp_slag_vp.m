% Calculate Gravitation load on the joints for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:10
% EndTime: 2019-12-05 15:09:11
% DurationCPUTime: 0.14s
% Computational Cost: add. (124->35), mult. (172->50), div. (0->0), fcn. (161->6), ass. (0->23)
t44 = sin(pkin(7));
t45 = cos(pkin(7));
t51 = g(1) * t45 + g(2) * t44;
t61 = MDP(11) + MDP(13);
t60 = MDP(12) - MDP(15);
t43 = pkin(8) + qJ(3);
t41 = sin(t43);
t57 = g(3) * t41;
t46 = sin(qJ(4));
t56 = t44 * t46;
t47 = cos(qJ(4));
t55 = t44 * t47;
t54 = t45 * t46;
t53 = t45 * t47;
t52 = MDP(16) + MDP(2);
t50 = pkin(4) * t47 + qJ(5) * t46 + pkin(3);
t42 = cos(t43);
t34 = t42 * t56 + t53;
t36 = t42 * t54 - t55;
t30 = g(1) * t36 + g(2) * t34 + t46 * t57;
t37 = t42 * t53 + t56;
t35 = t42 * t55 - t54;
t1 = [(-MDP(1) - t52) * g(3); t52 * (-g(1) * t44 + g(2) * t45); ((-t51 * pkin(6) - g(3) * t50) * t42 + (-g(3) * pkin(6) + t51 * t50) * t41) * MDP(16) + (MDP(5) - MDP(14)) * (t51 * t42 + t57) + (-t60 * t46 + t61 * t47 + MDP(4)) * (-g(3) * t42 + t51 * t41); (-g(1) * (-t36 * pkin(4) + t37 * qJ(5)) - g(2) * (-t34 * pkin(4) + t35 * qJ(5)) - (-pkin(4) * t46 + qJ(5) * t47) * t57) * MDP(16) + t60 * (g(1) * t37 + g(2) * t35 + t47 * t57) + t61 * t30; -t30 * MDP(16);];
taug = t1;
