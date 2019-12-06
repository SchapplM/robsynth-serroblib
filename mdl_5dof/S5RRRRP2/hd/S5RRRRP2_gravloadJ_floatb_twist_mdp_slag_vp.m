% Calculate Gravitation load on the joints for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:48:13
% EndTime: 2019-12-05 18:48:13
% DurationCPUTime: 0.08s
% Computational Cost: add. (166->37), mult. (148->50), div. (0->0), fcn. (111->8), ass. (0->21)
t51 = qJ(1) + qJ(2);
t45 = sin(t51);
t47 = cos(t51);
t38 = g(2) * t45 - g(3) * t47;
t50 = qJ(3) + qJ(4);
t44 = sin(t50);
t46 = cos(t50);
t56 = -g(1) * t46 - t38 * t44;
t61 = t56 * MDP(19) + (g(1) * t44 - t38 * t46) * MDP(20);
t54 = cos(qJ(3));
t60 = t54 * pkin(3) + pkin(4) * t46;
t40 = pkin(2) + t60;
t49 = -qJ(5) - pkin(8) - pkin(7);
t59 = -t47 * t40 + t45 * t49;
t39 = g(2) * t47 + g(3) * t45;
t58 = -t45 * t40 - t47 * t49;
t52 = sin(qJ(3));
t57 = (MDP(21) - MDP(6)) * t38 + (t54 * MDP(12) - t52 * MDP(13) + MDP(19) * t46 - MDP(20) * t44 + MDP(5)) * t39;
t55 = cos(qJ(1));
t53 = sin(qJ(1));
t1 = [(g(2) * t55 + g(3) * t53) * MDP(2) + (-g(2) * t53 + g(3) * t55) * MDP(3) + (-g(2) * (-t55 * pkin(1) + t59) - g(3) * (-t53 * pkin(1) + t58)) * MDP(22) + t57; (-g(2) * t59 - g(3) * t58) * MDP(22) + t57; (-g(1) * t54 - t38 * t52) * MDP(12) + (g(1) * t52 - t38 * t54) * MDP(13) + (-g(1) * t60 + t38 * (-t52 * pkin(3) - pkin(4) * t44)) * MDP(22) + t61; t56 * MDP(22) * pkin(4) + t61; -t39 * MDP(22);];
taug = t1;
