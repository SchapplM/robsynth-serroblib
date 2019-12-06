% Calculate Gravitation load on the joints for
% S5RRPRR5
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
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:34:04
% EndTime: 2019-12-05 18:34:04
% DurationCPUTime: 0.07s
% Computational Cost: add. (177->33), mult. (139->43), div. (0->0), fcn. (106->10), ass. (0->18)
t50 = qJ(1) + qJ(2);
t47 = sin(t50);
t48 = cos(t50);
t39 = g(2) * t47 - g(3) * t48;
t49 = pkin(9) + qJ(4);
t46 = qJ(5) + t49;
t42 = sin(t46);
t43 = cos(t46);
t58 = (-g(1) * t43 - t39 * t42) * MDP(23) + (g(1) * t42 - t39 * t43) * MDP(24);
t57 = -t47 * pkin(2) + t48 * qJ(3);
t40 = g(2) * t48 + g(3) * t47;
t56 = -t48 * pkin(2) - t47 * qJ(3);
t44 = sin(t49);
t45 = cos(t49);
t55 = (-MDP(6) + MDP(9)) * t39 + (t45 * MDP(16) - t44 * MDP(17) + t43 * MDP(23) - t42 * MDP(24) + MDP(7) * cos(pkin(9)) - MDP(8) * sin(pkin(9)) + MDP(5)) * t40;
t54 = cos(qJ(1));
t53 = sin(qJ(1));
t1 = [(g(2) * t54 + g(3) * t53) * MDP(2) + (-g(2) * t53 + g(3) * t54) * MDP(3) + (-g(2) * (-t54 * pkin(1) + t56) - g(3) * (-t53 * pkin(1) + t57)) * MDP(10) + t55; (-g(2) * t56 - g(3) * t57) * MDP(10) + t55; -t40 * MDP(10); (-g(1) * t45 - t39 * t44) * MDP(16) + (g(1) * t44 - t39 * t45) * MDP(17) + t58; t58;];
taug = t1;
