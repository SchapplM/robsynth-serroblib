% Calculate Gravitation load on the joints for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:15:19
% EndTime: 2020-01-03 12:15:19
% DurationCPUTime: 0.09s
% Computational Cost: add. (195->27), mult. (150->36), div. (0->0), fcn. (120->10), ass. (0->18)
t46 = qJ(3) + qJ(4);
t45 = qJ(5) + t46;
t39 = sin(t45);
t40 = cos(t45);
t47 = qJ(1) + qJ(2);
t42 = sin(t47);
t44 = cos(t47);
t53 = g(2) * t42 - g(3) * t44;
t56 = (-g(1) * t40 + t53 * t39) * MDP(26) + (g(1) * t39 + t53 * t40) * MDP(27);
t41 = sin(t46);
t43 = cos(t46);
t55 = (-g(1) * t43 + t53 * t41) * MDP(19) + (g(1) * t41 + t53 * t43) * MDP(20) + t56;
t48 = sin(qJ(3));
t50 = cos(qJ(3));
t52 = t53 * MDP(6) + (-t50 * MDP(12) + t48 * MDP(13) - t43 * MDP(19) + t41 * MDP(20) - t40 * MDP(26) + t39 * MDP(27) - MDP(5)) * (g(2) * t44 + g(3) * t42);
t51 = cos(qJ(1));
t49 = sin(qJ(1));
t1 = [(-g(2) * t51 - g(3) * t49) * MDP(2) + (g(2) * t49 - g(3) * t51) * MDP(3) + t52; t52; (-g(1) * t50 + t53 * t48) * MDP(12) + (g(1) * t48 + t53 * t50) * MDP(13) + t55; t55; t56;];
taug = t1;
