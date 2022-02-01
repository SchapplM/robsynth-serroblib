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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 12:08:23
% EndTime: 2022-01-20 12:08:23
% DurationCPUTime: 0.08s
% Computational Cost: add. (195->27), mult. (150->36), div. (0->0), fcn. (120->10), ass. (0->18)
t45 = qJ(3) + qJ(4);
t44 = qJ(5) + t45;
t38 = sin(t44);
t39 = cos(t44);
t46 = qJ(1) + qJ(2);
t41 = sin(t46);
t43 = cos(t46);
t53 = g(1) * t43 + g(2) * t41;
t55 = (-g(3) * t39 + t38 * t53) * MDP(26) + (g(3) * t38 + t39 * t53) * MDP(27);
t40 = sin(t45);
t42 = cos(t45);
t54 = (-g(3) * t42 + t40 * t53) * MDP(19) + (g(3) * t40 + t42 * t53) * MDP(20) + t55;
t47 = sin(qJ(3));
t49 = cos(qJ(3));
t51 = t53 * MDP(6) + (MDP(12) * t49 - MDP(13) * t47 + MDP(19) * t42 - MDP(20) * t40 + MDP(26) * t39 - MDP(27) * t38 + MDP(5)) * (g(1) * t41 - g(2) * t43);
t50 = cos(qJ(1));
t48 = sin(qJ(1));
t1 = [(g(1) * t48 - g(2) * t50) * MDP(2) + (g(1) * t50 + g(2) * t48) * MDP(3) + t51; t51; (-g(3) * t49 + t47 * t53) * MDP(12) + (g(3) * t47 + t49 * t53) * MDP(13) + t54; t54; t55;];
taug = t1;
