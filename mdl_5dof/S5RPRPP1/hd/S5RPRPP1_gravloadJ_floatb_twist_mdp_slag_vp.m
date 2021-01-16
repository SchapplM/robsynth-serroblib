% Calculate Gravitation load on the joints for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:14
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:14:04
% EndTime: 2021-01-15 11:14:06
% DurationCPUTime: 0.25s
% Computational Cost: add. (168->48), mult. (156->58), div. (0->0), fcn. (117->8), ass. (0->23)
t70 = MDP(12) + MDP(16);
t69 = MDP(13) - MDP(18);
t48 = qJ(3) + pkin(8);
t42 = sin(t48);
t44 = cos(t48);
t68 = t44 * pkin(4) + t42 * qJ(5);
t49 = qJ(1) + pkin(7);
t43 = sin(t49);
t45 = cos(t49);
t61 = g(1) * t45 + g(2) * t43;
t53 = cos(qJ(3));
t46 = t53 * pkin(3);
t41 = t46 + pkin(2);
t54 = cos(qJ(1));
t64 = t54 * pkin(1) + t45 * t41;
t62 = MDP(15) + MDP(19);
t60 = g(1) * t43 - g(2) * t45;
t50 = -qJ(4) - pkin(6);
t52 = sin(qJ(1));
t58 = -t52 * pkin(1) - t45 * t50;
t51 = sin(qJ(3));
t34 = -g(3) * t44 + t61 * t42;
t1 = [(g(1) * t54 + g(2) * t52) * MDP(3) + (-g(1) * (-t43 * t41 + t58) - g(2) * (-t43 * t50 + t64)) * MDP(15) + (-g(1) * t58 - g(2) * (t68 * t45 + t64) + (-g(1) * (-t41 - t68) + g(2) * t50) * t43) * MDP(19) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t52 - g(2) * t54) - (MDP(14) + MDP(17)) * t61 + (t53 * MDP(10) - t51 * MDP(11) - t69 * t42 + t70 * t44) * t60; (-MDP(4) - t62) * g(3); (g(3) * t51 + t61 * t53) * MDP(11) + (-g(3) * (t46 + t68) + t61 * (pkin(3) * t51 + pkin(4) * t42 - qJ(5) * t44)) * MDP(19) + (pkin(3) * MDP(15) + MDP(10)) * (-g(3) * t53 + t61 * t51) + t69 * (g(3) * t42 + t61 * t44) + t70 * t34; -t62 * t60; -t34 * MDP(19);];
taug = t1;
