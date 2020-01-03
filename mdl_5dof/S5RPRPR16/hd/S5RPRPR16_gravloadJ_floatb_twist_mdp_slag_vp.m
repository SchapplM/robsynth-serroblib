% Calculate Gravitation load on the joints for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR16_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR16_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR16_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:33
% EndTime: 2019-12-31 18:39:34
% DurationCPUTime: 0.22s
% Computational Cost: add. (76->45), mult. (167->63), div. (0->0), fcn. (143->6), ass. (0->22)
t50 = sin(qJ(3));
t53 = cos(qJ(3));
t56 = t50 * pkin(3) - t53 * qJ(4);
t51 = sin(qJ(1));
t54 = cos(qJ(1));
t71 = -g(1) * t51 + g(2) * t54;
t70 = -MDP(12) + MDP(15);
t69 = MDP(13) - MDP(16);
t65 = g(3) * t50;
t62 = t51 * t53;
t49 = sin(qJ(5));
t61 = t54 * t49;
t52 = cos(qJ(5));
t60 = t54 * t52;
t59 = t54 * pkin(1) + t51 * qJ(2);
t46 = t54 * qJ(2);
t40 = -t49 * t62 + t60;
t39 = -t52 * t62 - t61;
t38 = -t51 * t52 - t53 * t61;
t37 = t51 * t49 - t53 * t60;
t36 = -t53 * t71 - t65;
t1 = [(-g(1) * (-t51 * pkin(1) + t46) - g(2) * t59) * MDP(6) + (-g(1) * (t56 * t54 + t46) - g(2) * (t54 * pkin(6) + t59) + (-g(1) * (-pkin(1) - pkin(6)) - g(2) * t56) * t51) * MDP(17) + (-g(1) * t38 - g(2) * t40) * MDP(23) + (-g(1) * t37 - g(2) * t39) * MDP(24) - (MDP(2) - MDP(4) + MDP(14)) * t71 + (t70 * t50 - t69 * t53 + MDP(3) - MDP(5)) * (g(1) * t54 + g(2) * t51); -(-MDP(17) - MDP(6)) * t71; (g(3) * t56 + t71 * (pkin(3) * t53 + qJ(4) * t50)) * MDP(17) + t70 * t36 + (-MDP(23) * t49 - MDP(24) * t52 + t69) * (g(3) * t53 - t50 * t71); t36 * MDP(17); (-g(1) * t39 + g(2) * t37 - t52 * t65) * MDP(23) + (g(1) * t40 - g(2) * t38 + t49 * t65) * MDP(24);];
taug = t1;
