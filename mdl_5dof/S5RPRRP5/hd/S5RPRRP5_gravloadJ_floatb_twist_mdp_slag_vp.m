% Calculate Gravitation load on the joints for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:57
% EndTime: 2019-12-31 18:40:58
% DurationCPUTime: 0.09s
% Computational Cost: add. (190->36), mult. (149->44), div. (0->0), fcn. (113->8), ass. (0->20)
t54 = sin(qJ(4));
t56 = cos(qJ(4));
t61 = t56 * pkin(4) + t54 * qJ(5);
t72 = -pkin(3) - t61;
t71 = MDP(13) + MDP(15);
t70 = MDP(14) - MDP(17);
t53 = qJ(1) + pkin(8);
t52 = qJ(3) + t53;
t50 = sin(t52);
t51 = cos(t52);
t44 = g(1) * t51 + g(2) * t50;
t69 = g(1) * t50;
t64 = t50 * pkin(7) - t72 * t51;
t59 = (-MDP(16) + MDP(7)) * t44 + (-t70 * t54 + t71 * t56 + MDP(6)) * (-g(2) * t51 + t69);
t58 = t72 * t69;
t57 = cos(qJ(1));
t55 = sin(qJ(1));
t48 = t51 * pkin(7);
t33 = -g(3) * t56 + t44 * t54;
t1 = [(g(1) * t57 + g(2) * t55) * MDP(3) + (-g(1) * (-pkin(2) * sin(t53) - t55 * pkin(1) + t48) - g(2) * (pkin(2) * cos(t53) + t57 * pkin(1) + t64) - t58) * MDP(18) + t59 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t55 - g(2) * t57); (-MDP(18) - MDP(4)) * g(3); (-g(1) * t48 - g(2) * t64 - t58) * MDP(18) + t59; (-g(3) * t61 + t44 * (pkin(4) * t54 - qJ(5) * t56)) * MDP(18) + t70 * (g(3) * t54 + t44 * t56) + t71 * t33; -t33 * MDP(18);];
taug = t1;
