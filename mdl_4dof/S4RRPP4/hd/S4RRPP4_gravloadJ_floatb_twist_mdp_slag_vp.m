% Calculate Gravitation load on the joints for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RRPP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:18
% EndTime: 2019-12-31 16:59:19
% DurationCPUTime: 0.18s
% Computational Cost: add. (77->39), mult. (162->51), div. (0->0), fcn. (129->4), ass. (0->23)
t66 = MDP(9) + MDP(11) + MDP(15);
t65 = MDP(10) - MDP(13) - MDP(16);
t50 = sin(qJ(1));
t52 = cos(qJ(1));
t38 = g(1) * t52 + g(2) * t50;
t49 = sin(qJ(2));
t64 = t38 * t49;
t43 = t49 * qJ(3);
t51 = cos(qJ(2));
t56 = t51 * pkin(2) + t43;
t62 = pkin(2) * t49;
t61 = g(1) * t50;
t58 = t51 * pkin(3);
t57 = t51 * t52;
t55 = qJ(3) * t51;
t54 = pkin(2) * t57 + t50 * pkin(5) + (pkin(1) + t43) * t52;
t37 = -g(2) * t52 + t61;
t53 = -pkin(1) - t56;
t46 = t52 * pkin(5);
t41 = t52 * t55;
t39 = t50 * t55;
t33 = -g(3) * t51 + t64;
t1 = [(-g(1) * t46 - g(2) * t54 - t53 * t61) * MDP(14) + (-g(1) * (-t52 * qJ(4) + t46) - g(2) * (pkin(3) * t57 + t54) + (-g(1) * (t53 - t58) + g(2) * qJ(4)) * t50) * MDP(18) + (MDP(3) - MDP(12) + MDP(17)) * t38 + (-t65 * t49 + t66 * t51 + MDP(2)) * t37; (-g(1) * (-t52 * t62 + t41) - g(2) * (-t50 * t62 + t39) - g(3) * t56) * MDP(14) + (-g(1) * t41 - g(2) * t39 - g(3) * (t56 + t58) + (pkin(2) + pkin(3)) * t64) * MDP(18) + t65 * (g(3) * t49 + t38 * t51) + t66 * t33; (-MDP(14) - MDP(18)) * t33; t37 * MDP(18);];
taug = t1;
