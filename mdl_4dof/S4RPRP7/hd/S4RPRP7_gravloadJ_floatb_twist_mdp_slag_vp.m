% Calculate Gravitation load on the joints for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RPRP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:17
% EndTime: 2019-12-31 16:47:18
% DurationCPUTime: 0.13s
% Computational Cost: add. (53->30), mult. (109->38), div. (0->0), fcn. (83->4), ass. (0->12)
t37 = sin(qJ(3));
t39 = cos(qJ(3));
t41 = t37 * pkin(3) - t39 * qJ(4);
t38 = sin(qJ(1));
t40 = cos(qJ(1));
t50 = -g(1) * t38 + g(2) * t40;
t49 = MDP(12) + MDP(14);
t48 = MDP(13) - MDP(16);
t44 = t40 * pkin(1) + t38 * qJ(2);
t33 = t40 * qJ(2);
t26 = -g(3) * t37 - t50 * t39;
t1 = [(-g(1) * (-t38 * pkin(1) + t33) - g(2) * t44) * MDP(6) + (-g(1) * (t41 * t40 + t33) - g(2) * (t40 * pkin(5) + t44) + (-g(1) * (-pkin(1) - pkin(5)) - g(2) * t41) * t38) * MDP(17) - (MDP(2) - MDP(4) + MDP(15)) * t50 + (-t49 * t37 - t48 * t39 + MDP(3) - MDP(5)) * (g(1) * t40 + g(2) * t38); -(-MDP(17) - MDP(6)) * t50; (g(3) * t41 + t50 * (pkin(3) * t39 + qJ(4) * t37)) * MDP(17) - t49 * t26 + t48 * (g(3) * t39 - t37 * t50); t26 * MDP(17);];
taug = t1;
