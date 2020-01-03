% Calculate Gravitation load on the joints for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S4RPRR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:29
% EndTime: 2019-12-31 16:56:29
% DurationCPUTime: 0.09s
% Computational Cost: add. (44->30), mult. (102->46), div. (0->0), fcn. (92->6), ass. (0->17)
t33 = cos(qJ(3));
t41 = g(3) * t33;
t29 = sin(qJ(4));
t31 = sin(qJ(1));
t40 = t31 * t29;
t32 = cos(qJ(4));
t39 = t31 * t32;
t34 = cos(qJ(1));
t38 = t34 * t29;
t37 = t34 * t32;
t26 = g(1) * t31 - g(2) * t34;
t30 = sin(qJ(3));
t25 = t30 * t37 - t40;
t24 = t30 * t38 + t39;
t23 = t30 * t39 + t38;
t22 = -t30 * t40 + t37;
t1 = [(-g(1) * (-t31 * pkin(1) + t34 * qJ(2)) - g(2) * (t34 * pkin(1) + t31 * qJ(2))) * MDP(6) + (-g(1) * t25 - g(2) * t23) * MDP(19) + (g(1) * t24 - g(2) * t22) * MDP(20) + (MDP(2) - MDP(4)) * t26 + (-t30 * MDP(12) - t33 * MDP(13) + MDP(3) - MDP(5)) * (g(1) * t34 + g(2) * t31); -t26 * MDP(6); (t26 * t30 + t41) * MDP(13) + (-MDP(19) * t32 + MDP(20) * t29 - MDP(12)) * (-g(3) * t30 + t26 * t33); (-g(1) * t22 - g(2) * t24 + t29 * t41) * MDP(19) + (g(1) * t23 - g(2) * t25 + t32 * t41) * MDP(20);];
taug = t1;
