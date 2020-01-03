% Calculate Gravitation load on the joints for
% S4RPRP6
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
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:14
% EndTime: 2019-12-31 16:46:14
% DurationCPUTime: 0.07s
% Computational Cost: add. (39->26), mult. (75->32), div. (0->0), fcn. (53->4), ass. (0->11)
t22 = sin(qJ(3));
t28 = pkin(3) * t22;
t23 = sin(qJ(1));
t25 = cos(qJ(1));
t27 = t25 * pkin(1) + t23 * qJ(2);
t16 = g(1) * t25 + g(2) * t23;
t15 = g(1) * t23 - g(2) * t25;
t24 = cos(qJ(3));
t21 = -qJ(4) - pkin(5);
t18 = t25 * qJ(2);
t1 = [(-g(1) * (-t23 * pkin(1) + t18) - g(2) * t27) * MDP(6) + (-g(1) * (t25 * t28 + t18 + (-pkin(1) + t21) * t23) - g(2) * (-t25 * t21 + t23 * t28 + t27)) * MDP(15) + (MDP(2) - MDP(4) + MDP(14)) * t15 + (-t22 * MDP(12) - t24 * MDP(13) + MDP(3) - MDP(5)) * t16; (-MDP(15) - MDP(6)) * t15; (g(3) * t24 + t15 * t22) * MDP(13) + (MDP(15) * pkin(3) + MDP(12)) * (g(3) * t22 - t15 * t24); -t16 * MDP(15);];
taug = t1;
