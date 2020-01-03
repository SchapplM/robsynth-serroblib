% Calculate Gravitation load on the joints for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4RPRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:13
% EndTime: 2019-12-31 16:48:13
% DurationCPUTime: 0.04s
% Computational Cost: add. (66->15), mult. (56->22), div. (0->0), fcn. (40->6), ass. (0->10)
t19 = qJ(1) + pkin(7) + qJ(3);
t17 = sin(t19);
t18 = cos(t19);
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t26 = g(1) * t18 + g(2) * t17;
t27 = t26 * MDP(7) + (t22 * MDP(13) - t20 * MDP(14) + MDP(6)) * (g(1) * t17 - g(2) * t18);
t23 = cos(qJ(1));
t21 = sin(qJ(1));
t1 = [(g(1) * t23 + g(2) * t21) * MDP(3) + t27 + (MDP(4) * pkin(1) + MDP(2)) * (g(1) * t21 - g(2) * t23); -g(3) * MDP(4); t27; (-g(3) * t22 + t26 * t20) * MDP(13) + (g(3) * t20 + t26 * t22) * MDP(14);];
taug = t1;
