% Calculate Gravitation load on the joints for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:30
% EndTime: 2019-12-31 17:33:31
% DurationCPUTime: 0.07s
% Computational Cost: add. (48->20), mult. (97->30), div. (0->0), fcn. (102->6), ass. (0->11)
t31 = cos(qJ(3));
t30 = sin(qJ(3));
t29 = cos(pkin(7));
t28 = sin(pkin(7));
t27 = MDP(2) + MDP(8);
t17 = -t28 * t30 - t29 * t31;
t18 = -t28 * t31 + t29 * t30;
t26 = g(1) * t18 - g(2) * t17;
t24 = cos(qJ(5));
t23 = sin(qJ(5));
t1 = [(-MDP(1) - t27) * g(3); t27 * (-g(1) * t28 + g(2) * t29); (-g(1) * (-t18 * pkin(3) - t17 * qJ(4)) - g(2) * (t17 * pkin(3) - t18 * qJ(4))) * MDP(8) + (t23 * MDP(14) + t24 * MDP(15) - MDP(5) + MDP(7)) * (g(1) * t17 + g(2) * t18) + (MDP(4) - MDP(6)) * t26; -t26 * MDP(8); (-g(3) * t23 - t26 * t24) * MDP(14) + (-g(3) * t24 + t26 * t23) * MDP(15);];
taug = t1;
