% Calculate Gravitation load on the joints for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:29:53
% EndTime: 2019-03-08 18:29:53
% DurationCPUTime: 0.04s
% Computational Cost: add. (88->24), mult. (58->29), div. (0->0), fcn. (36->6), ass. (0->11)
t27 = qJ(1) + pkin(6);
t26 = qJ(3) + t27;
t24 = sin(t26);
t25 = cos(t26);
t33 = t25 * pkin(3) + t24 * qJ(4);
t32 = -t24 * pkin(3) + t25 * qJ(4);
t19 = g(1) * t24 - g(2) * t25;
t31 = (MDP(7) - MDP(9)) * (g(1) * t25 + g(2) * t24) + (MDP(6) + MDP(8)) * t19;
t29 = cos(qJ(1));
t28 = sin(qJ(1));
t1 = [(g(1) * t29 + g(2) * t28) * MDP(3) + (-g(1) * (-pkin(2) * sin(t27) - t28 * pkin(1) + t32) - g(2) * (pkin(2) * cos(t27) + t29 * pkin(1) + t33)) * MDP(10) + t31 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t28 - g(2) * t29); (-MDP(10) - MDP(4)) * g(3); (-g(1) * t32 - g(2) * t33) * MDP(10) + t31; -t19 * MDP(10);];
taug  = t1;
