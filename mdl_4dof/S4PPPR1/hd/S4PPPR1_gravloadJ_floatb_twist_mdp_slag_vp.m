% Calculate Gravitation load on the joints for
% S4PPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta1]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PPPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S4PPPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:09:04
% EndTime: 2019-03-08 18:09:04
% DurationCPUTime: 0.02s
% Computational Cost: add. (13->9), mult. (26->17), div. (0->0), fcn. (22->4), ass. (0->8)
t15 = MDP(2) + MDP(3);
t14 = cos(qJ(4));
t13 = sin(qJ(4));
t12 = cos(pkin(5));
t11 = sin(pkin(5));
t9 = t11 * t14 + t12 * t13;
t8 = -t11 * t13 + t12 * t14;
t1 = [(-MDP(1) - t15) * g(3); t15 * (-g(1) * t11 + g(2) * t12); (-g(1) * t12 - g(2) * t11) * MDP(3); (-g(1) * t8 - g(2) * t9) * MDP(5) + (g(1) * t9 - g(2) * t8) * MDP(6);];
taug  = t1;
