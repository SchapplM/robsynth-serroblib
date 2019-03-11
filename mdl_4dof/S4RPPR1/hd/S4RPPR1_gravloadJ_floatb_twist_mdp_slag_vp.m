% Calculate Gravitation load on the joints for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:27:31
% EndTime: 2019-03-08 18:27:31
% DurationCPUTime: 0.08s
% Computational Cost: add. (58->23), mult. (59->34), div. (0->0), fcn. (50->6), ass. (0->12)
t25 = qJ(1) + pkin(6);
t23 = sin(t25);
t24 = cos(t25);
t31 = sin(qJ(4));
t32 = cos(qJ(4));
t17 = -t23 * t31 - t24 * t32;
t18 = -t23 * t32 + t24 * t31;
t33 = -(g(1) * t17 + g(2) * t18) * MDP(10) + (g(1) * t18 - g(2) * t17) * MDP(9);
t27 = cos(qJ(1));
t26 = sin(qJ(1));
t19 = g(1) * t23 - g(2) * t24;
t1 = [(g(1) * t27 + g(2) * t26) * MDP(3) + t19 * MDP(5) + (-g(1) * t24 - g(2) * t23) * MDP(6) + (-g(1) * (-t26 * pkin(1) - t23 * pkin(2) + t24 * qJ(3)) - g(2) * (t27 * pkin(1) + t24 * pkin(2) + t23 * qJ(3))) * MDP(7) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t26 - g(2) * t27) - t33; (-MDP(4) - MDP(7)) * g(3); -t19 * MDP(7); t33;];
taug  = t1;
