% Calculate Gravitation load on the joints for
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PPRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PPRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:00
% EndTime: 2019-03-08 18:17:00
% DurationCPUTime: 0.02s
% Computational Cost: add. (32->11), mult. (21->14), div. (0->0), fcn. (12->4), ass. (0->8)
t16 = pkin(6) + qJ(3);
t15 = qJ(4) + t16;
t11 = sin(t15);
t12 = cos(t15);
t17 = (g(1) * t11 - g(2) * t12) * MDP(7) + (g(1) * t12 + g(2) * t11) * MDP(8);
t14 = cos(t16);
t13 = sin(t16);
t1 = [(-MDP(1) - MDP(2)) * g(2); -g(3) * MDP(2); (g(1) * t13 - g(2) * t14) * MDP(4) + (g(1) * t14 + g(2) * t13) * MDP(5) + t17; t17;];
taug  = t1;
