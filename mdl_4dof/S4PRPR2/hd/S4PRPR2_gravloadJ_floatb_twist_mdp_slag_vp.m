% Calculate Gravitation load on the joints for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:49
% EndTime: 2019-03-08 18:21:49
% DurationCPUTime: 0.02s
% Computational Cost: add. (30->11), mult. (26->15), div. (0->0), fcn. (14->4), ass. (0->7)
t13 = qJ(2) + pkin(6) + qJ(4);
t11 = sin(t13);
t12 = cos(t13);
t17 = (g(1) * t11 - g(2) * t12) * MDP(7) + (g(1) * t12 + g(2) * t11) * MDP(8);
t15 = cos(qJ(2));
t14 = sin(qJ(2));
t1 = [(-MDP(1) - MDP(5)) * g(2); (g(1) * t15 + g(2) * t14) * MDP(4) + t17 + (MDP(5) * pkin(2) + MDP(3)) * (g(1) * t14 - g(2) * t15); -g(3) * MDP(5); t17;];
taug  = t1;
