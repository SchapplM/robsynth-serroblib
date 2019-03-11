% Calculate Gravitation load on the joints for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% MDP [7x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3PRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [3x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S3PRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [7 1]), ...
  'S3PRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [7x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:04:00
% EndTime: 2019-03-08 18:04:00
% DurationCPUTime: 0.02s
% Computational Cost: add. (19->9), mult. (19->13), div. (0->0), fcn. (12->4), ass. (0->7)
t12 = qJ(2) + qJ(3);
t10 = sin(t12);
t11 = cos(t12);
t15 = (g(1) * t10 - g(2) * t11) * MDP(6) + (g(1) * t11 + g(2) * t10) * MDP(7);
t14 = cos(qJ(2));
t13 = sin(qJ(2));
t1 = [-g(2) * MDP(1); (g(1) * t13 - g(2) * t14) * MDP(3) + (g(1) * t14 + g(2) * t13) * MDP(4) + t15; t15;];
taug  = t1;
