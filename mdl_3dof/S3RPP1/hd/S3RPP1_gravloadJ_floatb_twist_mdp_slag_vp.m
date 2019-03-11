% Calculate Gravitation load on the joints for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [3x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S3RPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:01
% EndTime: 2019-03-08 18:05:01
% DurationCPUTime: 0.04s
% Computational Cost: add. (25->18), mult. (43->20), div. (0->0), fcn. (28->2), ass. (0->7)
t16 = sin(qJ(1));
t17 = cos(qJ(1));
t18 = t17 * pkin(1) + t16 * qJ(2);
t13 = t17 * qJ(2);
t11 = g(1) * t17 + g(2) * t16;
t10 = g(1) * t16 - g(2) * t17;
t1 = [(-g(1) * (-t16 * pkin(1) + t13) - g(2) * t18) * MDP(6) + (-g(1) * (t13 + (-pkin(1) - qJ(3)) * t16) - g(2) * (t17 * qJ(3) + t18)) * MDP(9) + (MDP(3) - MDP(5) - MDP(7)) * t11 + (MDP(2) - MDP(4) + MDP(8)) * t10; (-MDP(6) - MDP(9)) * t10; -t11 * MDP(9);];
taug  = t1;
