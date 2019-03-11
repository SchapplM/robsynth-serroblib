% Calculate Gravitation load on the joints for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [3x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S3RRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:06:53
% EndTime: 2019-03-08 18:06:53
% DurationCPUTime: 0.04s
% Computational Cost: add. (55->19), mult. (49->25), div. (0->0), fcn. (32->4), ass. (0->10)
t24 = qJ(1) + qJ(2);
t22 = sin(t24);
t23 = cos(t24);
t29 = t23 * pkin(2) + t22 * qJ(3);
t28 = -t22 * pkin(2) + t23 * qJ(3);
t17 = g(1) * t22 - g(2) * t23;
t27 = (MDP(6) - MDP(8)) * (g(1) * t23 + g(2) * t22) + (MDP(5) + MDP(7)) * t17;
t26 = cos(qJ(1));
t25 = sin(qJ(1));
t1 = [(g(1) * t25 - g(2) * t26) * MDP(2) + (g(1) * t26 + g(2) * t25) * MDP(3) + (-g(1) * (-t25 * pkin(1) + t28) - g(2) * (t26 * pkin(1) + t29)) * MDP(9) + t27; (-g(1) * t28 - g(2) * t29) * MDP(9) + t27; -t17 * MDP(9);];
taug  = t1;
