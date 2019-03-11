% Calculate Gravitation load on the joints for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [3x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S3RPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:57
% EndTime: 2019-03-08 18:05:57
% DurationCPUTime: 0.04s
% Computational Cost: add. (27->16), mult. (50->24), div. (0->0), fcn. (46->4), ass. (0->9)
t21 = sin(qJ(1));
t22 = cos(qJ(1));
t25 = sin(qJ(3));
t26 = cos(qJ(3));
t14 = -t21 * t25 - t22 * t26;
t15 = -t21 * t26 + t22 * t25;
t27 = (g(1) * t15 - g(2) * t14) * MDP(8) - (g(1) * t14 + g(2) * t15) * MDP(9);
t16 = g(1) * t21 - g(2) * t22;
t1 = [(-g(1) * (-t21 * pkin(1) + t22 * qJ(2)) - g(2) * (t22 * pkin(1) + t21 * qJ(2))) * MDP(6) + (MDP(3) - MDP(5)) * (g(1) * t22 + g(2) * t21) + (MDP(2) + MDP(4)) * t16 - t27; -t16 * MDP(6); t27;];
taug  = t1;
