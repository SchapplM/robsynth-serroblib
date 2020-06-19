% Calculate Gravitation load on the joints for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S2RR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S2RR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S2RR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:27
% EndTime: 2020-06-19 09:14:27
% DurationCPUTime: 0.04s
% Computational Cost: add. (18->8), mult. (18->12), div. (0->0), fcn. (12->4), ass. (0->7)
t12 = qJ(1) + qJ(2);
t10 = sin(t12);
t11 = cos(t12);
t15 = (g(1) * t10 - g(2) * t11) * MDP(5) + (g(1) * t11 + g(2) * t10) * MDP(6);
t14 = cos(qJ(1));
t13 = sin(qJ(1));
t1 = [(g(1) * t13 - g(2) * t14) * MDP(2) + (g(1) * t14 + g(2) * t13) * MDP(3) + t15; t15;];
taug = t1;
