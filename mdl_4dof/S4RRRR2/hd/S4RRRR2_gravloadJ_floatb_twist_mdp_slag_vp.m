% Calculate Gravitation load on the joints for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S4RRRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:19
% EndTime: 2019-12-31 17:23:19
% DurationCPUTime: 0.05s
% Computational Cost: add. (98->20), mult. (94->28), div. (0->0), fcn. (74->8), ass. (0->14)
t31 = qJ(3) + qJ(4);
t27 = sin(t31);
t29 = cos(t31);
t32 = qJ(1) + qJ(2);
t28 = sin(t32);
t30 = cos(t32);
t39 = g(1) * t30 + g(2) * t28;
t40 = (-g(3) * t29 + t39 * t27) * MDP(19) + (g(3) * t27 + t39 * t29) * MDP(20);
t33 = sin(qJ(3));
t35 = cos(qJ(3));
t37 = t39 * MDP(6) + (t35 * MDP(12) - t33 * MDP(13) + t29 * MDP(19) - t27 * MDP(20) + MDP(5)) * (g(1) * t28 - g(2) * t30);
t36 = cos(qJ(1));
t34 = sin(qJ(1));
t1 = [(g(1) * t34 - g(2) * t36) * MDP(2) + (g(1) * t36 + g(2) * t34) * MDP(3) + t37; t37; (-g(3) * t35 + t39 * t33) * MDP(12) + (g(3) * t33 + t39 * t35) * MDP(13) + t40; t40;];
taug = t1;
