% Calculate Gravitation load on the joints for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:52:20
% EndTime: 2021-01-15 20:52:21
% DurationCPUTime: 0.18s
% Computational Cost: add. (148->46), mult. (277->60), div. (0->0), fcn. (264->6), ass. (0->25)
t72 = sin(qJ(2));
t74 = cos(qJ(2));
t79 = sin(qJ(4));
t80 = cos(qJ(4));
t56 = -t72 * t79 - t74 * t80;
t57 = t72 * t80 - t74 * t79;
t73 = sin(qJ(1));
t75 = cos(qJ(1));
t76 = g(1) * t75 + g(2) * t73;
t85 = g(3) * t56 + t76 * t57;
t87 = MDP(21) + MDP(23);
t88 = MDP(20) + MDP(22);
t92 = -t88 * t85 - t87 * (-g(3) * t57 + t76 * t56);
t90 = MDP(9) + MDP(11);
t89 = MDP(10) - MDP(13);
t78 = t74 * pkin(2) + t72 * qJ(3);
t63 = t80 * pkin(4) + pkin(2) + pkin(3);
t64 = t79 * pkin(4) + qJ(3);
t77 = t63 * t74 + t64 * t72;
t61 = g(1) * t73 - g(2) * t75;
t71 = qJ(5) - pkin(6) + pkin(7);
t58 = pkin(1) + t78;
t51 = -g(3) * t74 + t76 * t72;
t49 = pkin(1) + t77;
t1 = [(-g(1) * (t75 * pkin(6) - t58 * t73) - g(2) * (t73 * pkin(6) + t58 * t75)) * MDP(14) + (-g(1) * (-t49 * t73 - t71 * t75) - g(2) * (t49 * t75 - t71 * t73)) * MDP(25) + (MDP(3) - MDP(12) + MDP(24)) * t76 + (-t88 * t56 + t87 * t57 - t89 * t72 + t90 * t74 + MDP(2)) * t61; (-g(3) * t78 - t76 * (-t72 * pkin(2) + t74 * qJ(3))) * MDP(14) + (-g(3) * t77 - t76 * (-t63 * t72 + t64 * t74)) * MDP(25) + t89 * (g(3) * t72 + t76 * t74) + t90 * t51 - t92; (-MDP(14) - MDP(25)) * t51; -t85 * MDP(25) * pkin(4) + t92; t61 * MDP(25);];
taug = t1;
