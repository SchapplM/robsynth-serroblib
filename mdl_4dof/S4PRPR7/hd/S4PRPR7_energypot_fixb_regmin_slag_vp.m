% Calculate minimal parameter regressor of potential energy for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% U_reg [1x14]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:54
% EndTime: 2019-12-31 16:25:54
% DurationCPUTime: 0.09s
% Computational Cost: add. (27->22), mult. (58->34), div. (0->0), fcn. (58->6), ass. (0->14)
t75 = cos(qJ(2));
t80 = g(3) * t75;
t72 = sin(qJ(4));
t73 = sin(qJ(2));
t79 = t72 * t73;
t74 = cos(qJ(4));
t78 = t73 * t74;
t70 = sin(pkin(6));
t71 = cos(pkin(6));
t77 = g(1) * t71 + g(2) * t70;
t76 = pkin(2) * t75 + qJ(3) * t73 + pkin(1);
t69 = g(3) * t73 + t75 * t77;
t68 = t73 * t77 - t80;
t1 = [-g(3) * qJ(1), 0, -t69, t68, t69, -t68, -g(3) * (t73 * pkin(2) - t75 * qJ(3) + qJ(1)) + (g(2) * pkin(4) - g(1) * t76) * t71 + (-g(1) * pkin(4) - g(2) * t76) * t70, 0, 0, 0, 0, 0, -g(1) * (t70 * t74 + t71 * t79) - g(2) * (t70 * t79 - t71 * t74) + t72 * t80, -g(1) * (-t70 * t72 + t71 * t78) - g(2) * (t70 * t78 + t71 * t72) + t74 * t80;];
U_reg = t1;
