% Calculate minimal parameter regressor of potential energy for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% U_reg [1x12]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PPRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:39
% EndTime: 2019-12-31 16:18:39
% DurationCPUTime: 0.06s
% Computational Cost: add. (29->17), mult. (38->27), div. (0->0), fcn. (38->6), ass. (0->15)
t71 = pkin(7) + qJ(3);
t69 = sin(t71);
t82 = g(3) * t69;
t81 = g(3) * qJ(1);
t72 = sin(pkin(6));
t74 = sin(qJ(4));
t80 = t72 * t74;
t75 = cos(qJ(4));
t79 = t72 * t75;
t73 = cos(pkin(6));
t78 = t73 * t74;
t77 = t73 * t75;
t76 = g(1) * t73 + g(2) * t72;
t70 = cos(t71);
t1 = [-t81, -g(1) * (pkin(1) * t73 + qJ(2) * t72) - g(2) * (pkin(1) * t72 - qJ(2) * t73) - t81, 0, -t76 * t70 - t82, -g(3) * t70 + t76 * t69, 0, 0, 0, 0, 0, -g(1) * (t70 * t77 + t80) - g(2) * (t70 * t79 - t78) - t75 * t82, -g(1) * (-t70 * t78 + t79) - g(2) * (-t70 * t80 - t77) + t74 * t82;];
U_reg = t1;
