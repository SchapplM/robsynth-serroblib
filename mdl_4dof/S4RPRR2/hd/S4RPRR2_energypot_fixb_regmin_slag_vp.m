% Calculate minimal parameter regressor of potential energy for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x14]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:12
% EndTime: 2019-12-31 16:48:12
% DurationCPUTime: 0.02s
% Computational Cost: add. (27->10), mult. (23->14), div. (0->0), fcn. (20->6), ass. (0->10)
t77 = qJ(1) + pkin(7) + qJ(3);
t75 = sin(t77);
t76 = cos(t77);
t83 = g(1) * t76 + g(2) * t75;
t79 = sin(qJ(1));
t81 = cos(qJ(1));
t82 = -g(1) * t81 - g(2) * t79;
t80 = cos(qJ(4));
t78 = sin(qJ(4));
t1 = [0, t82, g(1) * t79 - g(2) * t81, -g(3) * (qJ(2) + pkin(4)) + t82 * pkin(1), 0, -t83, g(1) * t75 - g(2) * t76, 0, 0, 0, 0, 0, -g(3) * t78 - t83 * t80, -g(3) * t80 + t83 * t78;];
U_reg = t1;
