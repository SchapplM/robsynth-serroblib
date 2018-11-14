% Calculate inertial parameters regressor of potential energy for
% S4PRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:10
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PRPP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP4_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP4_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:09:29
% EndTime: 2018-11-14 14:09:29
% DurationCPUTime: 0.03s
% Computational Cost: add. (38->23), mult. (30->22), div. (0->0), fcn. (20->4), ass. (0->12)
t30 = g(3) * (-qJ(3) - pkin(4));
t26 = cos(qJ(2));
t29 = t26 * pkin(2) + pkin(1);
t28 = g(1) * qJ(1);
t25 = sin(qJ(2));
t27 = t25 * pkin(2) + qJ(1);
t23 = qJ(2) + pkin(5);
t20 = cos(t23);
t19 = sin(t23);
t18 = g(1) * t20 - g(2) * t19;
t17 = -g(1) * t19 - g(2) * t20;
t1 = [0, 0, 0, 0, 0, 0, -g(2), -g(3), -g(1), -t28, 0, 0, 0, 0, 0, 0, -g(1) * t25 - g(2) * t26, -g(1) * t26 + g(2) * t25, g(3), -g(2) * pkin(1) + g(3) * pkin(4) - t28, 0, 0, 0, 0, 0, 0, t17, -t18, g(3), -g(1) * t27 - g(2) * t29 - t30, 0, 0, 0, 0, 0, 0, t17, g(3), t18, -g(1) * (t19 * pkin(3) - t20 * qJ(4) + t27) - g(2) * (t20 * pkin(3) + t19 * qJ(4) + t29) - t30;];
U_reg  = t1;
