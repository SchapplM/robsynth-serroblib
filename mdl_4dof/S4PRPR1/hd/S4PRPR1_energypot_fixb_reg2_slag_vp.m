% Calculate inertial parameters regressor of potential energy for
% S4PRPR1
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
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PRPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:42:22
% EndTime: 2018-11-14 13:42:22
% DurationCPUTime: 0.04s
% Computational Cost: add. (62->27), mult. (50->32), div. (0->0), fcn. (44->6), ass. (0->17)
t26 = pkin(4) + qJ(1);
t32 = g(3) * t26;
t23 = pkin(6) + qJ(2);
t19 = sin(t23);
t20 = cos(t23);
t25 = cos(pkin(6));
t31 = t25 * pkin(1) + t20 * pkin(2) + t19 * qJ(3);
t24 = sin(pkin(6));
t30 = t24 * pkin(1) + t19 * pkin(2) - t20 * qJ(3);
t29 = -g(1) * t25 - g(2) * t24;
t28 = cos(qJ(4));
t27 = sin(qJ(4));
t15 = -g(1) * t20 - g(2) * t19;
t14 = g(1) * t19 - g(2) * t20;
t13 = t19 * t28 - t20 * t27;
t12 = -t19 * t27 - t20 * t28;
t1 = [0, 0, 0, 0, 0, 0, t29, g(1) * t24 - g(2) * t25, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, t15, t14, -g(3), t29 * pkin(1) - t32, 0, 0, 0, 0, 0, 0, t15, -g(3), -t14, -g(1) * t31 - g(2) * t30 - t32, 0, 0, 0, 0, 0, 0, g(1) * t12 - g(2) * t13, -g(1) * t13 - g(2) * t12, g(3), -g(1) * (t20 * pkin(3) + t31) - g(2) * (t19 * pkin(3) + t30) - g(3) * (-pkin(5) + t26);];
U_reg  = t1;
