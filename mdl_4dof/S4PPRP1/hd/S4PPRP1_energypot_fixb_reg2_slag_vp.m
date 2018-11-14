% Calculate inertial parameters regressor of potential energy for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PPRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:39:14
% EndTime: 2018-11-14 13:39:14
% DurationCPUTime: 0.04s
% Computational Cost: add. (45->26), mult. (70->30), div. (0->0), fcn. (72->4), ass. (0->17)
t33 = g(3) * (-pkin(4) + qJ(1));
t32 = cos(qJ(3));
t31 = sin(qJ(3));
t30 = g(3) * qJ(1);
t22 = sin(pkin(5));
t23 = cos(pkin(5));
t29 = t23 * pkin(1) + t22 * qJ(2);
t28 = t23 * pkin(2) + t29;
t27 = t22 * pkin(1) - t23 * qJ(2);
t26 = t22 * pkin(2) + t27;
t10 = -t22 * t31 - t23 * t32;
t11 = -t22 * t32 + t23 * t31;
t25 = g(1) * t11 - g(2) * t10;
t13 = -g(1) * t23 - g(2) * t22;
t12 = g(1) * t22 - g(2) * t23;
t9 = g(1) * t10 + g(2) * t11;
t1 = [0, 0, 0, 0, 0, 0, t13, t12, -g(3), -t30, 0, 0, 0, 0, 0, 0, t13, -g(3), -t12, -g(1) * t29 - g(2) * t27 - t30, 0, 0, 0, 0, 0, 0, t9, t25, g(3), -g(1) * t28 - g(2) * t26 - t33, 0, 0, 0, 0, 0, 0, t9, g(3), -t25, -g(1) * (-t10 * pkin(3) + t11 * qJ(4) + t28) - g(2) * (-t11 * pkin(3) - t10 * qJ(4) + t26) - t33;];
U_reg  = t1;
