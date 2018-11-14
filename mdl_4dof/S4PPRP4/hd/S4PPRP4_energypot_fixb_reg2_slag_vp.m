% Calculate inertial parameters regressor of potential energy for
% S4PPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:06
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PPRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP4_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP4_energypot_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:06:13
% EndTime: 2018-11-14 14:06:13
% DurationCPUTime: 0.03s
% Computational Cost: add. (22->18), mult. (20->14), div. (0->0), fcn. (10->2), ass. (0->9)
t18 = -pkin(4) - pkin(1);
t17 = g(1) * qJ(1);
t16 = pkin(2) + qJ(1);
t15 = cos(qJ(3));
t14 = sin(qJ(3));
t13 = g(2) * qJ(2);
t12 = -g(1) * t15 + g(2) * t14;
t11 = g(1) * t14 + g(2) * t15;
t1 = [0, 0, 0, 0, 0, 0, g(3), -g(2), -g(1), -t17, 0, 0, 0, 0, 0, 0, -g(1), -g(3), g(2), g(3) * pkin(1) + t13 - t17, 0, 0, 0, 0, 0, 0, t12, t11, g(3), -g(1) * t16 - g(3) * t18 + t13, 0, 0, 0, 0, 0, 0, t12, t11, g(3), -g(1) * (t15 * pkin(3) + t16) - g(2) * (-t14 * pkin(3) - qJ(2)) - g(3) * (-qJ(4) + t18);];
U_reg  = t1;
