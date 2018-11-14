% Calculate inertial parameters regressor of potential energy for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4RRPP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:52:30
% EndTime: 2018-11-14 13:52:30
% DurationCPUTime: 0.03s
% Computational Cost: add. (54->24), mult. (42->24), div. (0->0), fcn. (32->4), ass. (0->13)
t30 = pkin(5) + pkin(4);
t34 = g(3) * t30;
t27 = qJ(1) + qJ(2);
t23 = sin(t27);
t24 = cos(t27);
t29 = cos(qJ(1));
t33 = t29 * pkin(1) + t24 * pkin(2) + t23 * qJ(3);
t28 = sin(qJ(1));
t32 = t28 * pkin(1) + t23 * pkin(2) - t24 * qJ(3);
t31 = -g(1) * t29 - g(2) * t28;
t18 = -g(1) * t24 - g(2) * t23;
t17 = g(1) * t23 - g(2) * t24;
t1 = [0, 0, 0, 0, 0, 0, t31, g(1) * t28 - g(2) * t29, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, t18, t17, -g(3), t31 * pkin(1) - t34, 0, 0, 0, 0, 0, 0, t18, -g(3), -t17, -g(1) * t33 - g(2) * t32 - t34, 0, 0, 0, 0, 0, 0, t18, -t17, g(3), -g(1) * (t24 * pkin(3) + t33) - g(2) * (t23 * pkin(3) + t32) - g(3) * (-qJ(4) + t30);];
U_reg  = t1;
