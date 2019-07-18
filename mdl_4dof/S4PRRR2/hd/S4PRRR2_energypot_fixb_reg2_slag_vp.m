% Calculate inertial parameters regressor of potential energy for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_energypot_fixb_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:23
% EndTime: 2019-07-18 13:27:23
% DurationCPUTime: 0.03s
% Computational Cost: add. (27->14), mult. (26->20), div. (0->0), fcn. (18->6), ass. (0->11)
t21 = qJ(2) + qJ(3);
t23 = sin(qJ(2));
t24 = cos(qJ(2));
t25 = g(1) * t23 - g(3) * t24;
t22 = g(2) * qJ(1);
t20 = qJ(4) + t21;
t19 = cos(t21);
t18 = sin(t21);
t17 = cos(t20);
t16 = sin(t20);
t1 = [0, 0, 0, 0, 0, 0, -g(3), g(1), g(2), t22, 0, 0, 0, 0, 0, 0, t25, g(1) * t24 + g(3) * t23, g(2), t22, 0, 0, 0, 0, 0, 0, g(1) * t18 - g(3) * t19, g(1) * t19 + g(3) * t18, g(2), t25 * pkin(1) + t22, 0, 0, 0, 0, 0, 0, g(1) * t16 - g(3) * t17, g(1) * t17 + g(3) * t16, g(2), -g(1) * (-t23 * pkin(1) - pkin(2) * t18) + t22 - g(3) * (t24 * pkin(1) + pkin(2) * t19);];
U_reg  = t1;
