% Calculate minimal parameter regressor of potential energy for
% S4PRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
% 
% Output:
% U_reg [1x10]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:42
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PRPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:41:10
% EndTime: 2018-11-14 13:41:10
% DurationCPUTime: 0.02s
% Computational Cost: add. (46->20), mult. (33->19), div. (0->0), fcn. (26->4), ass. (0->9)
t44 = pkin(4) + qJ(1);
t41 = pkin(5) + qJ(2);
t37 = sin(t41);
t38 = cos(t41);
t43 = t38 * pkin(2) + t37 * qJ(3) + cos(pkin(5)) * pkin(1);
t42 = t37 * pkin(2) - t38 * qJ(3) + sin(pkin(5)) * pkin(1);
t32 = g(1) * t38 + g(2) * t37;
t31 = g(1) * t37 - g(2) * t38;
t1 = [-g(3) * qJ(1), 0, -t32, t31, t32, -t31, -g(1) * t43 - g(2) * t42 - g(3) * t44, -t31, -t32, -g(1) * (t38 * qJ(4) + t43) - g(2) * (t37 * qJ(4) + t42) - g(3) * (pkin(3) + t44);];
U_reg  = t1;
