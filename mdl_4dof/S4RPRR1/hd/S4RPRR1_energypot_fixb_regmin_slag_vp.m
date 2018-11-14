% Calculate minimal parameter regressor of potential energy for
% S4RPRR1
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
% U_reg [1x10]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:50:34
% EndTime: 2018-11-14 13:50:34
% DurationCPUTime: 0.02s
% Computational Cost: add. (29->10), mult. (17->14), div. (0->0), fcn. (14->6), ass. (0->10)
t68 = qJ(1) + pkin(7) + qJ(3);
t69 = sin(qJ(1));
t70 = cos(qJ(1));
t71 = -g(1) * t70 - g(2) * t69;
t67 = qJ(4) + t68;
t66 = cos(t68);
t65 = sin(t68);
t64 = cos(t67);
t63 = sin(t67);
t1 = [0, t71, g(1) * t69 - g(2) * t70, -g(3) * (qJ(2) + pkin(4)) + t71 * pkin(1), 0, -g(1) * t66 - g(2) * t65, g(1) * t65 - g(2) * t66, 0, -g(1) * t64 - g(2) * t63, g(1) * t63 - g(2) * t64;];
U_reg  = t1;
