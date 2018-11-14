% Calculate minimal parameter regressor of potential energy for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% U_reg [1x10]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:55
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4RRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:54:30
% EndTime: 2018-11-14 13:54:30
% DurationCPUTime: 0.02s
% Computational Cost: add. (33->17), mult. (21->21), div. (0->0), fcn. (18->6), ass. (0->9)
t59 = qJ(1) + qJ(2);
t61 = cos(qJ(1));
t60 = sin(qJ(1));
t58 = qJ(3) + t59;
t57 = cos(t59);
t56 = sin(t59);
t55 = cos(t58);
t54 = sin(t58);
t1 = [0, -g(1) * t61 - g(2) * t60, g(1) * t60 - g(2) * t61, 0, -g(1) * t57 - g(2) * t56, g(1) * t56 - g(2) * t57, 0, -g(1) * t55 - g(2) * t54, g(1) * t54 - g(2) * t55, -g(1) * (t61 * pkin(1) + pkin(2) * t57 + pkin(3) * t55) - g(2) * (t60 * pkin(1) + pkin(2) * t56 + pkin(3) * t54) - g(3) * (qJ(4) + pkin(6) + pkin(5) + pkin(4));];
U_reg  = t1;
