% Calculate minimal parameter regressor of potential energy for
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
% U_reg [1x10]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:42:22
% EndTime: 2018-11-14 13:42:22
% DurationCPUTime: 0.02s
% Computational Cost: add. (38->16), mult. (30->22), div. (0->0), fcn. (30->6), ass. (0->10)
t50 = cos(qJ(4));
t49 = sin(qJ(4));
t48 = pkin(6) + qJ(2);
t47 = cos(t48);
t46 = sin(t48);
t45 = -g(1) * t47 - g(2) * t46;
t44 = g(1) * t46 - g(2) * t47;
t43 = t46 * t50 - t47 * t49;
t42 = -t46 * t49 - t47 * t50;
t1 = [-g(3) * qJ(1), 0, t45, t44, t45, -t44, -g(1) * (t47 * pkin(2) + t46 * qJ(3) + cos(pkin(6)) * pkin(1)) - g(2) * (t46 * pkin(2) - t47 * qJ(3) + sin(pkin(6)) * pkin(1)) - g(3) * (pkin(4) + qJ(1)) 0, g(1) * t42 - g(2) * t43, -g(1) * t43 - g(2) * t42;];
U_reg  = t1;
