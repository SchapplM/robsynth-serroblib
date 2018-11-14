% Calculate minimal parameter regressor of potential energy for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
% 
% Output:
% U_reg [1x8]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:13
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:12:29
% EndTime: 2018-11-14 14:12:29
% DurationCPUTime: 0.02s
% Computational Cost: add. (19->14), mult. (16->16), div. (0->0), fcn. (12->4), ass. (0->6)
t46 = cos(qJ(2));
t45 = sin(qJ(2));
t44 = qJ(2) + qJ(3);
t43 = cos(t44);
t42 = sin(t44);
t1 = [-g(1) * qJ(1), 0, -g(1) * t45 - g(2) * t46, -g(1) * t46 + g(2) * t45, 0, -g(1) * t42 - g(2) * t43, -g(1) * t43 + g(2) * t42, -g(1) * (t45 * pkin(2) + pkin(3) * t42 + qJ(1)) - g(2) * (t46 * pkin(2) + pkin(3) * t43 + pkin(1)) - g(3) * (-qJ(4) - pkin(5) - pkin(4));];
U_reg  = t1;
