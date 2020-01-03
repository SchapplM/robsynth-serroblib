% Calculate minimal parameter regressor of potential energy for
% S4PRPR4
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
% U_reg [1x14]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:59
% EndTime: 2019-12-31 16:21:59
% DurationCPUTime: 0.03s
% Computational Cost: add. (32->15), mult. (28->18), div. (0->0), fcn. (24->6), ass. (0->8)
t46 = pkin(6) + qJ(2);
t44 = sin(t46);
t45 = cos(t46);
t42 = g(1) * t44 - g(2) * t45;
t48 = cos(qJ(4));
t47 = sin(qJ(4));
t43 = g(1) * t45 + g(2) * t44;
t1 = [-g(3) * qJ(1), 0, -t43, t42, t43, -t42, -g(1) * (t45 * pkin(2) + t44 * qJ(3) + cos(pkin(6)) * pkin(1)) - g(2) * (t44 * pkin(2) - t45 * qJ(3) + sin(pkin(6)) * pkin(1)) - g(3) * (pkin(4) + qJ(1)), 0, 0, 0, 0, 0, -g(3) * t48 - t42 * t47, g(3) * t47 - t42 * t48;];
U_reg = t1;
