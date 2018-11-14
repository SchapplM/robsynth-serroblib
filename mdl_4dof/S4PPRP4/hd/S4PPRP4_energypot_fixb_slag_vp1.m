% Calculate potential energy for
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
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:06
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4PPRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP4_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP4_energypot_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP4_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP4_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:05:52
% EndTime: 2018-11-14 14:05:52
% DurationCPUTime: 0.08s
% Computational Cost: add. (35->32), mult. (30->28), div. (0->0), fcn. (10->2), ass. (0->6)
t13 = -pkin(1) - pkin(4);
t12 = rSges(5,1) + pkin(3);
t11 = pkin(2) + qJ(1);
t10 = cos(qJ(3));
t9 = sin(qJ(3));
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (qJ(1) + rSges(2,3)) + g(2) * rSges(2,2) - g(3) * rSges(2,1)) - m(3) * (g(1) * (qJ(1) + rSges(3,1)) + g(2) * (-qJ(2) - rSges(3,3)) + g(3) * (-pkin(1) + rSges(3,2))) - m(4) * (g(1) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t11) + g(2) * (-t9 * rSges(4,1) - t10 * rSges(4,2) - qJ(2)) + g(3) * (-rSges(4,3) + t13)) - m(5) * (g(1) * (-t9 * rSges(5,2) + t12 * t10 + t11) + g(2) * (-t10 * rSges(5,2) - t12 * t9 - qJ(2)) + g(3) * (-qJ(4) - rSges(5,3) + t13));
U  = t1;
