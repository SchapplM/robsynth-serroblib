% Calculate potential energy for
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
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4RPRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR1_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:50:22
% EndTime: 2018-11-14 13:50:22
% DurationCPUTime: 0.08s
% Computational Cost: add. (76->41), mult. (53->34), div. (0->0), fcn. (28->8), ass. (0->17)
t41 = qJ(2) + pkin(4);
t34 = qJ(1) + pkin(7);
t29 = sin(t34);
t35 = sin(qJ(1));
t40 = t35 * pkin(1) + pkin(2) * t29;
t30 = cos(t34);
t36 = cos(qJ(1));
t39 = t36 * pkin(1) + pkin(2) * t30;
t38 = pkin(5) + t41;
t37 = -m(3) * pkin(1) - mrSges(2,1);
t31 = qJ(3) + t34;
t28 = qJ(4) + t31;
t27 = cos(t31);
t26 = sin(t31);
t23 = cos(t28);
t22 = sin(t28);
t1 = (-mrSges(1,3) - m(2) * pkin(4) - mrSges(2,3) - m(3) * t41 - mrSges(3,3) - m(4) * t38 - mrSges(4,3) - m(5) * (pkin(6) + t38) - mrSges(5,3)) * g(3) + (-mrSges(1,2) - t36 * mrSges(2,2) - t29 * mrSges(3,1) - t30 * mrSges(3,2) - m(4) * t40 - t26 * mrSges(4,1) - t27 * mrSges(4,2) - m(5) * (pkin(3) * t26 + t40) - t22 * mrSges(5,1) - t23 * mrSges(5,2) + t37 * t35) * g(2) + (-mrSges(1,1) + t35 * mrSges(2,2) - t30 * mrSges(3,1) + t29 * mrSges(3,2) - m(4) * t39 - t27 * mrSges(4,1) + t26 * mrSges(4,2) - m(5) * (pkin(3) * t27 + t39) - t23 * mrSges(5,1) + t22 * mrSges(5,2) + t37 * t36) * g(1);
U  = t1;
