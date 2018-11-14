% Calculate potential energy for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4RPPR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR2_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:24
% EndTime: 2018-11-14 13:47:25
% DurationCPUTime: 0.13s
% Computational Cost: add. (65->39), mult. (79->36), div. (0->0), fcn. (62->6), ass. (0->15)
t41 = m(3) + m(4) + m(5);
t35 = sin(pkin(6));
t40 = m(5) * pkin(3) * t35 - mrSges(2,2) + mrSges(3,3);
t36 = cos(pkin(6));
t39 = -m(4) * pkin(2) - m(5) * (t36 * pkin(3) + pkin(2)) - mrSges(2,1) - mrSges(3,1);
t38 = cos(qJ(1));
t37 = sin(qJ(1));
t34 = pkin(6) + qJ(4);
t30 = cos(t34);
t29 = sin(t34);
t27 = -t38 * t35 + t37 * t36;
t26 = -t37 * t35 - t38 * t36;
t25 = -t38 * t29 + t37 * t30;
t24 = -t37 * t29 - t38 * t30;
t1 = (-mrSges(1,3) - mrSges(2,3) - mrSges(3,2) + m(4) * qJ(3) + mrSges(4,3) - m(5) * (-pkin(5) - qJ(3)) + mrSges(5,3) + (-m(2) - t41) * pkin(4)) * g(3) + (-t27 * mrSges(4,1) - t25 * mrSges(5,1) - t26 * mrSges(4,2) - t24 * mrSges(5,2) - mrSges(1,2) + (qJ(2) * t41 + t40) * t38 + (-t41 * pkin(1) + t39) * t37) * g(2) + (t26 * mrSges(4,1) + t24 * mrSges(5,1) - t27 * mrSges(4,2) - t25 * mrSges(5,2) - t40 * t37 + t39 * t38 - mrSges(1,1) - t41 * (t38 * pkin(1) + t37 * qJ(2))) * g(1);
U  = t1;
