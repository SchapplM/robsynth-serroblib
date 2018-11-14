% Calculate potential energy for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4RRPP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP1_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:51:28
% EndTime: 2018-11-14 13:51:28
% DurationCPUTime: 0.11s
% Computational Cost: add. (77->35), mult. (55->27), div. (0->0), fcn. (30->6), ass. (0->14)
t41 = -m(4) - m(5);
t40 = pkin(5) + pkin(4);
t32 = qJ(1) + qJ(2);
t37 = -m(3) * pkin(1) - mrSges(2,1);
t36 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t35 = m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3);
t34 = cos(qJ(1));
t33 = sin(qJ(1));
t28 = cos(t32);
t27 = sin(t32);
t26 = pkin(6) + t32;
t23 = cos(t26);
t22 = sin(t26);
t1 = (-mrSges(1,3) - m(2) * pkin(4) - mrSges(2,3) - m(3) * t40 - mrSges(3,3) - mrSges(4,3) - mrSges(5,2) + t41 * (qJ(3) + t40)) * g(3) + (-t27 * mrSges(3,1) - t34 * mrSges(2,2) - t28 * mrSges(3,2) + t36 * t22 + t35 * t23 + t37 * t33 - mrSges(1,2) + t41 * (t33 * pkin(1) + pkin(2) * t27)) * g(2) + (-t28 * mrSges(3,1) + t33 * mrSges(2,2) + t27 * mrSges(3,2) - t35 * t22 + t36 * t23 + t37 * t34 - mrSges(1,1) + t41 * (t34 * pkin(1) + pkin(2) * t28)) * g(1);
U  = t1;
