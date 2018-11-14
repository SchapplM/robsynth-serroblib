% Calculate potential energy for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2018-11-14 13:50
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4RPRP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP2_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:49:21
% EndTime: 2018-11-14 13:49:21
% DurationCPUTime: 0.12s
% Computational Cost: add. (57->34), mult. (79->28), div. (0->0), fcn. (62->4), ass. (0->12)
t37 = m(3) + m(4) + m(5);
t40 = mrSges(4,1) + mrSges(5,1);
t39 = -mrSges(4,2) - mrSges(5,2);
t31 = sin(qJ(3));
t36 = m(5) * pkin(3) * t31 - mrSges(2,2) + mrSges(3,3);
t33 = cos(qJ(3));
t35 = -m(4) * pkin(2) - m(5) * (t33 * pkin(3) + pkin(2)) - mrSges(2,1) - mrSges(3,1);
t34 = cos(qJ(1));
t32 = sin(qJ(1));
t24 = -t34 * t31 + t32 * t33;
t23 = -t32 * t31 - t34 * t33;
t1 = (-mrSges(1,3) - mrSges(2,3) - mrSges(3,2) + m(4) * pkin(5) + mrSges(4,3) - m(5) * (-qJ(4) - pkin(5)) + mrSges(5,3) + (-m(2) - t37) * pkin(4)) * g(3) + (-mrSges(1,2) - t40 * t24 + t39 * t23 + (t37 * qJ(2) + t36) * t34 + (-t37 * pkin(1) + t35) * t32) * g(2) + (t40 * t23 + t39 * t24 - t36 * t32 + t35 * t34 - mrSges(1,1) - t37 * (t34 * pkin(1) + t32 * qJ(2))) * g(1);
U  = t1;
