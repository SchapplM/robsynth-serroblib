% Calculate potential energy for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-28 15:34
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR2_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-28 15:33:55
% EndTime: 2019-05-28 15:33:55
% DurationCPUTime: 0.11s
% Computational Cost: add. (75->36), mult. (65->31), div. (0->0), fcn. (44->6), ass. (0->14)
t44 = m(4) + m(5);
t43 = mrSges(3,2) - mrSges(4,3);
t40 = -m(3) * pkin(1) - mrSges(2,1);
t39 = -m(5) * pkin(3) - mrSges(3,1) - mrSges(4,1);
t37 = cos(qJ(1));
t36 = cos(qJ(4));
t35 = sin(qJ(1));
t34 = sin(qJ(4));
t33 = qJ(1) + qJ(2);
t30 = cos(t33);
t29 = sin(t33);
t25 = t29 * t36 - t30 * t34;
t24 = -t29 * t34 - t30 * t36;
t1 = (-m(2) * pkin(4) - mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + mrSges(5,3) + (-m(3) - t44) * (pkin(5) + pkin(4))) * g(3) + (-mrSges(1,2) - t37 * mrSges(2,2) - t25 * mrSges(5,1) - t24 * mrSges(5,2) + t40 * t35 + (qJ(3) * t44 - t43) * t30 + t39 * t29 - t44 * (t35 * pkin(1) + t29 * pkin(2))) * g(2) + (t24 * mrSges(5,1) + t35 * mrSges(2,2) - t25 * mrSges(5,2) + t43 * t29 + t39 * t30 + t40 * t37 - mrSges(1,1) - t44 * (t37 * pkin(1) + t30 * pkin(2) + t29 * qJ(3))) * g(1);
U  = t1;
