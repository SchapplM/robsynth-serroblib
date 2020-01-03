% Calculate potential energy for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP5_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:04
% EndTime: 2019-12-31 18:16:05
% DurationCPUTime: 0.28s
% Computational Cost: add. (83->55), mult. (124->43), div. (0->0), fcn. (93->4), ass. (0->17)
t63 = -mrSges(4,1) - mrSges(5,1);
t62 = -mrSges(4,2) + mrSges(6,2);
t61 = -m(5) - m(6);
t43 = sin(qJ(3));
t45 = cos(qJ(3));
t60 = t63 * t43 + t62 * t45 + mrSges(2,2) - mrSges(3,3);
t59 = m(6) * qJ(5) - mrSges(2,1) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t58 = pkin(2) + pkin(5);
t57 = pkin(3) * t43;
t44 = sin(qJ(1));
t39 = t44 * pkin(1);
t54 = t44 * pkin(6) + t39;
t46 = cos(qJ(1));
t53 = t46 * pkin(1) + t44 * qJ(2);
t52 = qJ(4) * t45;
t50 = t46 * pkin(6) + t53;
t1 = (-m(4) * t58 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t61 * (t45 * pkin(3) + t43 * qJ(4) + t58) + (-m(2) - m(3)) * pkin(5) + (-m(6) * pkin(4) - mrSges(6,1) + t63) * t45 + (-mrSges(5,3) - t62) * t43) * g(3) + (-m(3) * t39 - m(4) * t54 - mrSges(1,2) + t61 * (t46 * t52 + t54) + (m(5) * t57 - t45 * mrSges(5,3) - (m(6) * (-pkin(3) - pkin(4)) - mrSges(6,1)) * t43 + (m(3) + m(4) - t61) * qJ(2) - t60) * t46 + t59 * t44) * g(2) + (-m(3) * t53 - m(4) * t50 - mrSges(1,1) + t61 * (t44 * t57 + t50) + t59 * t46 + (-(-m(5) * qJ(4) - mrSges(5,3)) * t45 - m(6) * (pkin(4) * t43 - t52) - t43 * mrSges(6,1) + t60) * t44) * g(1);
U = t1;
