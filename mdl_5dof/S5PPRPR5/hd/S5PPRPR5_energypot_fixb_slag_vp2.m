% Calculate potential energy for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRPR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR5_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:22
% EndTime: 2019-12-31 17:33:22
% DurationCPUTime: 0.23s
% Computational Cost: add. (95->46), mult. (139->38), div. (0->0), fcn. (132->6), ass. (0->18)
t66 = m(5) + m(6);
t64 = -mrSges(2,1) - mrSges(3,1);
t63 = mrSges(2,2) - mrSges(3,3);
t62 = m(6) * pkin(6) + mrSges(4,1) - mrSges(5,2) + mrSges(6,3);
t49 = sin(qJ(5));
t50 = cos(qJ(5));
t61 = t49 * mrSges(6,1) + t50 * mrSges(6,2) + t66 * qJ(4) - mrSges(4,2) + mrSges(5,3);
t60 = cos(qJ(3));
t59 = sin(qJ(3));
t47 = cos(pkin(7));
t57 = sin(pkin(7));
t58 = t47 * pkin(1) + t57 * qJ(2);
t56 = t47 * pkin(2) + t58;
t55 = t57 * pkin(1) - t47 * qJ(2);
t53 = t57 * pkin(2) + t55;
t38 = t47 * t59 - t57 * t60;
t37 = -t47 * t60 - t57 * t59;
t1 = (m(6) * pkin(4) + t50 * mrSges(6,1) - t49 * mrSges(6,2) + mrSges(5,1) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + (-m(4) - t66) * (-pkin(5) + qJ(1)) + (-m(2) - m(3)) * qJ(1)) * g(3) + (-m(3) * t55 - m(4) * t53 - mrSges(1,2) + t64 * t57 - t66 * (-t38 * pkin(3) + t53) - t63 * t47 + t62 * t38 + t61 * t37) * g(2) + (-m(3) * t58 - m(4) * t56 - mrSges(1,1) + t63 * t57 - t66 * (-t37 * pkin(3) + t56) + t64 * t47 - t61 * t38 + t62 * t37) * g(1);
U = t1;
