% Calculate potential energy for
% S4RRPP5
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP5_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:07
% EndTime: 2019-12-31 17:00:08
% DurationCPUTime: 0.21s
% Computational Cost: add. (64->36), mult. (104->30), div. (0->0), fcn. (79->4), ass. (0->14)
t41 = sin(qJ(2));
t43 = cos(qJ(2));
t61 = pkin(2) * t43 + qJ(3) * t41;
t60 = -m(5) * qJ(4) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t59 = mrSges(3,2) - mrSges(5,2) - mrSges(4,3);
t58 = -m(4) - m(5);
t42 = sin(qJ(1));
t57 = t61 * t42;
t55 = t59 * t41 + t60 * t43 - mrSges(2,1);
t54 = mrSges(4,1) + mrSges(5,1) + mrSges(3,3) - mrSges(2,2);
t39 = t42 * pkin(1);
t44 = cos(qJ(1));
t50 = -t44 * pkin(5) + t39;
t1 = (-mrSges(1,3) - mrSges(2,3) + t58 * (t41 * pkin(2) - t43 * qJ(3) + pkin(4)) + (-m(2) - m(3)) * pkin(4) - t59 * t43 + t60 * t41) * g(3) + (-mrSges(1,2) - m(3) * t50 - m(4) * (t50 + t57) - m(5) * (t39 + t57) + (-m(5) * (-pkin(3) - pkin(5)) + t54) * t44 + t55 * t42) * g(2) + (-mrSges(1,1) + (-m(5) * pkin(3) - t54) * t42 + (-m(3) + t58) * (t44 * pkin(1) + t42 * pkin(5)) + (t58 * t61 + t55) * t44) * g(1);
U = t1;
