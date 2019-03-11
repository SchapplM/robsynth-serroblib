% Calculate potential energy for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S3RRP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_energypot_floatb_twist_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S3RRP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_energypot_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_energypot_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRP1_energypot_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:06:50
% EndTime: 2019-03-08 18:06:50
% DurationCPUTime: 0.09s
% Computational Cost: add. (54->28), mult. (44->19), div. (0->0), fcn. (20->4), ass. (0->12)
t14 = -m(3) - m(4);
t15 = pkin(1) * t14 - mrSges(2,1);
t13 = pkin(3) + r_base(3);
t11 = -m(1) - m(2) + t14;
t10 = -m(4) * pkin(2) - mrSges(3,1) - mrSges(4,1);
t9 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t8 = cos(qJ(1));
t7 = sin(qJ(1));
t6 = qJ(1) + qJ(2);
t2 = cos(t6);
t1 = sin(t6);
t3 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t13 - mrSges(2,3) - mrSges(3,3) - mrSges(4,2) + t14 * (pkin(4) + t13)) * g(3) + (-mrSges(2,2) * t8 + t10 * t1 + t11 * r_base(2) + t15 * t7 + t9 * t2 - mrSges(1,2)) * g(2) + (t7 * mrSges(2,2) - t9 * t1 + t10 * t2 + t11 * r_base(1) + t15 * t8 - mrSges(1,1)) * g(1);
U  = t3;
