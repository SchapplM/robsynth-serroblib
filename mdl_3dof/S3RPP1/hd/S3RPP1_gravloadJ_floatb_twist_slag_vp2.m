% Calculate Gravitation load on the joints for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [3x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S3RPP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:13:43
% EndTime: 2018-11-14 10:13:43
% DurationCPUTime: 0.09s
% Computational Cost: add. (25->15), mult. (44->19), div. (0->0), fcn. (28->2), ass. (0->6)
t10 = m(3) + m(4);
t9 = mrSges(2,1) + mrSges(4,3) - mrSges(3,2);
t8 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2);
t5 = cos(qJ(1));
t4 = sin(qJ(1));
t1 = [(-t10 * (t5 * pkin(1) + t4 * qJ(2)) + (-m(4) * qJ(3) - t9) * t5 + t8 * t4) * g(2) + ((m(3) * pkin(1) - m(4) * (-pkin(1) - qJ(3)) + t9) * t4 + (-t10 * qJ(2) + t8) * t5) * g(1), t10 * (-g(1) * t4 + g(2) * t5) (-g(1) * t5 - g(2) * t4) * m(4)];
taug  = t1(:);
