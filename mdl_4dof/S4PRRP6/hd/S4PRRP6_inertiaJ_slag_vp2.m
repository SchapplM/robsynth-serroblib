% Calculate joint inertia matrix for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP6_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP6_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:13
% EndTime: 2019-12-31 16:30:14
% DurationCPUTime: 0.15s
% Computational Cost: add. (74->49), mult. (196->61), div. (0->0), fcn. (96->4), ass. (0->17)
t15 = cos(qJ(3));
t11 = t15 ^ 2;
t13 = sin(qJ(3));
t20 = t13 ^ 2 + t11;
t25 = -m(5) * pkin(3) - mrSges(5,1);
t24 = (mrSges(5,2) + mrSges(4,3)) * t20;
t14 = sin(qJ(2));
t23 = t20 * pkin(5) * t14;
t22 = t20 * pkin(5) ^ 2;
t18 = (m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t15 + (-mrSges(4,1) + t25) * t13;
t16 = cos(qJ(2));
t12 = t16 ^ 2;
t10 = t14 ^ 2;
t4 = -t15 * mrSges(4,1) + t13 * mrSges(4,2);
t3 = -t15 * mrSges(5,1) - t13 * mrSges(5,3);
t2 = -t15 * pkin(3) - t13 * qJ(4) - pkin(2);
t1 = [m(2) + m(3) * (t10 + t12) + 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (t20 * t10 + t12); (mrSges(3,1) - t3 - t4) * t16 + m(4) * (t16 * pkin(2) + t23) + m(5) * (-t2 * t16 + t23) + (-mrSges(3,2) + t24) * t14; -0.2e1 * pkin(2) * t4 + 0.2e1 * t2 * t3 + Ifges(3,3) + m(5) * (t2 ^ 2 + t22) + m(4) * (pkin(2) ^ 2 + t22) + (Ifges(4,2) + Ifges(5,3)) * t11 + ((Ifges(5,1) + Ifges(4,1)) * t13 + 0.2e1 * (Ifges(4,4) - Ifges(5,5)) * t15) * t13 + 0.2e1 * pkin(5) * t24; t18 * t14; (qJ(4) * mrSges(5,2) + Ifges(4,6) - Ifges(5,6)) * t15 + (-pkin(3) * mrSges(5,2) + Ifges(5,4) + Ifges(4,5)) * t13 + t18 * pkin(5); Ifges(5,2) + Ifges(4,3) + 0.2e1 * pkin(3) * mrSges(5,1) + 0.2e1 * qJ(4) * mrSges(5,3) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2); m(5) * t13 * t14; (m(5) * pkin(5) + mrSges(5,2)) * t13; t25; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
