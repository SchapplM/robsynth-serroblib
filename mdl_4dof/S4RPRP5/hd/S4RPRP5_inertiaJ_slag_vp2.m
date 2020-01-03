% Calculate joint inertia matrix for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP5_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP5_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:46
% EndTime: 2019-12-31 16:44:47
% DurationCPUTime: 0.17s
% Computational Cost: add. (155->61), mult. (305->78), div. (0->0), fcn. (258->4), ass. (0->23)
t29 = -m(5) * pkin(3) - mrSges(5,1);
t18 = cos(pkin(6));
t16 = t18 ^ 2;
t13 = -t18 * pkin(2) - pkin(1);
t28 = 0.2e1 * t13;
t27 = cos(qJ(3));
t26 = mrSges(5,2) + mrSges(4,3);
t25 = pkin(5) + qJ(2);
t17 = sin(pkin(6));
t24 = t17 ^ 2 + t16;
t11 = t25 * t18;
t19 = sin(qJ(3));
t22 = t25 * t17;
t3 = t19 * t11 + t27 * t22;
t5 = t27 * t11 - t19 * t22;
t23 = t3 ^ 2 + t5 ^ 2;
t21 = -t18 * mrSges(3,1) + t17 * mrSges(3,2);
t10 = t27 * t17 + t19 * t18;
t9 = t19 * t17 - t27 * t18;
t7 = t10 * mrSges(4,2);
t6 = t9 * mrSges(5,1);
t1 = t9 * pkin(3) - t10 * qJ(4) + t13;
t2 = [Ifges(2,3) + t7 * t28 + 0.2e1 * t1 * t6 - 0.2e1 * pkin(1) * t21 + Ifges(3,2) * t16 + (Ifges(3,1) * t17 + 0.2e1 * Ifges(3,4) * t18) * t17 + 0.2e1 * t24 * qJ(2) * mrSges(3,3) + m(4) * (t13 ^ 2 + t23) + m(5) * (t1 ^ 2 + t23) + m(3) * (t24 * qJ(2) ^ 2 + pkin(1) ^ 2) + (mrSges(4,1) * t28 + (Ifges(5,3) + Ifges(4,2)) * t9 - 0.2e1 * t26 * t5) * t9 + (-0.2e1 * t1 * mrSges(5,3) + (Ifges(5,1) + Ifges(4,1)) * t10 + 0.2e1 * (-Ifges(4,4) + Ifges(5,5)) * t9 + 0.2e1 * t26 * t3) * t10; -m(3) * pkin(1) + m(4) * t13 + m(5) * t1 + t9 * mrSges(4,1) - t10 * mrSges(5,3) + t21 + t6 + t7; m(3) + m(4) + m(5); (-qJ(4) * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * t9 + (-pkin(3) * mrSges(5,2) + Ifges(5,4) + Ifges(4,5)) * t10 + (m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t5 + (-mrSges(4,1) + t29) * t3; 0; Ifges(5,2) + Ifges(4,3) + 0.2e1 * pkin(3) * mrSges(5,1) + 0.2e1 * qJ(4) * mrSges(5,3) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2); m(5) * t3 + t10 * mrSges(5,2); 0; t29; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
