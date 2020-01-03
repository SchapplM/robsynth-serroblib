% Calculate joint inertia matrix for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:41
% EndTime: 2019-12-31 18:10:42
% DurationCPUTime: 0.23s
% Computational Cost: add. (156->85), mult. (292->94), div. (0->0), fcn. (147->4), ass. (0->26)
t19 = sin(qJ(3));
t20 = cos(qJ(3));
t37 = t19 ^ 2 + t20 ^ 2;
t36 = 0.2e1 * t37;
t35 = -mrSges(4,2) + mrSges(5,3);
t34 = -m(5) * pkin(3) - mrSges(5,1);
t18 = cos(pkin(7));
t10 = -t18 * pkin(1) - pkin(2);
t11 = t19 * qJ(4);
t26 = t20 * pkin(3) + t11;
t2 = t10 - t26;
t33 = -0.2e1 * t2;
t32 = -2 * mrSges(6,3);
t31 = m(5) + m(6);
t17 = sin(pkin(7));
t9 = t17 * pkin(1) + pkin(6);
t30 = t37 * t9 ^ 2;
t29 = mrSges(5,2) - mrSges(6,3);
t28 = -qJ(5) + t9;
t27 = t20 * mrSges(6,1) + t19 * mrSges(6,2);
t22 = qJ(4) ^ 2;
t21 = -pkin(3) - pkin(4);
t4 = t28 * t20;
t3 = t28 * t19;
t1 = t20 * pkin(4) - t2;
t5 = [Ifges(2,3) + Ifges(3,3) + 0.2e1 * t1 * t27 + m(6) * (t1 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(5) * (t2 ^ 2 + t30) + m(4) * (t10 ^ 2 + t30) + (-0.2e1 * t10 * mrSges(4,1) + mrSges(5,1) * t33 + t4 * t32 + (Ifges(6,2) + Ifges(5,3) + Ifges(4,2)) * t20) * t20 + (0.2e1 * t10 * mrSges(4,2) + mrSges(5,3) * t33 + t3 * t32 + (Ifges(5,1) + Ifges(4,1) + Ifges(6,1)) * t19 + 0.2e1 * (Ifges(4,4) - Ifges(6,4) - Ifges(5,5)) * t20) * t19 + (mrSges(5,2) + mrSges(4,3)) * t9 * t36 + (0.2e1 * t18 * mrSges(3,1) - 0.2e1 * t17 * mrSges(3,2) + m(3) * (t17 ^ 2 + t18 ^ 2) * pkin(1)) * pkin(1); m(6) * (t4 * t19 - t3 * t20); m(3) + (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t36; m(6) * (qJ(4) * t4 + t21 * t3) - t3 * mrSges(6,1) + t4 * mrSges(6,2) + (t29 * qJ(4) + Ifges(4,6) - Ifges(5,6) + Ifges(6,6)) * t20 + (-pkin(3) * mrSges(5,2) - t21 * mrSges(6,3) + Ifges(5,4) + Ifges(4,5) - Ifges(6,5)) * t19 + ((m(5) * qJ(4) + t35) * t20 + (-mrSges(4,1) + t34) * t19) * t9; (mrSges(4,1) + mrSges(5,1)) * t20 + t35 * t19 + m(5) * t26 + m(6) * (-t21 * t20 + t11) + t27; 0.2e1 * pkin(3) * mrSges(5,1) - 0.2e1 * t21 * mrSges(6,1) + Ifges(5,2) + Ifges(4,3) + Ifges(6,3) + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * qJ(4) + m(5) * (pkin(3) ^ 2 + t22) + m(6) * (t21 ^ 2 + t22); m(6) * t3 + (m(5) * t9 + t29) * t19; -t31 * t20; m(6) * t21 - mrSges(6,1) + t34; t31; m(6) * t1 + t27; 0; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
