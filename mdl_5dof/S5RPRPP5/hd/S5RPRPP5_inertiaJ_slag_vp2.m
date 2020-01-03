% Calculate joint inertia matrix for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:08
% EndTime: 2019-12-31 18:16:09
% DurationCPUTime: 0.24s
% Computational Cost: add. (165->91), mult. (281->94), div. (0->0), fcn. (120->2), ass. (0->28)
t14 = sin(qJ(3));
t15 = cos(qJ(3));
t25 = t14 ^ 2 + t15 ^ 2;
t37 = 2 * mrSges(6,3);
t17 = -pkin(1) - pkin(6);
t36 = t17 * t25;
t9 = qJ(4) * t14;
t35 = m(5) * (pkin(3) * t15 + t9);
t34 = 2 * mrSges(5,1);
t33 = -2 * mrSges(6,1);
t32 = 2 * qJ(2);
t31 = m(5) + m(6);
t16 = -pkin(3) - pkin(4);
t30 = t25 * t17 ^ 2;
t28 = mrSges(5,1) + mrSges(6,1);
t27 = -mrSges(5,2) + mrSges(6,3);
t26 = (mrSges(5,3) + mrSges(6,2));
t24 = qJ(5) + t17;
t23 = m(5) / 0.2e1 + m(4) / 0.2e1;
t21 = t15 * qJ(4) - qJ(2);
t19 = qJ(2) ^ 2;
t18 = qJ(4) ^ 2;
t10 = t15 * mrSges(6,2);
t5 = t14 * pkin(3) - t21;
t4 = t24 * t15;
t3 = t24 * t14;
t1 = t16 * t14 + t21;
t2 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t32) + 0.2e1 * t1 * t10 + Ifges(3,1) + Ifges(2,3) + m(5) * (t5 ^ 2 + t30) + m(6) * (t1 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(4) * (t19 + t30) + (m(3) * (pkin(1) ^ 2 + t19)) + ((mrSges(4,1) * t32) + t5 * t34 + t1 * t33 + t3 * t37 + (Ifges(6,2) + Ifges(5,3) + Ifges(4,2)) * t14) * t14 + ((mrSges(4,2) * t32) - 0.2e1 * t5 * mrSges(5,3) + t4 * t37 + (Ifges(6,1) + Ifges(4,1) + Ifges(5,1)) * t15 + 0.2e1 * (-Ifges(4,4) + Ifges(6,4) + Ifges(5,5)) * t14) * t15 - 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * t36; -(m(3) * pkin(1)) + mrSges(3,2) + m(6) * (t14 * t3 + t15 * t4) + 0.2e1 * t23 * t36 + t25 * (-mrSges(4,3) + t27); m(3) + 0.2e1 * (m(6) / 0.2e1 + t23) * t25; m(6) * (qJ(4) * t3 - t16 * t4) + t4 * mrSges(6,1) + t3 * mrSges(6,2) + (-pkin(3) * mrSges(5,2) - t16 * mrSges(6,3) + Ifges(5,4) + Ifges(4,5) - Ifges(6,5)) * t15 + (t27 * qJ(4) - Ifges(4,6) + Ifges(5,6) - Ifges(6,6)) * t14 + (t35 + (mrSges(4,1) + mrSges(5,1)) * t15 + (-mrSges(4,2) + mrSges(5,3)) * t14) * t17; (mrSges(4,1) + t28) * t15 + t35 + m(6) * (-t16 * t15 + t9) + (-mrSges(4,2) + t26) * t14; pkin(3) * t34 + t16 * t33 + Ifges(5,2) + Ifges(4,3) + Ifges(6,3) + 0.2e1 * t26 * qJ(4) + m(5) * (pkin(3) ^ 2 + t18) + m(6) * (t16 ^ 2 + t18); -m(6) * t4 + (-m(5) * t17 - t27) * t15; -t31 * t15; -m(5) * pkin(3) + m(6) * t16 - t28; t31; m(6) * t1 - t14 * mrSges(6,1) + t10; 0; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
