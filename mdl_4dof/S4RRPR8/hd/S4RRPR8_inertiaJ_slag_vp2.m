% Calculate joint inertia matrix for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR8_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR8_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:48
% EndTime: 2019-12-31 17:07:49
% DurationCPUTime: 0.25s
% Computational Cost: add. (168->81), mult. (320->107), div. (0->0), fcn. (228->4), ass. (0->26)
t16 = sin(qJ(2));
t18 = cos(qJ(2));
t32 = t16 ^ 2 + t18 ^ 2;
t31 = -m(4) * pkin(2) - mrSges(4,1);
t24 = t16 * qJ(3) + pkin(1);
t8 = -t18 * pkin(2) - t24;
t30 = -0.2e1 * t8;
t29 = pkin(2) + pkin(3);
t28 = pkin(5) - pkin(6);
t15 = sin(qJ(4));
t17 = cos(qJ(4));
t6 = -t15 * qJ(3) - t17 * t29;
t27 = t6 * mrSges(5,1);
t7 = t17 * qJ(3) - t15 * t29;
t26 = t7 * mrSges(5,2);
t25 = t32 * pkin(5) ^ 2;
t22 = t17 * mrSges(5,1) - t15 * mrSges(5,2);
t10 = t28 * t18;
t9 = t28 * t16;
t1 = -t15 * t10 + t17 * t9;
t2 = t17 * t10 + t15 * t9;
t4 = -t16 * t15 - t18 * t17;
t5 = -t18 * t15 + t16 * t17;
t21 = t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t5 + Ifges(5,6) * t4;
t3 = t29 * t18 + t24;
t11 = [Ifges(2,3) + t4 * (Ifges(5,4) * t5 + Ifges(5,2) * t4) + t5 * (Ifges(5,1) * t5 + Ifges(5,4) * t4) + 0.2e1 * t3 * (-t4 * mrSges(5,1) + t5 * mrSges(5,2)) + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t30 + (Ifges(4,3) + Ifges(3,2)) * t18) * t18 + m(5) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(4) * (t8 ^ 2 + t25) + m(3) * (pkin(1) ^ 2 + t25) + 0.2e1 * (-t1 * t5 + t2 * t4) * mrSges(5,3) + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * pkin(5) * t32 + (-0.2e1 * pkin(1) * mrSges(3,2) + mrSges(4,3) * t30 + 0.2e1 * (Ifges(3,4) - Ifges(4,5)) * t18 + (Ifges(4,1) + Ifges(3,1)) * t16) * t16; m(5) * (t6 * t1 + t7 * t2) + (t7 * t4 - t6 * t5) * mrSges(5,3) + (qJ(3) * mrSges(4,2) + Ifges(3,6) - Ifges(4,6)) * t18 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t16 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t18 + (-mrSges(3,1) + t31) * t16) * pkin(5) - t21; 0.2e1 * pkin(2) * mrSges(4,1) - 0.2e1 * t27 + 0.2e1 * t26 + 0.2e1 * qJ(3) * mrSges(4,3) + Ifges(4,2) + Ifges(3,3) + Ifges(5,3) + m(4) * (pkin(2) ^ 2 + qJ(3) ^ 2) + m(5) * (t6 ^ 2 + t7 ^ 2); m(5) * (t17 * t1 + t15 * t2) + (m(4) * pkin(5) + mrSges(4,2)) * t16 + (t15 * t4 - t17 * t5) * mrSges(5,3); m(5) * (t15 * t7 + t17 * t6) - t22 + t31; m(4) + m(5) * (t15 ^ 2 + t17 ^ 2); t21; -Ifges(5,3) - t26 + t27; t22; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t11(1), t11(2), t11(4), t11(7); t11(2), t11(3), t11(5), t11(8); t11(4), t11(5), t11(6), t11(9); t11(7), t11(8), t11(9), t11(10);];
Mq = res;
