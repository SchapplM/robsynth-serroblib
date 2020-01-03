% Calculate joint inertia matrix for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:03
% EndTime: 2019-12-31 17:38:03
% DurationCPUTime: 0.14s
% Computational Cost: add. (156->66), mult. (296->95), div. (0->0), fcn. (207->6), ass. (0->27)
t18 = sin(pkin(8));
t19 = cos(pkin(8));
t21 = sin(qJ(2));
t23 = cos(qJ(2));
t2 = t21 * t18 + t23 * t19;
t32 = t2 ^ 2;
t22 = cos(qJ(5));
t17 = t22 ^ 2;
t31 = t19 * t2;
t20 = sin(qJ(5));
t10 = t22 * mrSges(6,1) - t20 * mrSges(6,2);
t30 = mrSges(5,1) + t10;
t24 = -pkin(2) - pkin(3);
t9 = t19 * qJ(3) + t18 * t24;
t29 = t20 ^ 2 + t17;
t28 = t18 * t29;
t27 = t29 * mrSges(6,3);
t26 = mrSges(5,2) - t27;
t25 = -mrSges(6,1) * t20 - mrSges(6,2) * t22;
t8 = -t18 * qJ(3) + t19 * t24;
t15 = t19 ^ 2;
t14 = t18 ^ 2;
t7 = -pkin(6) + t9;
t6 = pkin(4) - t8;
t4 = -t23 * t18 + t21 * t19;
t1 = t4 ^ 2;
t3 = [m(2) + m(5) * (t1 + t32) + m(6) * (t29 * t1 + t32) + 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * (t21 ^ 2 + t23 ^ 2); (mrSges(3,1) + mrSges(4,1)) * t23 + (-mrSges(3,2) + mrSges(4,3)) * t21 + t30 * t2 + t26 * t4 + m(4) * (t23 * pkin(2) + t21 * qJ(3)) + m(5) * (-t8 * t2 + t9 * t4) + m(6) * (t29 * t7 * t4 + t6 * t2); t17 * Ifges(6,2) + 0.2e1 * pkin(2) * mrSges(4,1) - 0.2e1 * t8 * mrSges(5,1) + 0.2e1 * t9 * mrSges(5,2) + 0.2e1 * qJ(3) * mrSges(4,3) + 0.2e1 * t6 * t10 + Ifges(4,2) + Ifges(3,3) + Ifges(5,3) + (Ifges(6,1) * t20 + 0.2e1 * Ifges(6,4) * t22) * t20 - 0.2e1 * t7 * t27 + m(6) * (t29 * t7 ^ 2 + t6 ^ 2) + m(4) * (pkin(2) ^ 2 + qJ(3) ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2); -m(4) * t23 + m(5) * (t18 * t4 - t31) + m(6) * (t4 * t28 - t31); -m(4) * pkin(2) - mrSges(4,1) - t30 * t19 + t26 * t18 + m(6) * (-t19 * t6 + t7 * t28) + m(5) * (t18 * t9 + t19 * t8); m(4) + m(5) * (t14 + t15) + m(6) * (t29 * t14 + t15); 0; 0; 0; m(6) * t29 + m(5); t25 * t4; (-mrSges(6,2) * t7 - Ifges(6,6)) * t22 + (-mrSges(6,1) * t7 - Ifges(6,5)) * t20; t25 * t18; t10; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
