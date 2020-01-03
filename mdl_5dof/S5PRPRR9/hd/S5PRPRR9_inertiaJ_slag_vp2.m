% Calculate joint inertia matrix for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR9_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR9_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:37
% EndTime: 2019-12-31 17:39:38
% DurationCPUTime: 0.18s
% Computational Cost: add. (158->63), mult. (273->81), div. (0->0), fcn. (149->4), ass. (0->24)
t18 = cos(qJ(5));
t32 = t18 ^ 2;
t16 = sin(qJ(5));
t6 = -t18 * mrSges(6,1) + t16 * mrSges(6,2);
t31 = -mrSges(5,1) + t6;
t27 = t16 ^ 2 + t32;
t25 = t27 * mrSges(6,3);
t30 = -0.2e1 * t6;
t17 = sin(qJ(4));
t19 = cos(qJ(4));
t20 = -pkin(2) - pkin(3);
t4 = -t17 * qJ(3) + t19 * t20;
t29 = t4 * mrSges(5,1);
t5 = t19 * qJ(3) + t17 * t20;
t28 = t5 * mrSges(5,2);
t3 = -pkin(7) + t5;
t26 = t27 * t3;
t24 = t27 * t17;
t23 = -mrSges(6,1) * t16 - mrSges(6,2) * t18;
t22 = Ifges(6,2) * t32 + Ifges(5,3) + (Ifges(6,1) * t16 + 0.2e1 * Ifges(6,4) * t18) * t16;
t15 = t19 ^ 2;
t13 = t17 ^ 2;
t2 = pkin(4) - t4;
t1 = [m(6) * t27 + m(2) + m(3) + m(4) + m(5); 0; (2 * pkin(2) * mrSges(4,1)) - 0.2e1 * t29 + 0.2e1 * t28 + 0.2e1 * qJ(3) * mrSges(4,3) + t2 * t30 + Ifges(4,2) + Ifges(3,3) - 0.2e1 * t3 * t25 + m(6) * (t27 * t3 ^ 2 + t2 ^ 2) + m(5) * (t4 ^ 2 + t5 ^ 2) + m(4) * ((pkin(2) ^ 2) + qJ(3) ^ 2) + t22; 0; -m(4) * pkin(2) - mrSges(4,1) + t31 * t19 + (mrSges(5,2) - t25) * t17 + m(6) * (-t19 * t2 + t3 * t24) + m(5) * (t17 * t5 + t19 * t4); m(4) + m(5) * (t13 + t15) + m(6) * (t27 * t13 + t15); 0; m(6) * (-pkin(4) * t2 + pkin(7) * t26) - t28 + t29 + (pkin(4) + t2) * t6 + (-t27 * pkin(7) + t26) * mrSges(6,3) - t22; -t17 * mrSges(5,2) + (m(6) * pkin(7) + mrSges(6,3)) * t24 + (m(6) * pkin(4) - t31) * t19; pkin(4) * t30 + m(6) * (t27 * pkin(7) ^ 2 + pkin(4) ^ 2) + 0.2e1 * pkin(7) * t25 + t22; t6; (-mrSges(6,2) * t3 - Ifges(6,6)) * t18 + (-mrSges(6,1) * t3 - Ifges(6,5)) * t16; t23 * t17; Ifges(6,5) * t16 + Ifges(6,6) * t18 + pkin(7) * t23; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
