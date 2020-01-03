% Calculate joint inertia matrix for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:25
% EndTime: 2019-12-31 17:33:25
% DurationCPUTime: 0.12s
% Computational Cost: add. (65->40), mult. (143->49), div. (0->0), fcn. (64->4), ass. (0->19)
t10 = sin(qJ(5));
t12 = cos(qJ(5));
t7 = t12 ^ 2;
t19 = t10 ^ 2 + t7;
t1 = m(6) * t19;
t20 = m(5) + t1;
t21 = t19 * mrSges(6,3) - mrSges(5,2);
t2 = t10 * mrSges(6,1) + t12 * mrSges(6,2);
t18 = mrSges(5,3) + t2;
t14 = -pkin(3) - pkin(6);
t17 = t19 * t14;
t16 = t12 * mrSges(6,1) - t10 * mrSges(6,2);
t15 = qJ(4) ^ 2;
t13 = cos(qJ(3));
t11 = sin(qJ(3));
t8 = t13 ^ 2;
t6 = t11 ^ 2;
t4 = qJ(4) * t11;
t3 = [m(2) + m(3) + m(4) + t20; 0; m(3) + m(6) * (t19 * t8 + t6) + 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (t6 + t8); 0; (mrSges(4,1) + t21) * t13 + m(5) * (pkin(3) * t13 + t4) + m(6) * (-t13 * t17 + t4) + (-mrSges(4,2) + t18) * t11; Ifges(6,1) * t7 - 0.2e1 * pkin(3) * mrSges(5,2) + Ifges(5,1) + Ifges(4,3) + (-0.2e1 * Ifges(6,4) * t12 + Ifges(6,2) * t10) * t10 + m(6) * (t19 * t14 ^ 2 + t15) + m(5) * (pkin(3) ^ 2 + t15) - 0.2e1 * mrSges(6,3) * t17 + 0.2e1 * t18 * qJ(4); 0; -t20 * t13; -m(5) * pkin(3) + t14 * t1 - t21; t20; t2; -t16 * t13; (mrSges(6,1) * t14 + Ifges(6,5)) * t12 + (-mrSges(6,2) * t14 - Ifges(6,6)) * t10; t16; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
