% Calculate joint inertia matrix for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:17
% EndTime: 2019-12-31 17:32:17
% DurationCPUTime: 0.17s
% Computational Cost: add. (141->58), mult. (329->88), div. (0->0), fcn. (286->6), ass. (0->27)
t20 = sin(pkin(8));
t21 = cos(pkin(8));
t22 = sin(qJ(5));
t30 = cos(qJ(5));
t26 = -t22 * t20 + t30 * t21;
t8 = t30 * t20 + t22 * t21;
t1 = -mrSges(6,1) * t26 + t8 * mrSges(6,2);
t11 = -t21 * mrSges(5,1) + t20 * mrSges(5,2);
t34 = -t1 - t11;
t17 = t21 ^ 2;
t33 = -2 * mrSges(6,3);
t32 = m(5) + m(6);
t29 = pkin(6) + qJ(4);
t28 = t20 ^ 2 + t17;
t27 = qJ(4) * t28;
t24 = cos(qJ(3));
t23 = sin(qJ(3));
t19 = t24 ^ 2;
t18 = t23 ^ 2;
t14 = -t21 * pkin(4) - pkin(3);
t12 = t29 * t21;
t10 = t29 * t20;
t5 = t26 * t23;
t4 = t8 * t23;
t3 = -t22 * t10 + t30 * t12;
t2 = -t30 * t10 - t22 * t12;
t6 = [m(2) + m(3) + m(4) + m(5) * t28 + m(6) * (t26 ^ 2 + t8 ^ 2); m(6) * (t26 * t4 - t5 * t8); m(3) + m(4) * (t18 + t19) + m(5) * (t28 * t18 + t19) + m(6) * (t4 ^ 2 + t5 ^ 2 + t19); m(6) * (-t2 * t26 - t3 * t8); (t26 * t5 + t4 * t8) * mrSges(6,3) + (mrSges(4,1) + t34) * t24 + (t28 * mrSges(5,3) - mrSges(4,2)) * t23 + m(5) * (pkin(3) * t24 + t23 * t27) + m(6) * (-t14 * t24 - t2 * t4 + t3 * t5); Ifges(5,2) * t17 - 0.2e1 * pkin(3) * t11 + 0.2e1 * t14 * t1 + Ifges(4,3) + (Ifges(5,1) * t20 + 0.2e1 * Ifges(5,4) * t21) * t20 + (Ifges(6,1) * t8 + t2 * t33) * t8 + 0.2e1 * mrSges(5,3) * t27 + m(6) * (t14 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t28 * qJ(4) ^ 2 + pkin(3) ^ 2) - (-0.2e1 * Ifges(6,4) * t8 - Ifges(6,2) * t26 + t3 * t33) * t26; 0; -t32 * t24; -m(5) * pkin(3) + m(6) * t14 - t34; t32; t1; -t4 * mrSges(6,1) - t5 * mrSges(6,2); t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t8 + Ifges(6,6) * t26; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t6(1), t6(2), t6(4), t6(7), t6(11); t6(2), t6(3), t6(5), t6(8), t6(12); t6(4), t6(5), t6(6), t6(9), t6(13); t6(7), t6(8), t6(9), t6(10), t6(14); t6(11), t6(12), t6(13), t6(14), t6(15);];
Mq = res;
