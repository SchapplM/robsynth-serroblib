% Calculate joint inertia matrix for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:06
% EndTime: 2019-12-31 17:45:07
% DurationCPUTime: 0.16s
% Computational Cost: add. (192->58), mult. (329->83), div. (0->0), fcn. (258->6), ass. (0->27)
t24 = sin(pkin(7));
t14 = t24 * pkin(1) + qJ(3);
t37 = t14 ^ 2;
t25 = cos(pkin(8));
t22 = t25 ^ 2;
t26 = cos(pkin(7));
t18 = -t26 * pkin(1) - pkin(2);
t13 = -qJ(4) + t18;
t35 = -pkin(6) + t13;
t34 = cos(qJ(5));
t33 = sin(qJ(5));
t23 = sin(pkin(8));
t32 = t23 * mrSges(5,1) + t25 * mrSges(5,2);
t31 = t23 ^ 2 + t22;
t7 = -t34 * t23 - t33 * t25;
t9 = -t33 * t23 + t34 * t25;
t30 = t7 ^ 2 + t9 ^ 2;
t10 = m(5) * t31;
t29 = t31 * mrSges(5,3);
t3 = -t7 * mrSges(6,1) + t9 * mrSges(6,2);
t28 = m(6) * t30 + m(4) + t10;
t11 = t23 * pkin(4) + t14;
t5 = t35 * t25;
t4 = t35 * t23;
t2 = t33 * t5 + t34 * t4;
t1 = -t33 * t4 + t34 * t5;
t6 = [Ifges(4,1) + Ifges(2,3) + Ifges(3,3) + 0.2e1 * t11 * t3 + 0.2e1 * t18 * mrSges(4,2) + Ifges(5,1) * t22 + (-0.2e1 * Ifges(5,4) * t25 + Ifges(5,2) * t23) * t23 + (-0.2e1 * t1 * mrSges(6,3) + Ifges(6,1) * t9) * t9 + (0.2e1 * t2 * mrSges(6,3) + 0.2e1 * Ifges(6,4) * t9 + Ifges(6,2) * t7) * t7 + m(6) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2) + m(4) * (t18 ^ 2 + t37) + m(5) * (t31 * t13 ^ 2 + t37) + m(3) * (t24 ^ 2 + t26 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (mrSges(4,3) + t32) * t14 + 0.2e1 * (t26 * mrSges(3,1) - t24 * mrSges(3,2)) * pkin(1) - 0.2e1 * t13 * t29; m(6) * (t1 * t7 + t2 * t9); m(3) + t28; mrSges(4,2) - t30 * mrSges(6,3) - t29 + m(6) * (t9 * t1 - t7 * t2) + m(4) * t18 + t13 * t10; 0; t28; m(5) * t14 + m(6) * t11 + t3 + t32; 0; 0; m(5) + m(6); t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t9 + Ifges(6,6) * t7; -t3; t9 * mrSges(6,1) + t7 * mrSges(6,2); 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t6(1), t6(2), t6(4), t6(7), t6(11); t6(2), t6(3), t6(5), t6(8), t6(12); t6(4), t6(5), t6(6), t6(9), t6(13); t6(7), t6(8), t6(9), t6(10), t6(14); t6(11), t6(12), t6(13), t6(14), t6(15);];
Mq = res;
