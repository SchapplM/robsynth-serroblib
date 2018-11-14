% Calculate joint inertia matrix for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:55
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4RRRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP1_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP1_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:54:20
% EndTime: 2018-11-14 13:54:20
% DurationCPUTime: 0.12s
% Computational Cost: add. (112->47), mult. (206->54), div. (0->0), fcn. (111->4), ass. (0->31)
t17 = cos(qJ(2));
t11 = t17 * pkin(1) + pkin(2);
t14 = sin(qJ(3));
t16 = cos(qJ(3));
t15 = sin(qJ(2));
t31 = pkin(1) * t15;
t6 = t16 * t11 - t14 * t31;
t4 = pkin(3) + t6;
t2 = t4 * mrSges(5,1);
t3 = t6 * mrSges(4,1);
t34 = t2 + t3;
t33 = mrSges(4,2) + mrSges(5,2);
t32 = m(5) * pkin(3);
t30 = t14 * pkin(2);
t29 = t16 * pkin(2);
t27 = Ifges(4,3) + Ifges(5,3);
t26 = Ifges(3,3) + t27;
t7 = t14 * t11 + t16 * t31;
t24 = t33 * t7;
t12 = mrSges(4,1) * t29;
t10 = pkin(3) + t29;
t9 = t10 * mrSges(5,1);
t23 = t12 + t9 + t27;
t22 = t33 * t30;
t21 = (t17 * mrSges(3,1) - t15 * mrSges(3,2)) * pkin(1);
t19 = pkin(2) ^ 2;
t18 = pkin(3) * mrSges(5,1);
t13 = t14 ^ 2 * t19;
t5 = t7 ^ 2;
t1 = t7 * t30;
t8 = [Ifges(2,3) + 0.2e1 * t2 + 0.2e1 * t3 + m(4) * (t6 ^ 2 + t5) + m(5) * (t4 ^ 2 + t5) + m(3) * (t15 ^ 2 + t17 ^ 2) * pkin(1) ^ 2 + t26 - 0.2e1 * t24 + 0.2e1 * t21; Ifges(3,3) + t21 + m(4) * (t6 * t29 + t1) + m(5) * (t10 * t4 + t1) + t23 + t33 * (-t7 - t30) + t34; 0.2e1 * t12 + 0.2e1 * t9 - 0.2e1 * t22 + m(4) * (t16 ^ 2 * t19 + t13) + m(5) * (t10 ^ 2 + t13) + t26; t4 * t32 + t18 - t24 + t27 + t34; t10 * t32 + t18 - t22 + t23; m(5) * pkin(3) ^ 2 + 0.2e1 * t18 + t27; 0; 0; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t8(1) t8(2) t8(4) t8(7); t8(2) t8(3) t8(5) t8(8); t8(4) t8(5) t8(6) t8(9); t8(7) t8(8) t8(9) t8(10);];
Mq  = res;
