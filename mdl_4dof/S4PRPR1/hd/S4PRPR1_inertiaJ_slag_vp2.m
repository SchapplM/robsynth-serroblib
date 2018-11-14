% Calculate joint inertia matrix for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PRPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR1_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR1_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:42:12
% EndTime: 2018-11-14 13:42:12
% DurationCPUTime: 0.07s
% Computational Cost: add. (39->25), mult. (55->27), div. (0->0), fcn. (24->2), ass. (0->9)
t3 = sin(qJ(4));
t4 = cos(qJ(4));
t5 = -pkin(2) - pkin(3);
t1 = -t3 * qJ(3) + t4 * t5;
t8 = t1 * mrSges(5,1);
t2 = t4 * qJ(3) + t3 * t5;
t7 = t2 * mrSges(5,2);
t6 = t4 * mrSges(5,1) - t3 * mrSges(5,2);
t9 = [m(2) + m(3) + m(4) + m(5); 0; (2 * pkin(2) * mrSges(4,1)) - 0.2e1 * t8 + 0.2e1 * t7 + 0.2e1 * qJ(3) * mrSges(4,3) + Ifges(4,2) + Ifges(3,3) + Ifges(5,3) + m(4) * ((pkin(2) ^ 2) + qJ(3) ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2); 0; -mrSges(4,1) - m(4) * pkin(2) + m(5) * (t1 * t4 + t2 * t3) - t6; m(4) + m(5) * (t3 ^ 2 + t4 ^ 2); 0; -Ifges(5,3) - t7 + t8; t6; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t9(1) t9(2) t9(4) t9(7); t9(2) t9(3) t9(5) t9(8); t9(4) t9(5) t9(6) t9(9); t9(7) t9(8) t9(9) t9(10);];
Mq  = res;
