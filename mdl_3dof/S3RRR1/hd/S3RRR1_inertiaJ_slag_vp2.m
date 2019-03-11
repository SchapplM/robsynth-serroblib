% Calculate joint inertia matrix for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [3x3]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S3RRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_inertiaJ_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_inertiaJ_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRR1_inertiaJ_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRR1_inertiaJ_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:07:57
% EndTime: 2019-03-08 18:07:57
% DurationCPUTime: 0.08s
% Computational Cost: add. (48->26), mult. (106->36), div. (0->0), fcn. (56->4), ass. (0->17)
t7 = sin(qJ(2));
t18 = pkin(1) * t7;
t9 = cos(qJ(2));
t4 = pkin(1) * t9 + pkin(2);
t6 = sin(qJ(3));
t8 = cos(qJ(3));
t3 = t8 * t18 + t4 * t6;
t17 = t3 * mrSges(4,2);
t16 = t6 * mrSges(4,2);
t15 = Ifges(3,3) + Ifges(4,3);
t14 = pkin(2) * t16;
t2 = -t6 * t18 + t4 * t8;
t1 = t2 * mrSges(4,1);
t13 = Ifges(4,3) + t1 - t17;
t12 = (t9 * mrSges(3,1) - t7 * mrSges(3,2)) * pkin(1);
t5 = t8 * pkin(2) * mrSges(4,1);
t10 = [-0.2e1 * t17 + Ifges(2,3) + 0.2e1 * t1 + 0.2e1 * t12 + m(4) * (t2 ^ 2 + t3 ^ 2) + m(3) * (t7 ^ 2 + t9 ^ 2) * pkin(1) ^ 2 + t15; Ifges(3,3) + t5 + t12 + (m(4) * (t2 * t8 + t3 * t6) - t16) * pkin(2) + t13; -0.2e1 * t14 + 0.2e1 * t5 + m(4) * (t6 ^ 2 + t8 ^ 2) * pkin(2) ^ 2 + t15; t13; Ifges(4,3) + t5 - t14; Ifges(4,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t10(1) t10(2) t10(4); t10(2) t10(3) t10(5); t10(4) t10(5) t10(6);];
Mq  = res;
