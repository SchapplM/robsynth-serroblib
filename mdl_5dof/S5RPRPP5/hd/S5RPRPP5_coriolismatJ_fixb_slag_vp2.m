% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:09
% EndTime: 2019-12-31 18:16:10
% DurationCPUTime: 0.32s
% Computational Cost: add. (468->94), mult. (821->102), div. (0->0), fcn. (474->2), ass. (0->53)
t34 = sin(qJ(3));
t27 = qJ(4) * t34;
t35 = cos(qJ(3));
t62 = pkin(3) + pkin(4);
t55 = t62 * t35;
t16 = -t27 - t55;
t26 = t35 * qJ(4);
t68 = -t62 * t34 + t26;
t40 = -t34 * pkin(3) + t26;
t67 = m(5) * t40;
t28 = t35 * mrSges(5,3);
t66 = t35 * mrSges(4,2) - t28;
t65 = Ifges(6,4) + Ifges(5,5) - Ifges(4,4);
t64 = m(6) / 0.2e1;
t63 = m(5) + m(6);
t21 = qJ(2) - t40;
t61 = m(5) * t21;
t13 = qJ(2) - t68;
t60 = m(6) * t13;
t59 = m(6) * t35;
t57 = t34 * mrSges(5,1);
t54 = -mrSges(5,2) + mrSges(6,3);
t32 = t34 ^ 2;
t33 = t35 ^ 2;
t53 = t33 + t32;
t29 = t35 * mrSges(6,2);
t41 = -t28 + t57 + t61;
t1 = t41 * t27 + (t29 - t60) * t16 + (-t16 * mrSges(6,1) - qJ(2) * mrSges(4,2) + t13 * mrSges(6,2) + t21 * mrSges(5,3) - t65 * t34) * t34 + (t41 * pkin(3) + qJ(2) * mrSges(4,1) + t21 * mrSges(5,1) + t13 * mrSges(6,1) + (-Ifges(4,1) - Ifges(5,1) - Ifges(6,1) + Ifges(4,2) + Ifges(6,2) + Ifges(5,3)) * t34 + t65 * t35) * t35;
t52 = t1 * qJD(1);
t39 = -t29 + t66;
t43 = mrSges(5,1) + mrSges(6,1) + mrSges(4,1);
t2 = mrSges(3,3) + (m(4) + m(3)) * qJ(2) + t61 + t60 + t43 * t34 + t39;
t51 = t2 * qJD(1);
t3 = t13 * t59 + (t34 * mrSges(6,1) - t29 + t41) * t35;
t50 = t3 * qJD(1);
t37 = -pkin(1) - pkin(6);
t47 = qJ(5) + t37;
t19 = t47 * t34;
t20 = t47 * t35;
t4 = m(6) * (t34 * t19 + t35 * t20) + t53 * mrSges(6,3);
t49 = t4 * qJD(1);
t5 = t35 * mrSges(6,1) + t34 * mrSges(6,2) + 0.2e1 * (t55 / 0.4e1 + t27 / 0.4e1 - t16 / 0.4e1) * m(6);
t48 = t5 * qJD(1);
t46 = qJD(3) * t34;
t14 = (t33 / 0.2e1 + t32 / 0.2e1 + 0.1e1 / 0.2e1) * m(6);
t45 = t14 * qJD(1);
t25 = t63 * qJ(4) + mrSges(6,2) + mrSges(5,3);
t44 = t25 * qJD(3);
t42 = qJD(1) * t59;
t24 = t63 * t34;
t15 = t53 * t64 - m(6) / 0.2e1;
t8 = m(6) * t19 + (m(5) * t37 + t54) * t34;
t6 = [t2 * qJD(2) + t1 * qJD(3) - t3 * qJD(4) + t4 * qJD(5), t15 * qJD(5) + t51, t52 + t8 * qJD(4) + (pkin(3) * mrSges(5,2) - mrSges(6,3) * t62 - Ifges(5,4) - Ifges(4,5) + Ifges(6,5)) * t46 + (-t19 * mrSges(6,1) + t20 * mrSges(6,2) + m(6) * (qJ(4) * t20 - t19 * t62) + (t54 * qJ(4) - Ifges(4,6) + Ifges(5,6) - Ifges(6,6)) * t35 + (-t34 * mrSges(4,1) - t57 - t66 + t67) * t37) * qJD(3), t8 * qJD(3) - t50, t15 * qJD(2) + t49; t14 * qJD(5) - t51, 0, t24 * qJD(4) - t43 * t46 + (0.2e1 * t68 * t64 - t39 + t67) * qJD(3), t24 * qJD(3), t45; t5 * qJD(5) - t52, 0, t25 * qJD(4), t44, t48; -qJD(5) * t59 + t50, 0, -t44, 0, -t42; -t14 * qJD(2) - t5 * qJD(3) + qJD(4) * t59 - t49, -t45, -t48, t42, 0;];
Cq = t6;
