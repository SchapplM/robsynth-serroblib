% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-28 15:34
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-28 15:34:03
% EndTime: 2019-05-28 15:34:03
% DurationCPUTime: 0.28s
% Computational Cost: add. (857->71), mult. (1345->91), div. (0->0), fcn. (863->4), ass. (0->50)
t42 = sin(qJ(4));
t44 = cos(qJ(4));
t56 = t42 * mrSges(5,1) + t44 * mrSges(5,2);
t72 = qJD(3) * t56;
t71 = t56 * qJD(4);
t69 = m(4) * pkin(1);
t43 = sin(qJ(2));
t68 = t43 * pkin(1);
t45 = cos(qJ(2));
t67 = t45 * pkin(1);
t60 = -pkin(2) - t67;
t38 = -pkin(3) + t60;
t40 = qJ(3) + t68;
t26 = t38 * t42 + t40 * t44;
t21 = t26 * t44;
t30 = (-t42 * t45 + t43 * t44) * pkin(1);
t66 = t30 * mrSges(5,1);
t46 = -pkin(2) - pkin(3);
t36 = qJ(3) * t44 + t42 * t46;
t31 = t36 * t44;
t25 = t38 * t44 - t40 * t42;
t51 = (t42 * t43 + t44 * t45) * pkin(1);
t29 = mrSges(5,2) * t51;
t49 = t29 - t66 + (-mrSges(3,1) - mrSges(4,1)) * t68 + (-mrSges(3,2) + mrSges(4,3)) * t67;
t3 = -m(5) * (t25 * t30 + t26 * t51) - (t40 * t45 + t43 * t60) * t69 - t49;
t65 = t3 * qJD(1);
t24 = t26 * mrSges(5,1);
t7 = t25 * mrSges(5,2) + t24;
t64 = t7 * qJD(1);
t63 = t7 * qJD(4);
t34 = t36 * mrSges(5,1);
t35 = -qJ(3) * t42 + t44 * t46;
t9 = t35 * mrSges(5,2) + t34;
t62 = t9 * qJD(4);
t52 = mrSges(4,3) + t56;
t11 = m(5) * (-t25 * t42 + t21) + m(4) * t40 + t52;
t61 = t11 * qJD(1);
t57 = -t29 / 0.2e1 + t66 / 0.2e1;
t50 = (t35 / 0.2e1 + t25 / 0.2e1) * mrSges(5,2) + t24 / 0.2e1 + t34 / 0.2e1;
t1 = t50 - t57;
t55 = t1 * qJD(1) + t9 * qJD(2);
t15 = m(4) * qJ(3) + m(5) * (-t35 * t42 + t31) + t52;
t47 = m(4) * (0.4e1 * qJ(3) + 0.2e1 * t68) / 0.4e1 + m(5) * (t21 + t31 + (-t25 - t35) * t42) / 0.2e1 + t52;
t48 = -m(5) * (t44 * t30 + t42 * t51) / 0.2e1 - m(4) * t68 / 0.2e1;
t5 = t47 + t48;
t54 = t5 * qJD(1) + t15 * qJD(2);
t53 = (-qJD(1) - qJD(2)) * t56;
t4 = t47 - t48;
t2 = t50 + t57;
t6 = [-qJD(2) * t3 + qJD(3) * t11 + t63, -t65 + (m(5) * (t35 * t30 + t36 * t51) + (-pkin(2) * t43 + qJ(3) * t45) * t69 + t49) * qJD(2) + t4 * qJD(3) + t2 * qJD(4), t4 * qJD(2) + t61, t2 * qJD(2) - t63 + t64; qJD(3) * t5 + qJD(4) * t1 + t65, qJD(3) * t15 + t62, t54, t55 - t62; -qJD(2) * t5 - t61 + t71, -t54 + t71, 0, -t53 - t71; -qJD(2) * t1 - t64 - t72, -t55 - t72, t53, 0;];
Cq  = t6;
