% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:41
% EndTime: 2019-12-05 15:06:42
% DurationCPUTime: 0.33s
% Computational Cost: add. (454->68), mult. (1130->92), div. (0->0), fcn. (956->6), ass. (0->36)
t45 = sin(qJ(4));
t42 = t45 ^ 2;
t47 = cos(qJ(4));
t43 = t47 ^ 2;
t60 = t42 + t43;
t62 = -qJ(5) - pkin(6);
t36 = t62 * t45;
t37 = t62 * t47;
t76 = m(6) * (-t36 * t45 - t37 * t47);
t72 = m(6) * pkin(4);
t74 = -mrSges(6,1) - t72;
t44 = sin(pkin(8));
t46 = sin(qJ(3));
t58 = cos(pkin(8));
t67 = cos(qJ(3));
t34 = t67 * t44 + t46 * t58;
t69 = m(6) * t34;
t68 = t45 * pkin(4);
t63 = Ifges(6,4) + Ifges(5,4);
t61 = t45 * mrSges(6,1) + t47 * mrSges(6,2);
t33 = t46 * t44 - t67 * t58;
t1 = 0.4e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * (0.1e1 - t60) * t34 * t33;
t59 = t1 * qJD(1);
t23 = -m(6) * t68 - t61;
t57 = t23 * qJD(3);
t52 = t47 * pkin(4) + pkin(3);
t51 = t45 * mrSges(5,1) + t47 * mrSges(5,2);
t50 = -t47 * mrSges(6,1) + t45 * mrSges(6,2);
t10 = t60 * mrSges(6,3) + t76;
t8 = 0.2e1 * (t42 / 0.4e1 + t43 / 0.4e1 - 0.1e1 / 0.4e1) * t69;
t49 = t8 * qJD(1) + t10 * qJD(3);
t2 = t52 * t61 - t50 * t68 + (pkin(3) * mrSges(5,1) + t63 * t45 + t52 * t72) * t45 + (pkin(3) * mrSges(5,2) + (-Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,2)) * t45 - t63 * t47) * t47;
t48 = t2 * qJD(3);
t7 = (t60 + 0.1e1) * t69 / 0.2e1;
t4 = (t51 / 0.2e1 + t61 / 0.2e1 + (mrSges(5,2) / 0.2e1 + mrSges(6,2) / 0.2e1) * t47 + (t72 + mrSges(5,1) / 0.2e1 + mrSges(6,1) / 0.2e1) * t45) * t33;
t3 = [t1 * qJD(3), 0, t4 * qJD(4) + t7 * qJD(5) + t59 + ((-m(5) * pkin(3) - m(6) * t52 - t47 * mrSges(5,1) + t45 * mrSges(5,2) - mrSges(4,1) + t50) * t34 + (mrSges(4,2) - t76 + (-m(5) * pkin(6) - mrSges(5,3) - mrSges(6,3)) * t60) * t33) * qJD(3), t4 * qJD(3) + ((mrSges(5,2) + mrSges(6,2)) * t45 + (-mrSges(5,1) + t74) * t47) * qJD(4) * t34, t7 * qJD(3); 0, 0, 0, (-t51 + t23) * qJD(4), 0; t8 * qJD(5) - t59, 0, -t2 * qJD(4) + t10 * qJD(5), -t48 + (-t36 * mrSges(6,2) + (mrSges(5,2) * pkin(6) - Ifges(5,6) - Ifges(6,6)) * t45 + (-mrSges(5,1) * pkin(6) - mrSges(6,3) * pkin(4) + Ifges(5,5) + Ifges(6,5)) * t47 - t74 * t37) * qJD(4), t49; 0, 0, t23 * qJD(5) + t48, 0, t57; -t8 * qJD(3), 0, -t23 * qJD(4) - t49, -t57, 0;];
Cq = t3;
