% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:08
% EndTime: 2019-12-31 17:45:08
% DurationCPUTime: 0.25s
% Computational Cost: add. (1038->50), mult. (1786->70), div. (0->0), fcn. (1698->6), ass. (0->39)
t42 = sin(pkin(8));
t43 = cos(pkin(8));
t44 = cos(qJ(5));
t67 = sin(qJ(5));
t33 = -t67 * t42 + t44 * t43;
t25 = t33 ^ 2;
t32 = -t44 * t42 - t67 * t43;
t57 = t32 ^ 2 + t25;
t75 = -t57 / 0.2e1;
t77 = m(6) * (t57 / 0.2e1 + t75);
t74 = qJD(3) * t77;
t73 = qJD(2) * t77;
t36 = -cos(pkin(7)) * pkin(1) - pkin(2) - qJ(4);
t68 = -pkin(6) + t36;
t64 = t42 ^ 2 + t43 ^ 2;
t62 = sin(pkin(7)) * pkin(1) + qJ(3);
t35 = t42 * pkin(4) + t62;
t53 = t32 * mrSges(6,1) - t33 * mrSges(6,2);
t10 = t42 * mrSges(5,1) + t43 * mrSges(5,2) + mrSges(4,3) + m(6) * t35 + 0.4e1 * (m(4) / 0.4e1 + m(5) / 0.4e1) * t62 - t53;
t61 = t10 * qJD(1);
t54 = m(5) * t64;
t48 = m(6) * t75 - t54 / 0.2e1;
t56 = -m(5) / 0.2e1 - m(6) / 0.2e1;
t12 = t48 + t56;
t60 = t12 * qJD(1);
t20 = t33 * mrSges(6,1) + t32 * mrSges(6,2);
t59 = t20 * qJD(1);
t58 = t20 * qJD(5);
t1 = -t25 * Ifges(6,4) + t35 * t20 + (Ifges(6,4) * t32 + (Ifges(6,1) - Ifges(6,2)) * t33) * t32;
t50 = t1 * qJD(1);
t26 = t68 * t42;
t27 = t68 * t43;
t18 = -t67 * t26 + t44 * t27;
t19 = t44 * t26 + t67 * t27;
t7 = m(6) * (-t33 * t18 + t19 * t32) + t64 * mrSges(5,3) + t57 * mrSges(6,3) - t36 * t54;
t49 = t7 * qJD(1);
t11 = t48 - t56;
t8 = qJD(5) * t77;
t2 = [t10 * qJD(3) + t7 * qJD(4) + t1 * qJD(5), 0, t11 * qJD(4) + t61, t11 * qJD(3) + t49, (-t19 * mrSges(6,1) - t18 * mrSges(6,2) + Ifges(6,5) * t32 - Ifges(6,6) * t33) * qJD(5) + t50; 0, 0, t8, 0, -t58 + t74; t12 * qJD(4) - t61, t8, 0, t60, t53 * qJD(5) + t73; -t12 * qJD(3) - t49 + t58, 0, -t60, 0, t59; -t20 * qJD(4) - t50, -t74, -t73, -t59, 0;];
Cq = t2;
