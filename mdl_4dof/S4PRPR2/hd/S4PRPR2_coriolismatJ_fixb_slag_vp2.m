% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:02
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4PRPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR2_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:02:16
% EndTime: 2018-11-14 14:02:16
% DurationCPUTime: 0.12s
% Computational Cost: add. (426->21), mult. (912->34), div. (0->0), fcn. (1032->6), ass. (0->19)
t25 = sin(pkin(6));
t26 = cos(pkin(6));
t37 = sin(qJ(2));
t38 = cos(qJ(2));
t20 = -t25 * t38 - t26 * t37;
t21 = -t25 * t37 + t26 * t38;
t27 = sin(qJ(4));
t28 = cos(qJ(4));
t13 = t28 * t20 - t27 * t21;
t30 = t27 * t20 + t28 * t21;
t2 = t13 * mrSges(5,1) - t30 * mrSges(5,2);
t39 = pkin(2) * t25;
t24 = t26 * pkin(2) + pkin(3);
t17 = t28 * t24 - t27 * t39;
t18 = t27 * t24 + t28 * t39;
t6 = t18 * mrSges(5,1) + t17 * mrSges(5,2);
t31 = t6 * qJD(4);
t29 = t6 * qJD(2);
t1 = [0, t2 * qJD(4) + (-t37 * mrSges(3,1) + t20 * mrSges(4,1) - t38 * mrSges(3,2) - t21 * mrSges(4,2) + m(5) * (t17 * t13 + t18 * t30) + m(4) * (t20 * t26 + t21 * t25) * pkin(2) + t2) * qJD(2), 0 (qJD(4) + qJD(2)) * t2; 0, -t31, 0, -t29 - t31; 0, 0, 0, 0; 0, t29, 0, 0;];
Cq  = t1;
