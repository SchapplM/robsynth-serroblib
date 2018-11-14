% Calculate matrix of centrifugal and coriolis load on the joints for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
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
% Cq [3x3]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S3PRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR1_coriolismatJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:11:57
% EndTime: 2018-11-14 10:11:57
% DurationCPUTime: 0.07s
% Computational Cost: add. (82->12), mult. (228->20), div. (0->0), fcn. (208->4), ass. (0->12)
t14 = sin(qJ(3));
t15 = sin(qJ(2));
t16 = cos(qJ(3));
t17 = cos(qJ(2));
t10 = -t14 * t15 + t16 * t17;
t9 = -t14 * t17 - t16 * t15;
t19 = t9 * mrSges(4,1) - t10 * mrSges(4,2);
t23 = t19 * qJD(3);
t11 = (mrSges(4,1) * t14 + mrSges(4,2) * t16) * pkin(2);
t18 = t11 * qJD(2);
t8 = t11 * qJD(3);
t1 = [0 (-t15 * mrSges(3,1) - t17 * mrSges(3,2) + t19 + m(4) * (t10 * t14 + t16 * t9) * pkin(2)) * qJD(2) + t23, qJD(2) * t19 + t23; 0, -t8, -t8 - t18; 0, t18, 0;];
Cq  = t1;
