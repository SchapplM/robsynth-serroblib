% Calculate vector of centrifugal and coriolis load on the joints for
% S4PPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:57
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4PPRP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP2_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:57:22
% EndTime: 2018-11-14 13:57:22
% DurationCPUTime: 0.21s
% Computational Cost: add. (245->16), mult. (264->36), div. (0->0), fcn. (130->2), ass. (0->15)
t43 = (rSges(5,1) + pkin(3));
t44 = -qJD(3) * t43 + 2 * qJD(4);
t28 = pkin(5) + qJ(3);
t27 = cos(t28);
t42 = rSges(5,3) + qJ(4);
t39 = t42 * t27;
t26 = sin(t28);
t41 = t42 * t26;
t12 = t26 * rSges(4,1) + t27 * rSges(4,2);
t33 = qJD(3) * t12;
t15 = t27 * rSges(4,1) - t26 * rSges(4,2);
t32 = qJD(3) * t15;
t2 = (-qJD(3) * t41 + t44 * t27) * qJD(3);
t1 = (qJD(3) * t39 + t44 * t26) * qJD(3);
t3 = [-m(4) * qJD(3) * t33 + m(5) * t1; 0; ((0.2e1 * t12 * t32 - t15 * t33) * qJD(3) - t32 * t33) * m(4) + (t2 * (-t43 * t26 + t39) + t1 * (t43 * t27 + t41)) * m(5); 0.2e1 * (-t1 * t27 / 0.2e1 + t2 * t26 / 0.2e1) * m(5);];
tauc  = t3(:);
