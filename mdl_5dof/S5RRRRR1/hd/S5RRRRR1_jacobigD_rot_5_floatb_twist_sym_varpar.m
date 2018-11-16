% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
%
% Output:
% JgD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JgD_rot = S5RRRRR1_jacobigD_rot_5_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_jacobigD_rot_5_floatb_twist_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_jacobigD_rot_5_floatb_twist_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_jacobigD_rot_5_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 14:52:30
% EndTime: 2018-11-16 14:52:30
% DurationCPUTime: 0.02s
% Computational Cost: add. (25->11), mult. (15->8), div. (0->0), fcn. (15->4), ass. (0->9)
t118 = qJD(2) + qJD(3) + qJD(4);
t119 = qJ(2) + qJ(3) + qJ(4);
t124 = cos(t119) * t118;
t120 = sin(qJ(1));
t123 = qJD(1) * t120;
t121 = cos(qJ(1));
t122 = qJD(1) * t121;
t116 = sin(t119);
t1 = [0, -t122, -t122, -t122, -t116 * t123 + t121 * t124; 0, -t123, -t123, -t123, t116 * t122 + t120 * t124; 0, 0, 0, 0, -t118 * t116;];
JgD_rot  = t1;
