% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRP7_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_jacobigD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:49:25
% EndTime: 2019-02-26 21:49:25
% DurationCPUTime: 0.03s
% Computational Cost: add. (15->7), mult. (48->12), div. (0->0), fcn. (48->6), ass. (0->12)
t169 = qJD(2) - qJD(4);
t160 = sin(qJ(1));
t168 = qJD(1) * t160;
t163 = cos(qJ(1));
t167 = qJD(1) * t163;
t158 = sin(qJ(4));
t159 = sin(qJ(2));
t161 = cos(qJ(4));
t162 = cos(qJ(2));
t166 = t158 * t162 - t159 * t161;
t164 = t169 * (t158 * t159 + t161 * t162);
t1 = [0, t167, 0, -t167, -t164 * t163 - t166 * t168, 0; 0, t168, 0, -t168, -t164 * t160 + t166 * t167, 0; 0, 0, 0, 0, t169 * t166, 0;];
JgD_rot  = t1;
