% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPPRR7_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:57
% EndTime: 2019-02-26 21:31:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->16), mult. (78->39), div. (0->0), fcn. (80->8), ass. (0->22)
t150 = sin(pkin(6));
t154 = sin(qJ(1));
t170 = t150 * t154;
t156 = cos(qJ(2));
t169 = t150 * t156;
t157 = cos(qJ(1));
t168 = t150 * t157;
t153 = sin(qJ(2));
t167 = t154 * t153;
t166 = t154 * t156;
t165 = t156 * t157;
t164 = t157 * t153;
t163 = qJD(1) * t150;
t152 = sin(qJ(5));
t162 = qJD(2) * t152;
t151 = cos(pkin(6));
t161 = t151 * t165 - t167;
t160 = t151 * t166 + t164;
t159 = t151 * t164 + t166;
t158 = -t151 * t167 + t165;
t155 = cos(qJ(5));
t1 = [0, t157 * t163, 0, 0, -t159 * qJD(1) - t160 * qJD(2) (-t152 * t170 + t160 * t155) * qJD(5) + t158 * t162 + (t161 * t152 + t155 * t168) * qJD(1); 0, t154 * t163, 0, 0, t158 * qJD(1) + t161 * qJD(2) (t152 * t168 - t161 * t155) * qJD(5) + t159 * t162 + (t160 * t152 + t155 * t170) * qJD(1); 0, 0, 0, 0, qJD(2) * t169, t150 * t153 * t162 + (-t151 * t152 - t155 * t169) * qJD(5);];
JgD_rot  = t1;
