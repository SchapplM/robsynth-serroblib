% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRP11_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_jacobigD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:14:59
% EndTime: 2019-02-26 22:14:59
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->16), mult. (78->39), div. (0->0), fcn. (80->8), ass. (0->22)
t151 = sin(pkin(6));
t154 = sin(qJ(2));
t171 = t151 * t154;
t155 = sin(qJ(1));
t170 = t151 * t155;
t158 = cos(qJ(1));
t169 = t151 * t158;
t168 = t154 * t155;
t167 = t154 * t158;
t157 = cos(qJ(2));
t166 = t155 * t157;
t165 = t158 * t157;
t164 = qJD(1) * t151;
t156 = cos(qJ(3));
t163 = qJD(2) * t156;
t152 = cos(pkin(6));
t162 = t152 * t165 - t168;
t161 = t152 * t166 + t167;
t160 = t152 * t167 + t166;
t159 = -t152 * t168 + t165;
t153 = sin(qJ(3));
t1 = [0, t158 * t164, t162 * qJD(1) + t159 * qJD(2), 0 (-t159 * t153 + t156 * t170) * qJD(3) - t161 * t163 + (t153 * t169 - t160 * t156) * qJD(1), 0; 0, t155 * t164, t161 * qJD(1) + t160 * qJD(2), 0 (-t160 * t153 - t156 * t169) * qJD(3) + t162 * t163 + (t153 * t170 + t159 * t156) * qJD(1), 0; 0, 0, qJD(2) * t171, 0, t151 * t157 * t163 + (t152 * t156 - t153 * t171) * qJD(3), 0;];
JgD_rot  = t1;
