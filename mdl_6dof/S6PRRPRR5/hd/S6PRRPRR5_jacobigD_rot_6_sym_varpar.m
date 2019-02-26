% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR5_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:06:28
% EndTime: 2019-02-26 20:06:28
% DurationCPUTime: 0.04s
% Computational Cost: add. (22->10), mult. (84->29), div. (0->0), fcn. (88->8), ass. (0->19)
t158 = sin(pkin(6));
t161 = sin(qJ(3));
t171 = t158 * t161;
t162 = sin(qJ(2));
t170 = t158 * t162;
t160 = cos(pkin(6));
t169 = t160 * t162;
t164 = cos(qJ(2));
t168 = t160 * t164;
t167 = qJD(2) * t161;
t157 = sin(pkin(11));
t159 = cos(pkin(11));
t166 = t157 * t164 + t159 * t169;
t165 = -t157 * t169 + t159 * t164;
t163 = cos(qJ(3));
t156 = t158 * t164 * t167 + (t160 * t161 + t163 * t170) * qJD(3);
t155 = (t157 * t171 + t165 * t163) * qJD(3) + (-t157 * t168 - t159 * t162) * t167;
t154 = (-t159 * t171 + t166 * t163) * qJD(3) + (-t157 * t162 + t159 * t168) * t167;
t1 = [0, 0, t165 * qJD(2), 0, t155, t155; 0, 0, t166 * qJD(2), 0, t154, t154; 0, 0, qJD(2) * t170, 0, t156, t156;];
JgD_rot  = t1;
