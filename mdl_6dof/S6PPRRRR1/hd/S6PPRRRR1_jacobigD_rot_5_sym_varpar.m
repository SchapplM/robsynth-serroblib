% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRRR1_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:43
% EndTime: 2019-02-26 19:42:43
% DurationCPUTime: 0.09s
% Computational Cost: add. (20->10), mult. (80->28), div. (0->0), fcn. (88->10), ass. (0->17)
t165 = sin(pkin(12));
t171 = cos(pkin(6));
t176 = t165 * t171;
t166 = sin(pkin(7));
t167 = sin(pkin(6));
t175 = t167 * t166;
t169 = cos(pkin(12));
t174 = t169 * t171;
t173 = cos(qJ(3));
t172 = sin(qJ(3));
t170 = cos(pkin(7));
t168 = cos(pkin(13));
t164 = sin(pkin(13));
t163 = (t166 * t171 * t172 + (t168 * t170 * t172 + t164 * t173) * t167) * qJD(3);
t162 = ((-t164 * t176 + t169 * t168) * t173 + ((-t169 * t164 - t168 * t176) * t170 + t165 * t175) * t172) * qJD(3);
t161 = ((t164 * t174 + t165 * t168) * t173 + ((-t165 * t164 + t168 * t174) * t170 - t169 * t175) * t172) * qJD(3);
t1 = [0, 0, 0, t162, t162, 0; 0, 0, 0, t161, t161, 0; 0, 0, 0, t163, t163, 0;];
JgD_rot  = t1;
