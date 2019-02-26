% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPPR6_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:06:23
% EndTime: 2019-02-26 22:06:23
% DurationCPUTime: 0.06s
% Computational Cost: add. (33->17), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->23)
t164 = sin(pkin(6));
t167 = sin(qJ(1));
t182 = t164 * t167;
t169 = cos(qJ(1));
t181 = t164 * t169;
t166 = sin(qJ(2));
t180 = t166 * t167;
t179 = t166 * t169;
t168 = cos(qJ(2));
t178 = t167 * t168;
t177 = t169 * t168;
t176 = qJD(1) * t164;
t163 = qJ(3) + pkin(11);
t162 = cos(t163);
t175 = qJD(2) * t162;
t174 = qJD(2) * t164;
t165 = cos(pkin(6));
t173 = t165 * t177 - t180;
t172 = t165 * t178 + t179;
t171 = t165 * t179 + t178;
t170 = -t165 * t180 + t177;
t161 = sin(t163);
t1 = [0, t169 * t176, t173 * qJD(1) + t170 * qJD(2), 0, 0 (-t170 * t161 + t162 * t182) * qJD(3) - t172 * t175 + (t161 * t181 - t171 * t162) * qJD(1); 0, t167 * t176, t172 * qJD(1) + t171 * qJD(2), 0, 0 (-t171 * t161 - t162 * t181) * qJD(3) + t173 * t175 + (t161 * t182 + t170 * t162) * qJD(1); 0, 0, t166 * t174, 0, 0, t168 * t162 * t174 + (-t161 * t164 * t166 + t162 * t165) * qJD(3);];
JgD_rot  = t1;
