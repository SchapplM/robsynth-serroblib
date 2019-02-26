% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR10_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_jacobigD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:18
% EndTime: 2019-02-26 21:59:18
% DurationCPUTime: 0.06s
% Computational Cost: add. (33->17), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->23)
t167 = sin(pkin(6));
t170 = sin(qJ(1));
t185 = t167 * t170;
t172 = cos(qJ(1));
t184 = t167 * t172;
t169 = sin(qJ(2));
t183 = t169 * t170;
t182 = t169 * t172;
t171 = cos(qJ(2));
t181 = t170 * t171;
t180 = t172 * t171;
t179 = qJD(1) * t167;
t166 = pkin(12) + qJ(4);
t164 = sin(t166);
t178 = qJD(2) * t164;
t177 = qJD(2) * t167;
t168 = cos(pkin(6));
t176 = t168 * t180 - t183;
t175 = t168 * t181 + t182;
t174 = t168 * t182 + t181;
t173 = -t168 * t183 + t180;
t165 = cos(t166);
t1 = [0, t172 * t179, 0, t176 * qJD(1) + t173 * qJD(2) (t164 * t185 + t173 * t165) * qJD(4) - t175 * t178 + (-t174 * t164 - t165 * t184) * qJD(1), 0; 0, t170 * t179, 0, t175 * qJD(1) + t174 * qJD(2) (-t164 * t184 + t174 * t165) * qJD(4) + t176 * t178 + (t173 * t164 - t165 * t185) * qJD(1), 0; 0, 0, 0, t169 * t177, t171 * t164 * t177 + (t165 * t167 * t169 + t164 * t168) * qJD(4), 0;];
JgD_rot  = t1;
