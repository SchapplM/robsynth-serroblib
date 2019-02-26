% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRP7_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_jacobigD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:39
% EndTime: 2019-02-26 22:12:39
% DurationCPUTime: 0.06s
% Computational Cost: add. (33->17), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->23)
t168 = sin(pkin(6));
t171 = sin(qJ(1));
t186 = t168 * t171;
t173 = cos(qJ(1));
t185 = t168 * t173;
t170 = sin(qJ(2));
t184 = t170 * t171;
t183 = t170 * t173;
t172 = cos(qJ(2));
t182 = t171 * t172;
t181 = t173 * t172;
t180 = qJD(1) * t168;
t167 = qJ(3) + pkin(11);
t165 = sin(t167);
t179 = qJD(2) * t165;
t178 = qJD(2) * t168;
t169 = cos(pkin(6));
t177 = t169 * t181 - t184;
t176 = t169 * t182 + t183;
t175 = t169 * t183 + t182;
t174 = -t169 * t184 + t181;
t166 = cos(t167);
t1 = [0, t173 * t180, t177 * qJD(1) + t174 * qJD(2), 0 (t165 * t186 + t174 * t166) * qJD(3) - t176 * t179 + (-t175 * t165 - t166 * t185) * qJD(1), 0; 0, t171 * t180, t176 * qJD(1) + t175 * qJD(2), 0 (-t165 * t185 + t175 * t166) * qJD(3) + t177 * t179 + (t174 * t165 - t166 * t186) * qJD(1), 0; 0, 0, t170 * t178, 0, t172 * t165 * t178 + (t166 * t168 * t170 + t165 * t169) * qJD(3), 0;];
JgD_rot  = t1;
