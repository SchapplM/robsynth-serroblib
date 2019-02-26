% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPPR4_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:00
% EndTime: 2019-02-26 20:00:00
% DurationCPUTime: 0.04s
% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
t167 = sin(pkin(6));
t170 = sin(qJ(3));
t180 = t167 * t170;
t171 = sin(qJ(2));
t179 = t167 * t171;
t169 = cos(pkin(6));
t178 = t169 * t171;
t173 = cos(qJ(2));
t177 = t169 * t173;
t176 = qJD(2) * t170;
t166 = sin(pkin(10));
t168 = cos(pkin(10));
t175 = t166 * t173 + t168 * t178;
t174 = -t166 * t178 + t168 * t173;
t172 = cos(qJ(3));
t1 = [0, 0, t174 * qJD(2), 0, 0 (-t166 * t180 - t174 * t172) * qJD(3) - (-t166 * t177 - t168 * t171) * t176; 0, 0, t175 * qJD(2), 0, 0 (t168 * t180 - t175 * t172) * qJD(3) - (-t166 * t171 + t168 * t177) * t176; 0, 0, qJD(2) * t179, 0, 0, -t167 * t173 * t176 + (-t169 * t170 - t172 * t179) * qJD(3);];
JgD_rot  = t1;
